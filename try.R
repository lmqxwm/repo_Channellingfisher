library(haven)
library(fastDummies)
library(stringr)
library(VGAM)
library(nnet)

data = read_dta("AER-2008-1240_R1_AbelerFalkGoetteHuffman2009_data_file.dta")
data = dummy_cols(data, c("controls_temperature", "controls_time_of_day")) # make dummy vars
colnames(data)[(ncol(data)-5):ncol(data)] = c("CT1", "CT2", "CT3", "CTD1", "CTD2", "CTD3")
write_dta(data, "DataAFGH.dta")

quiet <- function(x) { # supress messages
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

data = read_dta("DataAFGH.dta")
regs = readLines("regressions.txt")
cmds = c() # commands of regressions
table_reg = c() # (ntable, regs in this table)
ntable = 0
n = 0
testvars = list()
regtypes = c()
for (i in 1:length(regs)){
  if (i == length(regs)){
    table_reg = rbind(table_reg, c(ntable, n)) # last table
  }
  if (sum(str_detect(regs[i], regex("table", ignore_case=TRUE)))>0){
    if (ntable>0){
      table_reg = rbind(table_reg, c(ntable, n))
    }
    ntable = ntable + 1
    n = 0
  }else{
    n = n+1
    sent = str_split_1(regs[i], "[:blank:]|\\|")
    sent = sent[!sent%in%c("")]
    q1 = which(str_detect(sent, "\\)"))[1]
    sent = str_split_1(regs[i], "\\(|\\)|[:blank:]|\\|")
    sent = sent[!sent%in%c("")]
    testvar = sent[2:q1] # num of test vars
    testvars= c(testvars, list(testvar))
    q2 = which(str_detect(sent, "if"))[1]
    if (sum(str_detect(sent, "\\,"))>0){
      q3 = which(str_detect(sent, "\\,"))[1]
    }else{ q3 = 0 }
    ind = paste("temp[which(", 
                paste("temp['",
                      sent[(q2+1):ifelse(q3>0, q3-1, length(sent))], 
                      "']==1", sep="", collapse="|"), 
                "),]", sep="") # if (|)
    regtype = sent[q1+1]
    regtypes = c(regtypes, regtype)
    yvar = sent[q1+2]
    xvar = sent[(q1+3):(q2-1)]
    if (regtype=="reg"){
      cmd = paste("model=lm(",
                  paste(yvar, "~", paste(xvar, collapse="+"), sep=""),
                  ",data=", ind, ")", sep="")
    } else if (regtype=="tobit") {
      upp = ""
      low = ""
      if(sum(str_detect(sent, "ul"))>0){
           upp = paste("Upper=max(temp[which(", 
                    paste("temp['",
                          sent[(q2+1):ifelse(q3>0, q3-1, length(sent))], 
                          "']==1", sep="", collapse="|"), 
                    "),'", yvar, "'])", sep="")
      }
      if(sum(str_detect(sent, "ll"))>0){
           low = paste("Lower=min(temp[which(", 
                    paste("temp['",
                          sent[(q2+1):ifelse(q3>0, q3-1, length(sent))], 
                          "']==1", sep="", collapse="|"), 
                    "),'", yvar, "'])", sep="")
      }
      cen = paste(upp, low, sep="")
      cmd = paste("model=vglm(",
                  paste(yvar, "~", paste(xvar, collapse="+"), sep=""),
                  ",family=tobit(", cen,
                  "),data=", ind, ")", sep="")
    }else if (regtype=="mlogit"){
      cmd = paste("model=multinom(",
                  paste(yvar, "~", paste(xvar, collapse="+"), sep=""),
                  ",data=", ind, ")", sep="")
    }else{
      stop("Wrong name!")
      }
    cmds = c(cmds, cmd)
  }
}


nreg = length(cmds)
ntestvar = 0 
# total number of testing variables
# (multiple logistic regression will have two placeholders with two non-baseline categories)
for (i in 1:nreg){
  if (regtypes[i] == "mlogit"){
    ntestvar = ntestvar + 2*length(testvars[[i]])
  }else{
    ntestvar = ntestvar + length(testvars[[i]])
  }
}



reps = 10 # repeating times

### randomization tests
resultF = matrix(nrow=reps, ncol=nreg*3)
resultB = matrix(nrow=reps, ncol=ntestvar*2)
for (t in 1:reps){
  cat("Repeat time: ", t, "\n")
  treatvars = c("treat_hi","treat_nosal","treat_r","treat_lo","treat_sal")
  treatinds = colnames(data) %in% treatvars
  N = nrow(data)
  seed_t = t
  set.seed(seed_t)
  rand = sample(1:N)
  temp = cbind(data[,treatinds], data[rand,!treatinds])
  now_var = 0 # test var now
  for (i in 1:nreg){
    #cat("Processing regression ", i, "\n")
    quiet(eval(parse(text=cmds[i])))
    model_wo = quiet(update(model, paste("~ .-",paste(testvars[[i]], collapse="-"))))
    if (regtypes[i] == "reg"){
      sum_ = summary(model)
      sum_wo = anova(model_wo, model)
      resultF[t, (3*i-2):(3*i)] = c(sum_wo$`Pr(>F)`[2], sum_wo$Df[2]-length(testvars[[i]]), sum_wo$Res.Df[1])
      for (j in 1:length(testvars[[i]])){
        resultB[t, (2*now_var+1):(2*now_var+2)] = sum_$coefficients[j+1, 1:2]
        now_var = now_var + 1
      }
    } else if (regtypes[i] == "tobit"){
      sum_ = summary(model)
      sum_wo = anova(model_wo, model, type="I")
      resultF[t, (3*i-2):(3*i)] = c(sum_wo$`Pr(>Chi)`[2], sum_wo$Df[2]-length(testvars[[i]]), sum_wo$`Resid. Df`[2])
      for (j in 1:length(testvars[[i]])){
        resultB[t, (2*now_var+1)] = sum_@coef3[j+2, 1]
        resultB[t, (2*now_var+2)] = sqrt(vcov(model)[j+2, j+2])
        now_var = now_var + 1
      }
    } else if (regtypes[i] == "mlogit"){
      while (dim(summary(model)$coefficients)[1]<2){
        seed_t = seed_t + reps
        set.seed(seed_t)
        rand = sample(1:N)
        temp = cbind(data[,treatinds], data[rand,!treatinds])
        quiet(eval(parse(text=cmds[i])))
        model_wo = quiet(update(model, paste("~ .-",paste(testvars[[i]], collapse="-"))))
      }
      sum_ = summary(model)
      sum_wo = anova(model_wo, model)
      resultF[t, (3*i-2):(3*i)] = c(sum_wo$`Pr(Chi)`[2], sum_wo$`   Df`[2]-length(testvars[[i]])*2, sum_wo$`Resid. df`[2])
      for (j in 1:length(testvars[[i]])){
        resultB[t, (2*now_var+1)] = sum_$coefficients[1, j+1]
        resultB[t, (2*now_var+2)] = sum_$standard.errors[1, j+1]
        resultB[t, (2*now_var+3)] = sum_$coefficients[2, j+1]
        resultB[t, (2*now_var+4)] = sum_$standard.errors[2, j+1]
        now_var = now_var + 2
      }
    }
  }
}

### bootstrap
# resultF = matrix(nrow=reps, ncol=nreg*3)
# resultB = matrix(nrow=reps, ncol=ntestvar*2)
# for (t in 1:reps){
#   cat("Repeat time: ", t, "\n")
#   N = nrow(data)
#   seed_t = t
#   set.seed(seed_t)
#   rand = sample(1:N, replace=TRUE)
#   temp = data[rand,]
#   now_var = 0 # test var now
#   for (i in 1:nreg){
#     #cat("Processing regression ", i, "\n")
#     quiet(eval(parse(text=cmds[i])))
#     model_wo = quiet(update(model, paste("~ .-",paste(testvars[[i]], collapse="-"))))
#     if (regtypes[i] == "reg"){
#       sum_ = summary(model)
#       sum_wo = anova(model_wo, model)
#       resultF[t, (3*i-2):(3*i)] = c(sum_wo$`Pr(>F)`[2], sum_wo$Df[2]-length(testvars[[i]]), sum_wo$Res.Df[1])
#       for (j in 1:length(testvars[[i]])){
#         resultB[t, (2*now_var+1):(2*now_var+2)] = sum_$coefficients[j+1, 1:2]
#         now_var = now_var + 1
#       }
#     } else if (regtypes[i] == "tobit"){
#       sum_ = summary(model)
#       sum_wo = anova(model_wo, model, type="I")
#       resultF[t, (3*i-2):(3*i)] = c(sum_wo$`Pr(>Chi)`[2], sum_wo$Df[2]-length(testvars[[i]]), sum_wo$`Resid. Df`[2])
#       for (j in 1:length(testvars[[i]])){
#         resultB[t, (2*now_var+1)] = sum_@coef3[j+2, 1]
#         resultB[t, (2*now_var+2)] = sqrt(vcov(model)[j+2, j+2])
#         now_var = now_var + 1
#       }
#     } else if (regtypes[i] == "mlogit"){
#       while (dim(summary(model)$coefficients)[1]<2){
#         seed_t = seed_t + reps
#         rand = sample(1:N, replace=TRUE)
#         temp = data[rand,]
#         quiet(eval(parse(text=cmds[i])))
#         model_wo = quiet(update(model, paste("~ .-",paste(testvars[[i]], collapse="-"))))
#       }
#       sum_ = summary(model)
#       sum_wo = anova(model_wo, model)
#       resultF[t, (3*i-2):(3*i)] = c(sum_wo$`Pr(Chi)`[2], sum_wo$`   Df`[2]-length(testvars[[i]])*2, sum_wo$`Resid. df`[2])
#       for (j in 1:length(testvars[[i]])){
#         resultB[t, (2*now_var+1)] = sum_$coefficients[1, j+1]
#         resultB[t, (2*now_var+2)] = sum_$standard.errors[1, j+1]
#         resultB[t, (2*now_var+3)] = sum_$coefficients[2, j+1]
#         resultB[t, (2*now_var+4)] = sum_$standard.errors[2, j+1]
#         now_var = now_var + 2
#       }
#     }
#   }
# }

### jackknife
# N = nrow(data)
# resultF = matrix(nrow=N, ncol=nreg*3)
# resultB = matrix(nrow=N, ncol=ntestvar*2)
# for (t in 1:N){
#   cat("Repeat time: ", t, "\n")
#   rand = setdiff(1:N, t)
#   temp = data[rand,]
#   now_var = 0 # test var now
#   for (i in 1:nreg){
#     #cat("Processing regression ", i, "\n")
#     quiet(eval(parse(text=cmds[i])))
#     model_wo = quiet(update(model, paste("~ .-",paste(testvars[[i]], collapse="-"))))
#     if (regtypes[i] == "reg"){
#       sum_ = summary(model)
#       sum_wo = anova(model_wo, model)
#       resultF[t, (3*i-2):(3*i)] = c(sum_wo$`Pr(>F)`[2], sum_wo$Df[2]-length(testvars[[i]]), sum_wo$Res.Df[1])
#       for (j in 1:length(testvars[[i]])){
#         resultB[t, (2*now_var+1):(2*now_var+2)] = sum_$coefficients[j+1, 1:2]
#         now_var = now_var + 1
#       }
#     } else if (regtypes[i] == "tobit"){
#       sum_ = summary(model)
#       sum_wo = anova(model_wo, model, type="I")
#       resultF[t, (3*i-2):(3*i)] = c(sum_wo$`Pr(>Chi)`[2], sum_wo$Df[2]-length(testvars[[i]]), sum_wo$`Resid. Df`[2])
#       for (j in 1:length(testvars[[i]])){
#         resultB[t, (2*now_var+1)] = sum_@coef3[j+2, 1]
#         resultB[t, (2*now_var+2)] = sqrt(vcov(model)[j+2, j+2])
#         now_var = now_var + 1
#       }
#     } else if (regtypes[i] == "mlogit"){
#       sum_ = summary(model)
#       sum_wo = anova(model_wo, model)
#       resultF[t, (3*i-2):(3*i)] = c(sum_wo$`Pr(Chi)`[2], sum_wo$`   Df`[2]-length(testvars[[i]])*2, sum_wo$`Resid. df`[2])
#       for (j in 1:length(testvars[[i]])){
#         resultB[t, (2*now_var+1)] = sum_$coefficients[1, j+1]
#         resultB[t, (2*now_var+2)] = sum_$standard.errors[1, j+1]
#         resultB[t, (2*now_var+3)] = sum_$coefficients[2, j+1]
#         resultB[t, (2*now_var+4)] = sum_$standard.errors[2, j+1]
#         now_var = now_var + 2
#       }
#     }
#   }
# }

# original regression results
resultFF = matrix(nrow=1, ncol=nreg*3)
resultBB = matrix(nrow=1, ncol=ntestvar*2)
temp = data
now_var = 0 # test var now
for (i in 1:nreg){
  #cat("Processing regression ", i, "\n")
  quiet(eval(parse(text=cmds[i])))
  model_wo = quiet(update(model, paste("~ .-",paste(testvars[[i]], collapse="-"))))
  if (regtypes[i] == "reg"){
    sum_ = summary(model)
    sum_wo = anova(model_wo, model)
    resultFF[1, (3*i-2):(3*i)] = c(sum_wo$`Pr(>F)`[2], sum_wo$Df[2]-length(testvars[[i]]), sum_wo$Res.Df[1])
    for (j in 1:length(testvars[[i]])){
      resultBB[1, (2*now_var+1):(2*now_var+2)] = sum_$coefficients[j+1, 1:2]
      now_var = now_var + 1
    }
  } else if (regtypes[i] == "tobit"){
    sum_ = summary(model)
    sum_wo = anova(model_wo, model, type="I")
    resultFF[1, (3*i-2):(3*i)] = c(sum_wo$`Pr(>Chi)`[2], sum_wo$Df[2]-length(testvars[[i]]), sum_wo$`Resid. Df`[2])
    for (j in 1:length(testvars[[i]])){
      resultBB[1, (2*now_var+1)] = sum_@coef3[j+2, 1]
      resultBB[1, (2*now_var+2)] = sqrt(vcov(model)[j+2, j+2])
      now_var = now_var + 1
    }
  } else if (regtypes[i] == "mlogit"){
    while (dim(summary(model)$coefficients)[1]<2){
      seed_t = seed_t + reps
      set.seed(seed_t)
      rand = sample(1:N)
      temp = cbind(data[,treatinds], data[rand,!treatinds])
      quiet(eval(parse(text=cmds[i])))
      model_wo = quiet(update(model, paste("~ .-",paste(testvars[[i]], collapse="-"))))
    }
    sum_ = summary(model)
    sum_wo = anova(model_wo, model)
    resultFF[1, (3*i-2):(3*i)] = c(sum_wo$`Pr(Chi)`[2], sum_wo$`   Df`[2]-length(testvars[[i]])*2, sum_wo$`Resid. df`[2])
    for (j in 1:length(testvars[[i]])){
      resultBB[1, (2*now_var+1)] = sum_$coefficients[1, j+1]
      resultBB[1, (2*now_var+2)] = sum_$standard.errors[1, j+1]
      resultBB[1, (2*now_var+3)] = sum_$coefficients[2, j+1]
      resultBB[1, (2*now_var+4)] = sum_$standard.errors[2, j+1]
      now_var = now_var + 2
    }
  }
}


# randomization-t
randt = matrix(nrow=reps, ncol=ntestvar)
for (t in 1:reps){
  for (i in 1:ntestvar){
    randt[t, i] = abs(resultB[t, (2*i-1)]/resultB[t, (2*i)])>abs(resultBB[1, (2*i-1)]/resultBB[1, (2*i)])
  }
}
t_M = matrix(nrow=2, ncol=ntestvar)
for (t in 1:reps){
  t_M[1,] = colMeans(randt) # lower bound
  t_M[2,] = colMeans(rbind(randt, rep(1, ntestvar))) # upper bound
}

# randomization-c
randc = matrix(nrow=reps, ncol=ntestvar)
for (t in 1:reps){
  for (i in 1:ntestvar){
    randc[t, i] = abs(resultB[t, (2*i-1)])>abs(resultBB[1, (2*i-1)])
  }
}
c_M = matrix(nrow=2, ncol=ntestvar)
for (t in 1:reps){
  c_M[1,] = colMeans(randc) # lower bound
  c_M[2,] = colMeans(rbind(randc, rep(1, ntestvar))) # upper bound
}
