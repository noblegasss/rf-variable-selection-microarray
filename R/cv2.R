library(caret)

 Brain_10fold = createFolds(Brain_new$type, k = 10)
data = Brain_new
data$type = as.factor(data$type)
formula = type~.
data = Leukemia_new
fold = Leukemia_10fold

cv = function(formula, data, fold, k = 10, ntree = 500){
  formula = formula(formula)
  response = all.vars(formula)[1]
  class = unique(data[[response]])
  V1.sel = list()
  V2.sel = list()
  V3.sel = list()
  class.err1 = matrix(nrow = k, ncol = length(class))
  class.err2 = matrix(nrow = k, ncol = length(class))
  class.err3 = matrix(nrow = k, ncol = length(class))
  class.err4 = matrix(nrow = k, ncol = length(class))
  class.err5 = matrix(nrow = k, ncol = length(class))
  overall.err = matrix(nrow = k, ncol = 5)
  i = 0
  while(i < k){
    i = i+1
    training = data[-fold[[i]],]
    test = data[fold[[i]],]
    R1 = rfsrc(formula, data = training, ntree = ntree)
    pred1 = predict(R1, test)
    conf1 = caret::confusionMatrix(pred1$class, test[[response]])
    class.err1[i,] = 1 - diag(conf1$table)/ colSums(conf1$table)
    overall.err[i,1] = mean(pred1$class != test[[response]])
    V1 = VarSelImbal(formula, training, ntree = ntree, se = 0.001)
    R2 = V1$rf.select
    pred2 = predict.mulrf(R2, test)
    class.err2[i,] = pred2$class.error
    overall.err[i,2] = pred2$error
    V1.sel[[i]] = V1$variable.selected
    tran.sel1 = training[,colnames(training)%in%c(response,V1.sel[[i]])]
    R3 = rfsrc(formula, tran.sel1, ntree = ntree)
    pred3 = predict(R3, test)
    conf3 = caret::confusionMatrix(pred3$class, test[[response]])
    class.err3[i,] = 1 - diag(conf3$table)/ colSums(conf3$table)
    overall.err[i,3] = mean(pred3$class != test[[response]])
    V2 = var.select(formula, training, method = "vh", ntree = ntree, fast = TRUE)
    V2.sel[[i]] = V2$topvars
    R4 = V2$rfsrc.refit.obj
    pred4 = predict(R4, test)
    conf4 = caret::confusionMatrix(pred4$class, test[[response]])
    class.err4[i,] = 1 - diag(conf4$table)/ colSums(conf4$table)
    overall.err[i,4] = mean(pred4$class != test[[response]])
    Class = training[[response]]
    V3 = varSelRF(training, Class, ntree = ntree, c.sd = 0)
    V3.sel[[i]] = V3$selected.vars
    tran.sel3 = training[,colnames(training)%in%c(V3.sel[[i]])]
    R5 = rfsrc(formula, tran.sel3, ntree = ntree)
    pred5 = predict(R5, test)
    conf5 = caret::confusionMatrix(pred5$class, test[[response]])
    class.err5[i,] = 1 - diag(conf5$table)/ colSums(conf5$table)
    overall.err[i,5] = mean(pred5$class != test[[response]])
  }
  cv.err = colMeans(overall.err)
  cv.class.err1 = colMeans(class.err1)
  cv.class.err2 = colMeans(class.err2)
  cv.class.err3 = colMeans(class.err3)
  cv.class.err4 = colMeans(class.err4)
  cv.class.err5 = colMeans(class.err5)
  cv.class.err = rbind(cv.class.err1, cv.class.err2, cv.class.err3, cv.class.err4, cv.class.err5)
  colnames(cv.class.err) = class
  return(list(imbalanced.sel = V1.sel, RHVM.sel = V2.sel, VarSel.sel = V3.sel, overall.cv = cv.err, class.cv = cv.class.err))
}

cv2 = function(formula, data, fold, ntree = 500, i = i){
  formula = formula(formula)
  response = all.vars(formula)[1]
  class = unique(data[[response]])
  V1.sel = list()


 # class.err2 = matrix(nrow = k, ncol = length(class))

    training = data[-fold[[i]],]
    test = data[fold[[i]],]
    V1 = VarSelImbal(formula, training, ntree = ntree, se = 0.005)
    R2 = V1$rf.select
    pred2 = predict.mulrf(R2, test)
    class.err2 = pred2$class.error
    overall.err = pred2$error
    V1.sel = V1$variable.selected

  names(class.err2) = class
  return(list(imbalanced.sel = V1.sel, overall.cv = overall.err, class.cv = class.err2))
}


Brain_new$type = as.factor(Brain_new$type)
                                                                                                                                                                         
cvBrain = cv(formula = type~., data =  Brain_new, fold = Brain_10fold)

Leukemia_new = read.csv("Data/Leukemia_new.csv", row.names = 1, check.names = F)
data = Leukemia_new
Leukemia_new$type = as.factor(Leukemia_new$type)


Leukemia_10fold = createFolds(Leukemia$type, k = 10)
cvLeukemia = cv(formula = type~., data =  Leukemia_new, fold = Leukemia_10fold)

Leukemia_3fold = createFolds(Leukemia_new$type, k = 3)
Brain_3fold = createFolds(Brain_new$type, k = 3)

saveRDS(cvBrain1000,"Data/cvBrain1000.rds")
saveRDS(cvLeukemia1000, "DATA/cvLeukemia1000.rds")


cvBrain1000 = cv(formula = type~., data = Brain_new, fold = Brain_10fold, ntree = 1000)

cvLeukemia1000 = cv(formula = type~., data = Leukemia, fold = Leukemia_10fold, ntree = 1000)


cv.only = function(ntree = x,){
  formula = formula(type~.)
  data = Leukemia_new
  fold = Leukemia_3fold
  k = 3
  se = 0.001
  response = all.vars(formula)[1]
  class = unique(data[[response]])
  V1.sel = list()
  class.err1 = matrix(nrow = k, ncol = length(class))
  class.err2 = matrix(nrow = k, ncol = length(class))
  overall.err = matrix(nrow = k, ncol = 2)
  i = 0
  while(i < k){
    i = i+1
    training = data[-fold[[i]],]
    test = data[fold[[i]],]
    V1 = VarSelImbal(formula, training, ntree = ntree, se = se)
    R2 = V1$rf.select
    pred2 = predict.mulrf(R2, test)
    conf2 = caret::confusionMatrix(as.factor(pred2$pred.class), test[[response]])
    class.err1[i,] = 1 - diag(conf2$table)/ colSums(conf2$table)
    overall.err[i,1] = mean(pred2$pred.class != test[[response]])
    V1.sel[[i]] = V1$variable.selected
    tran.sel1 = training[,colnames(training)%in%c(response,V1.sel[[i]])]
    R3 = rfsrc(formula, tran.sel1, ntree = ntree)
    pred3 = predict(R3, test)
    conf3 = caret::confusionMatrix(pred3$class, test[[response]])
    class.err2[i,] = 1 - diag(conf3$table)/ colSums(conf3$table)
    overall.err[i,2] = mean(pred3$class != test[[response]])
  }
  cv.err = colMeans(overall.err)
  cv.class.err1 = colMeans(class.err1)
  cv.class.err2 = colMeans(class.err2)
  cv.class.err = rbind(cv.class.err1, cv.class.err2)
  colnames(cv.class.err) = class
  return(list(imbalanced.sel = V1.sel, overall.cv = cv.err, class.cv = cv.class.err))
}

j = 0
cv.Brain0.001 = NULL
for (i in c(500,1000,2000,5000)){
  j = j+1
  cv.Brain0.001[[j]] = cv.only(type~., data = Brain_new, fold = Brain_3fold, k = 3, ntree = i)
}


library(parallel)

n_core = detectCores() - 1
cl = makeCluster(n_core,setup_strategy = "sequential")

cv.Brain0.001 = parLapply(cl,c(500,1000,2000,5000), function(x) cv.only(ntree = x))

cv.Brain0.005 = parLapply(cl,c(500,1000,2000,5000), function(x) cv.only(ntree = x))

cv.Leukemia0.001 = parLapply(cl,c(500,1000,2000,5000), function(x) cv.only(ntree = x))

cv.Leukemia0.005 = parLapply(cl,c(500,1000,2000,5000), function(x) cv.only(ntree = x))



