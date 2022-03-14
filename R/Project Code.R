library(ModelMetrics)
library(randomForest)
library(randomForestSRC)
library(genefilter)
library(hgu133plus2frmavecs)
library(hgu133plus2.db)
library(hgu133a.db)
library(hgu133afrmavecs)
library(annotate)
library(varSelRF)
library(dplyr)
library(tibble)
library(parallel)
library(doParallel)
library(foreach)

#### Load dataset

Brain = read.csv("Brain_GSE50161.csv",check.names=F,row.names = 1)
Leukemia = read.csv("Leukemia_GSE28497.csv",check.names=F,row.names = 1)


#### Filter Probe

exp = na.omit(t(Brain))[-1,]

arrayIQR<-apply(exp,1,IQR)
probe<-rownames(exp)
uniqueGenes<-findLargest(as.vector(probe),arrayIQR,"hgu133plus2.db")

exp2<-exp[uniqueGenes,]
geneSymbol<-getSYMBOL(rownames(exp2),"hgu133plus2.db")

rownames(exp2)<-geneSymbol
write.csv(exp2, "Data/Brain_new.csv")

exp = t(Leukemia[,-1])

arrayIQR<-apply(exp,1,IQR)
probe<-rownames(exp)
uniqueGenes<-findLargest(as.vector(probe),arrayIQR,"hgu133a.db")

exp2<-exp[uniqueGenes,]
geneSymbol<-getSYMBOL(rownames(exp2),"hgu133a.db")

rownames(exp2)<-geneSymbol

Leukemia_new = data.frame(type = Leukemia$type, t(exp2))
write.csv(Leukemia_new, "Data/Leukemia_new.csv")


#### Algorithm Code

##### Multi-RFQ

m_RFQ = function(formula,data, ntree = 1000){
  formula = formula(formula)
  response = all.vars(formula)[1]
  class = unique(data[[response]])
  pairs = t(combn(class,2))
  n = nrow(pairs)
  row_names = rownames(data)
  oob = matrix(0,nrow = nrow(data),ncol = length(class))
  colnames(oob) = class
  rownames(oob) = row_names
  imrf = list()
  i = 0
  while(i<n){
    i = i+1
    two_class = data %>% rownames_to_column("sample") %>% 
      filter(get(response) %in% pairs[i,]) %>% column_to_rownames("sample")
    two_class[,response] = as.factor(as.character(two_class[,response]))   
    obs_name = rownames(two_class)
    imrf[[i]] = imbalanced(formula, two_class, ntree = ntree)
    oob_class0 = data.frame(imrf[[i]]$class.oob)
    oob_class = model.matrix(~.-1, oob_class0)
    colnames(oob_class) = pairs[i,]
    rownames(oob_class) = obs_name
    oob_class = oob_class[match(row_names , rownames(oob_class)),]
    rownames(oob_class) = row_names
    oob_class[is.na(oob_class)] = 0
    oob[,colnames(oob) %in% colnames(oob_class)] = oob[,colnames(oob) %in% colnames(oob_class)] + oob_class
  }
  oob_final = colnames(oob)[apply(oob, 1, which.max)]
  names(oob_final) = row_names
  confM = caret::confusionMatrix(factor(data[[response]],levels = class),factor(oob_final,levels = class))
  table = confM$table
  oob_info = confM$overall
  oob_class.error = 1 - diag(table)/colSums(table, na.rm = T)
  max_oob.error = max(oob_class.error, na.rm = T)
  oob_error = mean(data[[response]] != oob_final)
  return(list(formula = formula,
              confusionMatrix = table, oob.info = oob_info, oob.class.error = oob_class.error,
              oob.error = oob_error,  max_oob.error = max_oob.error,
              rf = imrf, oob_prob = oob_final, 
              class = class, class.pairs = pairs))
}


predict.mRFQ = function(object, newdata){
  row_names = rownames(newdata)
  class = object$class
  newx = matrix(0,nrow = nrow(newdata),ncol = length(class))
  colnames(newx) = as.character(class)
  rownames(newx) = row_names
  response = all.vars(object$formula)[1]
  pairs = object$class.pairs
  n = nrow(pairs)
  imrf = object$rf
  i = 0
  while(i<n){
    i = i+1
    two_newx = newdata %>% rownames_to_column("sample") %>% 
      filter(get(response) %in% as.character(pairs[i,])) %>% column_to_rownames("sample")
    two_newx[,response] = as.factor(as.character(two_newx[,response]))   
    obs_name = rownames(two_newx)
    p = predict(imrf[[i]], two_newx)$class
    p = data.frame(p)
    pred = model.matrix(~.-1, p)
    colnames(pred) = pairs[i,]
    rownames(pred) = obs_name
    pred = pred[match(row_names , rownames(pred)),]
    rownames(pred) = row_names
    pred[is.na(pred)] = 0
    newx[,colnames(newx) %in% colnames(pred)] = newx[,colnames(newx) %in% colnames(pred)] + pred
  }
  pred_class = colnames(newx)[apply(newx, 1, which.max)]
  names(pred_class) = row_names
  confM = caret::confusionMatrix(factor(newdata[[response]],levels = class),factor(pred_class,levels = class))
  table = confM$table
  info = confM$overall
  class.error = 1 - diag(table)/colSums(table)
  error = unname(1-info[1])
  return(list(confusionMatrix = table, info = info,error = error, class.error = class.error,
              pred.class = pred_class, pred_vote = newx))
}

##### Variable selectin with multi-RFQ

RFQ_VSmax = function(formula, data , dropout = 0.2, ntree = 1000, se = 0.005){
  formula = formula(formula)
  response = all.vars(formula)[1]
  rf_all = rfsrc(type~., data = data, ntree = ntree, importance = TRUE)
  VIMP = data.frame(rf_all$importance)
  VIMP = VIMP %>% filter(all > 0) %>% arrange(desc(all))
  VIMP_name = rownames(VIMP)
  train_select = data[,colnames(data)%in%c(response,VIMP_name)]
  VIMP_new = list()
  rf_select = list()
  l = length(VIMP_name)
  i = 0
  oob_error = NULL
  while(l > 1){
    i = i+1
    VIMP_new[[i]] = VIMP_name
    rf_select[[i]] = m_RFQ(formula, data = train_select, ntree = ntree)
    oob_error[i] = rf_select[[i]]$max_oob.error
    l = floor(l * (1-dropout))
    VIMP_name = VIMP_name[1:l]
    train_select = train_select[,colnames(train_select)%in%c(response,VIMP_name)]
  }
  oob.min = which(oob_error <= min(oob_error, na.rm = T) + se)
  index = oob.min[length(oob.min)]
  var_select = VIMP_new[[index]]
  rf_final = rf_select[[index]]
  oob.info = rf_final$oob.info
  return(list(variable.selected = var_select, oob.error.calculate = oob_error,
              selection = index, oob.info = oob.info, rf.select = rf_final))
}


RFQ_VS = function(formula, data , dropout = 0.2, ntree = 1000, se = 0.005){
  formula = formula(formula)
  response = all.vars(formula)[1]
  rf_all = rfsrc(type~., data = data, ntree = ntree, importance = TRUE)
  VIMP = data.frame(rf_all$importance)
  VIMP = VIMP %>% filter(all > 0) %>% arrange(desc(all))
  VIMP_name = rownames(VIMP)
  train_select = data[,colnames(data)%in%c(response,VIMP_name)]
  VIMP_new = list()
  rf_select = list()
  l = length(VIMP_name)
  i = 0
  oob_error = NULL
  while(l > 1){
    i = i+1
    VIMP_new[[i]] = VIMP_name
    rf_select[[i]] = m_RFQ(formula, data = train_select, ntree = ntree)
    oob_error[i] = rf_select[[i]]$oob_error
    l = floor(l * (1-dropout))
    VIMP_name = VIMP_name[1:l]
    train_select = train_select[,colnames(train_select)%in%c(response,VIMP_name)]
  }
  oob.min = which(oob_error <= min(oob_error)+se)
  var_select = VIMP_new[[max(oob.min)]]
  rf_final = rf_select[[max(oob.min)]]
  oob.info = rf_final$oob.info
  return(list(variable.selected = var_select, oob.error.calculate = oob_error,
              selection = max(oob.min), oob.info = oob.info, rf.select = rf_final))
}


#### Ten-fold cv

cv = function(formula, data, fold, k = 10, ntree = 500){
  formula = formula(formula)
  response = all.vars(formula)[1]
  class = unique(data[[response]])
  V1.sel = list()
  V1.sel2 = list()
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
    V1 = RFQ_VS(formula, training, ntree = ntree, se = 0.005)
    R2 = V1$rf.select
    pred2 =predict.mRFQ(R2, test)
    class.err2[i,] = pred2$class.error
    overall.err[i,2] = pred2$error
    V1.sel[[i]] = V1$variable.selected
    V1.2 = RFQ_VSmax(formula, training, ntree = ntree, se = 0.005)
    R3 = v1.2$rf.select
    pred3 = predict.mRFQ(R3, test)
    class.err3[i,] = pred3$class.error
    overall.err[i,3] = pred3$error
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

#### Three-fold cv

cv3= function(formula = type~., data, fold, k = 3, ntree = 100, se = 0.005){
  formula = formula(formula)
  response = all.vars(formula)[1]
  var.sel = list()
  class = unique(data[[response]])
  overall.err = matrix(nrow = k, ncol = 5)
  class.err1 = matrix(nrow = k, ncol = length(class))
  class.err2 = matrix(nrow = k, ncol = length(class))
  class.err3 = matrix(nrow = k, ncol = length(class))
  class.err4 = matrix(nrow = k, ncol = length(class))
  class.err5 = matrix(nrow = k, ncol = length(class))
  i = 0
  while(i<k){
    i = i+1
    training = data[-fold[[i]],]
    test = data[fold[[i]],]
    rf_all = rfsrc(formula, training, ntree = ntree)
    mrfqmax = RFQ_VS.max(formula, training, ntree = ntree)
    mrfq = RFQ_VS(formula, training, ntree = ntree)
    var.sel[[i]] = mrfq$variable.selected
    var.sel2[[i]] = mrfqmax$variable.selected
    new_sel = training[,colnames(training) %in% c(response, var.sel[[i]])]
    rf_sel = rfsrc(formula, new_sel, ntree = ntree)
    new_sel = training[,colnames(training) %in% c(response, var.sel2[[i]])]
    rf_sel2 = rfsrc(formula, new_sel, ntree = ntree)
    mqrf_sel = mrfq$rf.select
    mqrf_sel2 = mrfqmax$rf.select
    p1 = predict(rf_all,test)
    p2 = predict.mRFQ(mqrf_sel2,test)
    p3 = predict.mRFQ(mqrf_sel,test)
    p4 = predict(rf_sel,test)
    p5 = predict(rf_sel2,test)
    confM = confusionMatrix(test[[response]], p1$class)
    table = confM$table
    class.err1[i,] = 1 - diag(table)/apply(table,2,sum)
    overall.err[i,1] = unname(1-confM$overall[1])
    class.err2[i,] = p2$class.error
    overall.err[i,2] = p2$error
    class.err3[i,] = p3$class.error
    overall.err[i,3] = p3$error
    confM = confusionMatrix(test[[response]], p4$class)
    table = confM$table
    class.err4[i,] = 1 - diag(table)/apply(table,2,sum)
    overall.err[i,4] = unname(1-confM$overall[1])
    confM = confusionMatrix(test[[response]], p5$class)
    table = confM$table
    class.err5[i,] = 1 - diag(table)/apply(table,2,sum)
    overall.err[i,5] = unname(1-confM$overall[1])
  }
  overall.err = colMeans(overall.err, na.rm = T)
  class.err1 = colMeans(class.err1, na.rm = T)
  class.err2 = colMeans(class.err2, na.rm = T)
  class.err3 = colMeans(class.err3, na.rm = T)
  class.err4 = colMeans(class.err4, na.rm = T)
  class.err5 = colMeans(class.err5, na.rm = T)
  class.err = rbind(class.err1,class.err2,class.err3,class.err4,class.err5)
  colnames(class.err) = class
  return(list(var_select = var.sel, overall_err = overall.err, class_err = class.err))
}

##### Create fold

Brain_10fold = createFolds(Brain_new$type, k = 10)
Leukemia_10fold = createFolds(Leukemia_new$type, k = 10)

Brain_3fold = createFolds(Brain_new$type, k = 3)
Leukemia_3fold = createFolds(Leukemia_new$type, k = 3)

Leukemia3 = foreach(i = c(100,500,1000,2000),.packages=c("dplyr", "tidyr", "randomForestSRC", "caret", "dplyr", "tibble")) %dopar% {
  cv3(type~.,data = Leukemia_new,fold = Leukemia_3fold,k = 3, ntree = i, se = 0.005)
}

Brain3 = foreach(i = c(100,500,1000,2000),.packages=c("dplyr", "tidyr", "randomForestSRC", "caret", "dplyr", "tibble")) %dopar% {
  cv3(type~.,data = Brain_new,fold = Brain_3fold,k = 3, ntree = i, se = 0.005)
}

cvBrain1000 = cv(formula = type~., data = Brain_new, fold = Brain_10fold, ntree = 1000)

cvLeukemia1000 = cv(formula = type~., data = Leukemia_new, fold = Leukemia_10fold, ntree = 1000)

