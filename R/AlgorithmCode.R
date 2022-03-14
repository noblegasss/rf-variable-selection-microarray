library(ModelMetrics)
library(randomForest)
library(randomForestSRC)
library(varSelRF)
library(dplyr)
library(tibble)

Brain_new = read.csv("Data/Brain_new.csv", row.names = 1, check.names = F)

set.seed(510)
sample = sample.int(nrow(Brain_new), 110)
training = Brain_new[sample,]
y_train = Brain_new$type
test = Brain_new[-sample,]

training$type = as.factor(training$type)
data = training

## Multiclass imbalanced randomForest

multi_imrf = function(formula,data, ntree = 1000){
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

rfs = rfsrc(type~., data = data)
rfs

mulimrf = multi_imrf(type~., data = data, ntree = 1000)
object = mulimrf
newdata = test

im = imbalanced(type~., data = data, ntree = 100)

predict.mulrf = function(object, newdata){
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

newdata$type = as.factor(newdata$type)
pred = predict(rfs, newdata)
pred$class
caret::confusionMatrix(newdata$type, pred$class)

pred.im = predict.mulrf(mulimrf, newdata)
pred.im$pred_prob

data = training

## Variable selectin with multiclass imbalanced randomForest

VarSelImbal.max = function(formula, data , dropout = 0.2, ntree = 1000, se = 0.005){
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
    rf_select[[i]] = multi_imrf(formula, data = train_select, ntree = ntree)
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


VarSelImbal = function(formula, data , dropout = 0.2, ntree = 1000, se = 0.005){
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
    rf_select[[i]] = multi_imrf(formula, data = train_select, ntree = ntree)
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


set.seed(510)
sample = sample.int(nrow(Leukemia), 110)
training = Leukemia[sample,]
y_train = Leukemia$type
test = Leukemia[-sample,]


Var_sel = VarSelImbal.max(type~., data = data, dropout = 0.2, ntree = 100, se = 0.01)

train_s = training[,colnames(training)%in%c("type",Var_sel$variable.selected)]
rf_s = rfsrc(type~., train_s, ntree = 1000)
rf_s

rf_s2 = multi_imrf(type~., train_s2)
rf_s2$oob.info
rf_s2$class.error
rf_s2$avg.class.oob.error

test$type = as.factor(test$type)

pred = predict(rf_s, test)
pred

pred2 = predict.mulrf(Var_sel$rf.select, test)
pred2$oob.info[,11]

Var_sel2 = varSelRF(training[,-1], training$type, ntree = 500)
train_s2 = training[,colnames(training)%in%c("type",Var_sel2$selected.vars)]

rf_ss = rfsrc(type~., train_s2, ntree = 1000)
rf_ss
