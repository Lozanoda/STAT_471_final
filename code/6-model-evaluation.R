# load libraries
library(glmnetUtils)
library(tidyverse)

# load test data and gene info
ALS_test = read_tsv("data/clean/ALS_test.tsv")

#function to calculate metrics for many thresholds
get_metrics = function(yhat, y){
  err = mean(yhat != y)  
  positive = sum(y)
  negative = length(y) - positive
  fpr = sum((yhat == 1) & (y == 0)) / negative
  fnr = sum((yhat == 0) & (y == 1)) / positive
  return(c(err, fpr, fnr))
}


#load log fit object
load("results/glm_fit.Rda")

#get logistic fit probabilities
log_prob= predict(glm_fit,                        
                         newdata = ALS_test, 
                         type = "response") %>%            
  as.numeric()


#calculate log missclassification rates for many thresholds
glm_err = c(0, 0, 0, 0)
glm_fpr = c(0, 0, 0, 0)
glm_fnr = c(0, 0, 0, 0)
thresholds = c(0.1, 0.3, 0.5, 0.7)
for(i in 1:4) {
  glm_predictions_binary = as.numeric(log_prob >= thresholds[i])
  curr_metrics = get_metrics(glm_predictions_binary, ALS_test$diagnosis)
  glm_err[i] = curr_metrics[1] 
  glm_fpr[i] = curr_metrics[2]
  glm_fnr[i] = curr_metrics[3]
}

#make and save tibble with glm metrics for thresholds
glm_threshs=tibble(Threshold = thresholds, 'Misclassification error' = glm_err, 'False positive rate' = glm_fpr, 'False Negative rate' = glm_fnr)
write_tsv(glm_threshs,"results/glm_threshs.tsv")


# load ridge fit object
load("results/ridge_fit.Rda")

# ridge predictions for test data
ridge_prob= predict(ridge_fit,                        
                  newdata = ALS_test,
                  s = "lambda.1se",
                  type = "response") %>%            
  as.numeric()


#calculate ridge missclassification rates for many thresholds
ridge_err = c(0, 0, 0, 0)
ridge_fpr = c(0, 0, 0, 0)
ridge_fnr = c(0, 0, 0, 0)
thresholds = c(0.1, 0.3, 0.5, 0.7)
for(i in 1:4) {
  ridge_predictions_binary = as.numeric(ridge_prob >= thresholds[i])
  curr_metrics = get_metrics(ridge_predictions_binary, ALS_test$diagnosis)
  ridge_err[i] = curr_metrics[1] 
  ridge_fpr[i] = curr_metrics[2]
  ridge_fnr[i] = curr_metrics[3]
}

#make and save tibble with ridge metrics for thresholds
ridge_threshs=tibble(Threshold = thresholds, 'Misclassification error' = ridge_err, 'False positive rate' = ridge_fpr, 'False negative rate' = ridge_fnr)
write_tsv(ridge_threshs,"results/ridge_threshs.tsv")

# load lasso fit object
load("results/lasso_fit.Rda")

# lasso predictions for test data
lasso_prob= predict(lasso_fit,                        
                    newdata = ALS_test,
                    s = "lambda.1se",
                    type = "response") %>%            
  as.numeric()


#calculate lasso missclassification rates for many thresholds
lasso_err = c(0, 0, 0, 0)
lasso_fpr = c(0, 0, 0, 0)
lasso_fnr = c(0, 0, 0, 0)
thresholds = c(0.1, 0.3, 0.5, 0.7)
for(i in 1:4) {
  lasso_predictions_binary = as.numeric(lasso_prob >= thresholds[i])
  curr_metrics = get_metrics(lasso_predictions_binary, ALS_test$diagnosis)
  lasso_err[i] = curr_metrics[1] 
  lasso_fpr[i] = curr_metrics[2]
  lasso_fnr[i] = curr_metrics[3]
}

#make and save tibble with lasso metrics for thresholds
lasso_threshs=tibble(Threshold = thresholds, 'Misclassification error' = lasso_err, 'False positive rate' = lasso_fpr, 'False negative rate' = lasso_fnr)
write_tsv(lasso_threshs,"results/lasso_threshs.tsv")

# load elnet fit object
load("results/elnet_fit_best.Rda")
load("results/elnet_fit.Rda")
elnet_predictions = predict(elnet_fit, 
                            alpha = elnet_fit_best$alpha,
                            newdata = ALS_test,
                            type = "response",
                            s = "lambda.1se") %>% as.numeric()

#calculate elnet missclassification rates for many thresholds
elnet_err = c(0, 0, 0, 0)
elnet_fpr = c(0, 0, 0, 0)
elnet_fnr = c(0, 0, 0, 0)
thresholds = c(0.1, 0.3, 0.5, 0.7)
for(i in 1:4) {
  elnet_predictions_binary = as.numeric(elnet_predictions >= thresholds[i])
  curr_metrics = get_metrics(elnet_predictions_binary, ALS_test$diagnosis)
  elnet_err[i] = curr_metrics[1] 
  elnet_fpr[i] = curr_metrics[2]
  elnet_fnr[i] = curr_metrics[3]
}

#make and save tibble with lelnet metrics for thresholds
elnet_threshs=tibble(Threshold = thresholds, 'Misclassification error' = elnet_err, 'False positive rate' = elnet_fpr, 'False negative rate' = elnet_fnr)
write_tsv(elnet_threshs,"results/elnet_threshs.tsv")



#load random forest model 
load("results/rf_fit_tuned.Rda")

#compute test misclassification error of the random forest
rf_predictions = predict(rf_fit_tuned,
type = "response", 
newdata = ALS_test)

rf_metrics = c(0.5,get_metrics(rf_predictions , ALS_test$diagnosis))


#collect best missclassifcation rates for all models
glm_threshs_best=glm_threshs%>%
  filter(Threshold==0.5)

lasso_threshs_best=lasso_threshs%>%
  slice_min(`Misclassification error`)

elnet_threshs_best=elnet_threshs%>%
  slice_min(`Misclassification error`)

ridge_threshs_best=ridge_threshs%>%
  slice_min(`Misclassification error`)

#makes sure all names are the same for joining
names(glm_threshs_best) <- names(lasso_threshs_best)
names(glm_threshs_best) <- names(elnet_threshs_best)
names(glm_threshs_best) <- names(ridge_threshs_best)

model_evaluation=rbind(glm_threshs_best,lasso_threshs_best,elnet_threshs_best,ridge_threshs_best)%>%
  add_row(Threshold=rf_metrics[1], 'Misclassification error' = rf_metrics[2], 'False positive rate' = rf_metrics[3] ,'False negative rate' = rf_metrics[4])%>%
  add_column(Model =c('Logistic','Lasso','Elastic Net','Ridge','Random Forest'), .before = 'Threshold')

write_tsv(model_evaluation,"results/model_evaluation.tsv")

top_genes_mul_models=gene_info%>%
  filter(`Gene Acronym` %in% c('ATP5I','AIF1','CDC37','YARS','GABBR1','DEGS1','MRPS18C','TRDMT1','MMP9','SLC20A1','PRKAR1A'))
  
write_tsv(top_genes_mul_models,"results/top_genes_mul_models.tsv")


