 # load libraries
library(tidyverse)
library(glmnetUtils)                    # to run ridge and lasso
source("code/functions/plot_glmnet.R")            # for lasso/ridge trace plots

# read in the training data and gene info
ALS_train = read_tsv("data/clean/ALS_train.tsv")
gene_info=read_tsv("data/clean/gene_info.tsv")


#run logistic regression
set.seed(1)
glm_fit= glm(diagnosis ~. , family = "binomial", data = ALS_train,control = list(maxit = 50)) 


# save the  log fit object
save(glm_fit, file = "results/glm_fit.Rda")

#get top 10 negative coeffs
glm_coeffs_down<-as.tibble(summary(glm_fit)$coefficients)%>%
  mutate("Feature ID" = row.names(summary(glm_fit)$coefficients))%>%
  arrange(Estimate)%>%
  filter(`Feature ID`!= "(Intercept)")%>%
  head(10)%>%
  select(c("Estimate",`Feature ID`))%>%
  inner_join(gene_info)



#get top 10 positive coeffs
glm_coeffs_up<-as.tibble(summary(glm_fit)$coefficients)%>%
  mutate("Feature ID" = row.names(summary(glm_fit)$coefficients))%>%
  arrange(desc(Estimate))%>%
  filter(`Feature ID`!= "(Intercept)")%>%
  head(10)%>%
  select(c("Estimate",`Feature ID`))%>%
  inner_join(gene_info)

glm_coeffs_20=glm_coeffs_up%>%
  rbind(glm_coeffs_down)

#save table
write_tsv(glm_coeffs_20, file = "results/glm_coeffs_20.tsv")


# run ridge
set.seed(1)
ridge_fit = cv.glmnet(diagnosis ~ .,   
                      alpha = 0,
                      family = "binomial",
                      type.measure = "class",        
                      nfolds = 10,               
                      data = ALS_train)

# save the lasso fit object
save(ridge_fit, file = "results/ridge_fit.Rda") 

#save ridge feature plot 
p1 = plot_glmnet(ridge_fit, ALS_train, features_to_plot = 10)
ggsave(filename = "results/ridge-trace-plot.png", 
       plot = p1, 
       device = "png", 
       width = 9, 
       height = 4)

# save ridge CV plot
png(width = 6, 
    height = 4,
    res = 300,
    units = "in", 
    filename = "results/ridge-cv-plot.png")
plot(ridge_fit)
dev.off()

ridge_features=tibble('Feature ID'=p1[["data"]][["Feature"]][1:10])%>%
  inner_join(gene_info)
write_tsv(ridge_features,"results/ridge_features.tsv")


# run lasso 
set.seed(1)
lasso_fit = cv.glmnet(diagnosis ~ .,   
                      alpha = 1,
                      family = "binomial",
                      type.measure = "class",        
                      nfolds = 10,               
                      data = ALS_train)

# save the lasso fit object
save(lasso_fit, file = "results/lasso_fit.Rda")

# create lasso CV plot
png(width = 6, 
    height = 4,
    res = 300,
    units = "in", 
    filename = "results/lasso-cv-plot.png")
plot(lasso_fit)
dev.off()

# create lasso trace plot
p2 = plot_glmnet(lasso_fit, ALS_train, features_to_plot = 10)
ggsave(filename = "results/lasso-trace-plot.png", 
       plot = p2, 
       device = "png", 
       width = 9, 
       height = 4)


# find features in trace plot
lasso_trace_features=tibble('Feature ID'=p2[["data"]][["Feature"]][1:10])%>%
  inner_join(gene_info)
write_tsv(lasso_trace_features,"results/lasso_trace_features.tsv")

# extract features selected by lasso and their coefficients and related cellular data
beta_hat_std = extract_std_coefs(lasso_fit, ALS_train)
beta_hat_std =beta_hat_std %>%
  filter(coefficient != 0) %>%
  rename('Feature ID'=feature)%>%
  inner_join(gene_info)%>%
  arrange(desc(abs(coefficient)))%>%
  head(20)
write_tsv(beta_hat_std,"results/lasso-features-table.tsv")


# run elastic net 
set.seed(1)
elnet_fit = cva.glmnet(diagnosis ~ .,   
                      family = "binomial",
                      nfolds = 10,  
                      type.measure = "class",        
                      data = ALS_train)
#get best fit
elnet_fit_best = extract_best_elnet(elnet_fit)


# save the elnet fit and best fit object
save(elnet_fit_best, file = "results/elnet_fit_best.Rda")
save(elnet_fit, file = "results/elnet_fit.Rda")
# create  elnet CV plot
png(width = 6, 
    height = 4,
    res = 300,
    units = "in", 
    filename = "results/elnet-cv-plot.png")
plot(elnet_fit_best)
dev.off()

# create lasso trace plot
p3= plot_glmnet(elnet_fit_best, ALS_train, features_to_plot = 10)
ggsave(filename = "results/elnet-trace-plot.png", 
       plot = p3, 
       device = "png", 
       width = 9, 
       height = 4)

# find features in trace plot
elnet_trace_features=tibble('Feature ID'=p3[["data"]][["Feature"]][1:10])%>%
  inner_join(gene_info)
write_tsv(ridge_features,"results/elnet_trace_features.tsv")


