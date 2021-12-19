#import packages
library(randomForest)       
library(tidyverse)
library(gbm)         

#load training and gene info data
ALS_train = read_tsv("data/clean/ALS_train.tsv")
gene_info=read_tsv("data/clean/gene_info.tsv")


set.seed(1)
# tune the random forests with different values of m

#create different m values including sqrt(2500)
mvalues = sort(c(as.integer(seq.int(from = 1, to = 2500, length.out = 30)),sqrt(2500)))
oob_errors = numeric(length(mvalues))
ntree = 500
for(idx in 1:length(mvalues)){
  m = mvalues[idx]
  rf_fit = randomForest(factor(diagnosis) ~ ., mtry = m, data = ALS_train) 
  oob_errors[idx] = rf_fit$err.rate[ntree]
}

# plot OOB error versus m
m_and_oob_errors = tibble(m = mvalues, oob_err = oob_errors) 

p1=m_and_oob_errors %>%
  ggplot(aes(x = m, y = oob_err)) + 
  geom_line() + geom_point() + 
  scale_x_continuous(breaks = mvalues) +
  labs(x='M values',y='Out of Bag Error')+
  theme_bw()
  ggsave(filename = "results/rf_m_oob.png", 
         plot = p1a, 
         device = "png", 
         width = 13, 
         height = 4)  

# extract m corresponding to min value of OOB error
best_m = m_and_oob_errors %>% arrange(oob_errors) %>% head(1) %>% pull(m)

rf_fit_best = randomForest(factor(diagnosis) ~ ., mtry = best_m, data = ALS_train) 
# plot OOB error as a function of number of trees with best m

#create and save rf tree fit plot
p2=tibble(oob_error = rf_fit_best$err.rate[,"OOB"], trees = 1:1000) %>%
  ggplot(aes(x = trees, y = oob_error)) + geom_line() + 
  labs(x = "Number of trees", y = "OOB error") + theme_bw()
ggsave(filename = "results/rf_trees_oob.png", 
       plot = p2, 
       device = "png", 
       width = 9, 
       height = 4)


# run tuned random forest with variable importance
rf_fit_tuned = randomForest(factor(diagnosis) ~ ., mtry = best_m, ntree = 500,
                            importance = TRUE, data = ALS_train)

# save the random forest object
save(rf_fit_tuned, file = "results/rf_fit_tuned.Rda")



#create and save importance plot
Random_Forest_Variable_Importance=rf_fit_tuned
png(width = 8, 
    height = 7,
    res = 300,
    units = "in", 
    filename = "results/rf_varImpPlot.png")
varImpPlot(Random_Forest_Variable_Importance, n.var = 20)
dev.off()

#save features in  importance plot by  gini index
gini=as.tibble(varImpPlot(Random_Forest_Variable_Importance, n.var = 20),rownames=NA)%>%
  rownames_to_column('Feature ID')%>%
  arrange(desc(MeanDecreaseGini))%>%
  head(15)%>%
  inner_join(gene_info)
write_tsv(gini, file = "results/rf_gini.tsv")

#save features in  importance plot for accuracy 
acc=as.tibble(varImpPlot(Random_Forest_Variable_Importance, n.var = 20),rownames=NA)%>%
  rownames_to_column('Feature ID')%>%
  arrange(desc(MeanDecreaseAccuracy))%>%
  head(15)%>%
  inner_join(gene_info)
write_tsv(acc, file = "results/rf_acc.tsv")
