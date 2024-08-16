rm(list=ls(all=TRUE))
setwd("C:\\Users\\mbsaad\\Desktop\\Recovered Files\\Projects\\QASEM_NC\\Finalized_Qasem\\Modeling")
load('trained_models.RData')

library(personalized)
library(rms)
library(dplyr)
library(tidyverse)
library(metaheuristicOpt)
library(cluster)
library(factoextra)
library(magrittr)
library(repr)
library(ggplot2)

options(warn=-1)
set.seed(123)

source("fn_biomarkers\\fn_train_model.R")
source("fn_biomarkers\\fn_eval_model.R")
source("fn_biomarkers\\fn_ensemble.R")
source("fn_biomarkers\\summarize_subgroups.R")
source("fn_biomarkers\\fn_wes.R")

freq_num <-72 # frequency a features being selected


#-----------------------------------------
#            Loss A
#-----------------------------------------
loss_type <- "poisson_loss_lasso"
avg_iter_results(train_cate_A,valid_cate_A,train_sample_A,valid_sample_A)
subset <- Freq_A[Freq_A$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
valid_set <-c()
valid_cate <-c()

refit_cate_A <- refit_features(df_disc,valid_set,sel_varnames,style,loss_type)
refit_model_A <- refit_cate_A[[2]]
disc_cate_A <- refit_cate_A[[1]]
test_cate_A <-refit_external(refit_model_A,df_dana)
test_cate_A <-test_cate_A[[1]]
table <- avg_cate(disc_cate_A,valid_cate,test_cate_A,'Loss A')
knitr::kable(table, format = "markdown")

table <- Freq_A[order(Freq_A$Freq, decreasing = TRUE),]
rownames(table) <-NULL
table$features <- as.factor(table$features)

table %>%
  #mutate(name = fct_reorder(features, Freq)) %>%
  ggplot(aes(x=reorder(features, -desc(Freq)), y=Freq)) +
  geom_bar(stat="identity", fill="#076fa2", alpha=.7, width=.7) + 
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 60)) +
  coord_flip() +
  xlab("") +
  theme_bw()


#-----------------------------------------
#            Loss B
#-----------------------------------------
loss_type <- "sq_loss_gam"
avg_iter_results(train_cate_B,valid_cate_B,train_sample_B,valid_sample_B)

subset <- Freq_B[Freq_B$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
valid_set <-c()
valid_cate <-c()

refit_cate_B <- refit_features(df_disc,valid_set,sel_varnames,style,loss_type)
refit_model_B <- refit_cate_B[[2]]
disc_cate_B <- refit_cate_B[[1]]
test_cate_B <-refit_external(refit_model_B,df_dana)
test_cate_B <-test_cate_B[[1]]
table <- avg_cate(disc_cate_B,valid_cate,test_cate_B,'Loss B')
knitr::kable(table, format = "markdown")

table <- Freq_B[order(Freq_B$Freq, decreasing = TRUE),]
rownames(table) <-NULL
table$features <- as.factor(table$features)

table %>%
  #mutate(name = fct_reorder(features, Freq)) %>%
  ggplot(aes(x=reorder(features, -desc(Freq)), y=Freq)) +
  geom_bar(stat="identity", fill="#076fa2", alpha=.7, width=.7) + 
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 60)) +
  coord_flip() +
  xlab("") +
  theme_bw()


#-----------------------------------------
#            Loss C
#-----------------------------------------
loss_type <- "poisson_loss_gam"
avg_iter_results(train_cate_C,valid_cate_C,train_sample_C,valid_sample_C)

subset <- Freq_C[Freq_C$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
valid_set <-c()
valid_cate <-c()

refit_cate_C <- refit_features(df_disc,valid_set,sel_varnames,style,loss_type)
refit_model_C <- refit_cate_C[[2]]
disc_cate_C <- refit_cate_C[[1]]
test_cate_C <-refit_external(refit_model_C,df_dana)
test_cate_C <-test_cate_C[[1]]
table <- avg_cate(disc_cate_C,valid_cate,test_cate_C,'Loss C')
knitr::kable(table, format = "markdown")

table <- Freq_C[order(Freq_C$Freq, decreasing = TRUE),]
rownames(table) <-NULL
table$features <- as.factor(table$features)

table %>%
  #mutate(name = fct_reorder(features, Freq)) %>%
  ggplot(aes(x=reorder(features, -desc(Freq)), y=Freq)) +
  geom_bar(stat="identity", fill="#076fa2", alpha=.7, width=.7) + 
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 60)) +
  coord_flip() +
  xlab("") +
  theme_bw()




#-----------------------------------------
#            Loss D
#-----------------------------------------
loss_type <- "logistic_loss_gam"
avg_iter_results(train_cate_D,valid_cate_D,train_sample_D,valid_sample_D)

subset <- Freq_D[Freq_D$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
valid_set <-c()
valid_cate <-c()

refit_cate_D <- refit_features(df_disc,valid_set,sel_varnames,style,loss_type)
refit_model_D <- refit_cate_D[[2]]
disc_cate_D <- refit_cate_D[[1]]
test_cate_D <-refit_external(refit_model_D,df_dana)
test_cate_D <-test_cate_D[[1]]
table <- avg_cate(disc_cate_D,valid_cate,test_cate_D,'Loss D')
knitr::kable(table, format = "markdown")

table <- Freq_D[order(Freq_D$Freq, decreasing = TRUE),]
rownames(table) <-NULL
table$features <- as.factor(table$features)

table %>%
  #mutate(name = fct_reorder(features, Freq)) %>%
  ggplot(aes(x=reorder(features, -desc(Freq)), y=Freq)) +
  geom_bar(stat="identity", fill="#076fa2", alpha=.7, width=.7) + 
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 60)) +
  coord_flip() +
  xlab("") +
  theme_bw()



#-----------------------------------------
#            Loss E
#-----------------------------------------
loss_type <- "sq_loss_xgboost"
avg_iter_results(train_cate_E,valid_cate_E,train_sample_E,valid_sample_E)

subset <- Freq_E[Freq_E$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
valid_set <-c()
valid_cate <-c()

refit_cate_E <- refit_features(df_disc,valid_set,sel_varnames,style,loss_type)
refit_model_E <- refit_cate_E[[2]]
disc_cate_E <- refit_cate_E[[1]]
test_cate_E <-refit_external(refit_model_E,df_dana)
test_cate_E <-test_cate_E[[1]]
table <- avg_cate(disc_cate_E,valid_cate,test_cate_E,'Loss E')
knitr::kable(table, format = "markdown")


table <- Freq_E[order(Freq_E$Freq, decreasing = TRUE),]
rownames(table) <-NULL
table$features <- as.factor(table$features)

table %>%
  #mutate(name = fct_reorder(features, Freq)) %>%
  ggplot(aes(x=reorder(features, -desc(Freq)), y=Freq)) +
  geom_bar(stat="identity", fill="#076fa2", alpha=.7, width=.7) + 
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 60)) +
  coord_flip() +
  xlab("") +
  theme_bw()

#-----------------------------------------
#            Composite : hard Voting
#-----------------------------------------
T = 3
voting_train <- voting_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_disc,'disc',T)
voting_test <- voting_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_dana,'test',T)
table <- avg_cate(voting_train,valid_cate,voting_test,'VOTING')
knitr::kable(table, format = "markdown")


#-----------------------------------------
#            Composite : soft voting
#-----------------------------------------
#debug(avg_scheme)
avg_train <- avg_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_disc,'disc')
avg_test <- avg_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_dana,'test')

table <- avg_cate(avg_train,valid_set,avg_test,'AVERAGE')
knitr::kable(table, format = "markdown")


#-----------------------------------------
#            Composite : Attention
#-----------------------------------------

#debug(WeightedEnsemble)
#results <- WeightedEnsemble(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_disc,'GWO','disc')
x <-results[[1]]

#x <- data.frame(0.5,0.1,0.0,0.0,0.4)
#x <- data.frame(0.1,0.8,0.0,0.0,0.2)
#x <- data.frame(0.1,0.5,0.1,0.0,0.4)
#x <- t(x)
#colnames(x) <- "w"
#debug(avg_cate)
att <- att_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_disc,'disc',x)
attention_train <-att[[1]]

att <- att_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_dana,'test',x)
attention_test <-att[[1]]

table <- avg_cate(attention_train,valid_set,attention_test,'ATTENTION')
knitr::kable(table, format = "markdown")

