rm(list=ls(all=TRUE))
setwd("C:\\Users\\MBSaad\\Desktop\\Projects\\QASEM_NC\\Finalized_Qasem\\Modeling")

library(personalized)
library(rms)
library(dplyr)
library(tidyverse)
library(metaheuristicOpt)
library(cluster)
library(factoextra)
library(magrittr)

options(warn=-1)
set.seed(123)

source("fn_biomarkers\\fn_train_model.R")
source("fn_biomarkers\\fn_eval_model.R")
source("fn_biomarkers\\fn_ensemble.R")
source("fn_biomarkers\\summarize_subgroups.R")
source("fn_biomarkers\\fn_wes.R")



#------------------------------------------
#             Discovery
#-----------------------------------------
df_mda <- read.csv("Matched_MDA.csv")
df_mayo <- read.csv("Matched_MAYO.csv")
drops <- c('Liver.met','Brain.met','Met.status','OS','OS_events','PFS','PFS_events')
df_mda <- df_mda[ , !(names(df_mda) %in% drops)]
df_mayo <- df_mayo[ , !(names(df_mayo) %in% drops)]

df_natgen <-read.csv("Matched_NATGEN.csv")
drops <- c('OS','OS_events','PFS','PFS_events')
df_natgen <- df_natgen[ , !(names(df_natgen) %in% drops)]
df_disc <-rbind(df_mda,df_mayo,df_natgen)

df_dana <-read.csv("Matched_DANA.csv")
matched_ids <- read.csv("match_ids.csv")
#match_ids <- match_ids[1:277,]

x.varnames <- c("TP53",	"KRAS",	"EGFR",	"PIK3CA",	"CDKN2A",	"NF1",	"MET",	"ERBB2",	"BRAF",	"APC",	"ATM",	"BRCA2",	"KIT",	"PDGFRA",	"RB1",	"SMAD4",	"MYC",
                "NOTCH1",	"FGFR1",	"GNAS",	"PTEN",	"FBXW7","NTRK3","CCND1","CDK4",	"JAK2",	"ALK","PTPN11",
                "Gender",	"Age",	"Tobacco.Use",	"Pathology", "PD.L1.expression","Line.of.IO_conden")

estimator = "weighting"
style = 'all'
numFold <-5 #kfold CV
iter <-2 #repeated sampling number
freq_num <-3 # frequency a features being selected

#debug(cross_validation)
#-----------------------------------------
#            Loss A
#-----------------------------------------
dict_A <- list()
Freq_A <-data.frame()
train_cate_A <-data.frame()
valid_cate_A <-data.frame()
train_sample_A <-data.frame()
valid_sample_A <-data.frame()
for (run in 1:iter)
{
  print(sprintf("==== Sampling %d ====",run))
  loss_type <- "poisson_loss_lasso"
  matched_ids <-matched_ids[sample(nrow(matched_ids)), ]
  lossA <-cross_validation(df_disc,matched_ids,x.varnames,style,numFold,loss_type)
  dict_A[[run]] <- lossA
  sel_vars <- lossA[[6]]
  train_cate_A <-rbind(train_cate_A,dict_A[[run]][[1]])
  valid_cate_A <-rbind(valid_cate_A,dict_A[[run]][[2]])
  train_sample_A <-rbind(train_sample_A,dict_A[[run]][[7]])
  valid_sample_A <-rbind(valid_sample_A,dict_A[[run]][[8]])
  
  for (k in 1:(length(sel_vars)))
  {
    temp <- data.frame(sel_vars[[k]])
    colnames(temp) <- c('features')
    Freq_A <- rbind(Freq_A,temp)
    
  }
}

avg_iter_results(train_cate_A,valid_cate_A,train_sample_A,valid_sample_A)
  
Freq_A <- table(Freq_A)
Freq_A <- data.frame(Freq_A)
subset <- Freq_A[Freq_A$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
valid_set <-c()
valid_cate <-c()

refit_cate_A <- refit_features(df_disc,valid_set,sel_varnames,style,loss_type)
refit_model_A <- refit_cate_A[[2]]
disc_cate_A <- refit_cate_A[[1]]

#debug(avg_cate)
test_cate_A <-refit_external(refit_model_A,df_dana)
test_cate_A <-test_cate_A[[1]]
table <- avg_cate(disc_cate_A,valid_cate,test_cate_A,'Loss A')
knitr::kable(table, format = "markdown")


#debug(cross_validation)
#-----------------------------------------
#            Loss B
#-----------------------------------------
dict_B <- list()
Freq_B <-data.frame()
train_cate_B <-data.frame()
valid_cate_B <-data.frame()
train_sample_B <-data.frame()
valid_sample_B <-data.frame()

for (run in 1:iter)
{
  print(sprintf("==== Sampling %d ====",run))
  loss_type <- "sq_loss_gam"
  matched_ids <-matched_ids[sample(nrow(matched_ids)), ]
  lossB <-cross_validation(df_disc,matched_ids,x.varnames,style,numFold,loss_type)
  dict_B[[run]] <- lossB
  sel_vars <- lossB[[6]]
  train_cate_B <-rbind(train_cate_B,dict_B[[run]][[1]])
  valid_cate_B <-rbind(valid_cate_B,dict_B[[run]][[2]])
  train_sample_B <-rbind(train_sample_B,dict_B[[run]][[7]])
  valid_sample_B <-rbind(valid_sample_B,dict_B[[run]][[8]])
  
  for (k in 1:(length(sel_vars)))
  {
    temp <- data.frame(sel_vars[[k]])
    colnames(temp) <- c('features')
    Freq_B <- rbind(Freq_B,temp)
    
  }
}


avg_iter_results(train_cate_B,valid_cate_B,train_sample_B,valid_sample_B)

Freq_B <- table(Freq_B)
Freq_B <- data.frame(Freq_B)
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


#-----------------------------------------
#            Loss C
#-----------------------------------------
dict_C <- list()
Freq_C <-data.frame()
train_cate_C <-data.frame()
valid_cate_C <-data.frame()
train_sample_C <-data.frame()
valid_sample_C <-data.frame()

for (run in 1:iter)
{
  print(sprintf("==== Sampling %d ====",run))
  loss_type <- "poisson_loss_gam"
  matched_ids <-matched_ids[sample(nrow(matched_ids)), ]
  lossC <-cross_validation(df_disc,matched_ids,x.varnames,style,numFold,loss_type)
  dict_C[[run]] <- lossC
  sel_vars <- lossC[[6]]
  train_cate_C <-rbind(train_cate_C,dict_C[[run]][[1]])
  valid_cate_C <-rbind(valid_cate_C,dict_C[[run]][[2]])
  train_sample_C <-rbind(train_sample_C,dict_C[[run]][[7]])
  valid_sample_C <-rbind(valid_sample_C,dict_C[[run]][[8]])
  
  for (k in 1:(length(sel_vars)))
  {
    temp <- data.frame(sel_vars[[k]])
    colnames(temp) <- c('features')
    Freq_C <- rbind(Freq_C,temp)
    
  }
}

avg_iter_results(train_cate_C,valid_cate_C,train_sample_C,valid_sample_C)

Freq_C <- table(Freq_C)
Freq_C <- data.frame(Freq_C)
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


#-----------------------------------------
#            Loss D
#-----------------------------------------
dict_D <- list()
Freq_D <-data.frame()
train_cate_D <-data.frame()
valid_cate_D <-data.frame()
train_sample_D <-data.frame()
valid_sample_D <-data.frame()

for (run in 1:iter)
{
  print(sprintf("==== Sampling %d ====",run))
  loss_type <- "logistic_loss_gam"
  matched_ids <-matched_ids[sample(nrow(matched_ids)), ]
  lossD <-cross_validation(df_disc,matched_ids,x.varnames,style,numFold,loss_type)
  dict_D[[run]] <- lossD
  sel_vars <- lossD[[6]]
  train_cate_D <-rbind(train_cate_D,dict_D[[run]][[1]])
  valid_cate_D <-rbind(valid_cate_D,dict_D[[run]][[2]])
  train_sample_D <-rbind(train_sample_D,dict_D[[run]][[7]])
  valid_sample_D <-rbind(valid_sample_D,dict_D[[run]][[8]])
  
  for (k in 1:(length(sel_vars)))
  {
    temp <- data.frame(sel_vars[[k]])
    colnames(temp) <- c('features')
    Freq_D <- rbind(Freq_D,temp)
    
  }
}

avg_iter_results(train_cate_D,valid_cate_D,train_sample_D,valid_sample_D)

Freq_D <- table(Freq_D)
Freq_D <- data.frame(Freq_D)
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


#-----------------------------------------
#            Loss E
#-----------------------------------------
dict_E <- list()
Freq_E <-data.frame()
train_cate_E <-data.frame()
valid_cate_E <-data.frame()
train_sample_E <-data.frame()
valid_sample_E <-data.frame()

for (run in 1:iter)
{
  print(sprintf("==== Sampling %d ====",run))
  loss_type <- "sq_loss_xgboost"
  matched_ids <-matched_ids[sample(nrow(matched_ids)), ]
  lossE <-cross_validation(df_disc,matched_ids,x.varnames,style,numFold,loss_type)
  dict_E[[run]] <- lossE
  sel_vars <- lossE[[6]]
  train_cate_E <-rbind(train_cate_E,dict_E[[run]][[1]])
  valid_cate_E <-rbind(valid_cate_E,dict_E[[run]][[2]])
  train_sample_E <-rbind(train_sample_E,dict_E[[run]][[7]])
  valid_sample_E <-rbind(valid_sample_E,dict_E[[run]][[8]])
  
  for (k in 1:(length(sel_vars)))
  {
    temp <- data.frame(sel_vars[[k]])
    colnames(temp) <- c('features')
    Freq_E <- rbind(Freq_E,temp)
    
  }
}

avg_iter_results(train_cate_E,valid_cate_E,train_sample_E,valid_sample_E)

Freq_E <- table(Freq_E)
Freq_E <- data.frame(Freq_E)
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


#-----------------------------------------
#            Composite : Voting
#-----------------------------------------
T = 3
voting_train <- voting_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_disc,'disc',T)
voting_test <- voting_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_dana,'test',T)
table <- avg_cate(voting_train,valid_cate,voting_test,'VOTING')
knitr::kable(table, format = "markdown")


#-----------------------------------------
#            Composite : Averaging
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
results <- WeightedEnsemble(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_disc,'GA','disc')
x <-results[[1]]

#x <- data.frame(0.5,0.1,0.0,0.0,0.4)
#x <- data.frame(0.1,0.4,0.1,0.1,0.4)
#x <- data.frame(0.1,0.5,0.1,0.0,0.4)
#x <- t(x)
#colnames(x) <- "w"
#debug(avg_cate)
attention_train <- att_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_disc,'disc',x)
attention_train <-attention_train[[1]]
attention_test <- att_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_dana,'test',x)
attention_test <-attention_test[[1]]
table <- avg_cate(attention_train,valid_set,attention_test,'ATTENTION')
knitr::kable(table, format = "markdown")

