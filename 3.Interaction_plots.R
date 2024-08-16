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

#-----------------------------------------
#            Loss A
#-----------------------------------------
#--Discovery---
loss_type <- "poisson_loss_lasso"
plot(refit_model_A, type = "interaction")
trt <-refit_model_A$trts
y <-refit_model_A$y
pi.x <- refit_model_A$pi.x
wts <- 1 / (pi.x * (trt[1]) + (1 - pi.x) * (trt[2]))
resp <-y*wts
recom<-refit_model_A$recommended.trts
TX<-refit_model_A$trt.received

train_cohort <-data.frame(resp,TX,recom)
anov <- aov(resp ~ TX * recom, data = train_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Discovery p-int : %f',pvalue)
rm(pvalue)


#--Testing----
subset <- Freq_A[Freq_A$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
test_cohort <- refit_external(refit_model_A,df_dana)
test_cohort <- test_cohort[[2]]
temp <- calc_effects_table(test_cohort)
resp <-test_cohort$y * ((test_cohort$wts)/2)
TX <-test_cohort$TX
recom <-test_cohort$recom

test_cohort <-data.frame(resp,TX,recom)
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col = "#ebebeb") #  gray background
grid(nx = NULL, ny = NULL,col = "gray", lwd = 2) # grid lines
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="Recom",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2,add=TRUE)

anov <- aov(resp ~ TX * recom, data = test_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Testing p-int : %f',pvalue)
rm(pvalue)

#-----------------------------------------
#            Loss B
#-----------------------------------------
#--Discovery---
loss_type <- "sq_loss_gam"
plot(refit_model_B, type = "interaction")
trt <-refit_model_B$trts
y <-refit_model_B$y
pi.x <- refit_model_B$pi.x
wts <- 1 / (pi.x * (trt[1]) + (1 - pi.x) * (trt[2]))
resp <-y*wts
recom<-refit_model_B$recommended.trts
TX<-refit_model_B$trt.received

train_cohort <-data.frame(resp,TX,recom)
anov <- aov(resp ~ TX * recom, data = train_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Discovery p-int : %f',pvalue)
rm(pvalue)

#--Testing----
subset <- Freq_B[Freq_B$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
test_cohort <- refit_external(refit_model_B,df_dana)
test_cohort <- test_cohort[[2]]
temp <- calc_effects_table(test_cohort)
resp <-test_cohort$y * ((test_cohort$wts)/2)
TX <-test_cohort$TX
recom <-test_cohort$recom

test_cohort <-data.frame(resp,TX,recom)
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col = "#ebebeb") #  gray background
grid(nx = NULL, ny = NULL,col = "gray", lwd = 2) # grid lines
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2,add=TRUE)

anov <- aov(resp ~ TX * recom, data = test_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Testing p-int : %f',pvalue)
rm(pvalue)


#-----------------------------------------
#            Loss C
#-----------------------------------------
#--Discovery---
loss_type <- "poisson_loss_gam"
plot(refit_model_C, type = "interaction")
trt <-refit_model_C$trts
y <-refit_model_C$y
pi.x <- refit_model_C$pi.x
wts <- 1 / (pi.x * (trt[1]) + (1 - pi.x) * (trt[2]))
resp <-y*wts
recom<-refit_model_C$recommended.trts
TX<-refit_model_C$trt.received

train_cohort <-data.frame(resp,TX,recom)
anov <- aov(resp ~ TX * recom, data = train_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Discovery p-int : %f',pvalue)
rm(pvalue)


#--Testing----
subset <- Freq_C[Freq_C$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
test_cohort <- refit_external(refit_model_C,df_dana)
test_cohort <- test_cohort[[2]]
temp <- calc_effects_table(test_cohort)
resp <-test_cohort$y * ((test_cohort$wts)/2)
TX <-test_cohort$TX
recom <-test_cohort$recom

test_cohort <-data.frame(resp,TX,recom)
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col = "#ebebeb") #  gray background
grid(nx = NULL, ny = NULL,col = "gray", lwd = 2) # grid lines
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2,add=TRUE)

anov <- aov(resp ~ TX * recom, data = test_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Testing p-int : %f',pvalue)
rm(pvalue)

#-----------------------------------------
#            Loss D
#-----------------------------------------
#--Discovery---
loss_type <- "logistic_loss_gam"
plot(refit_model_D, type = "interaction")
trt <-refit_model_D$trts
y <-refit_model_D$y
pi.x <- refit_model_D$pi.x
wts <- 1 / (pi.x * (trt[1]) + (1 - pi.x) * (trt[2]))
resp <-y*wts
recom<-refit_model_D$recommended.trts
TX<-refit_model_D$trt.received

train_cohort <-data.frame(resp,TX,recom)
anov <- aov(resp ~ TX * recom, data = train_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Discovery p-int : %f',pvalue)
rm(pvalue)

#--Testing----
subset <- Freq_D[Freq_D$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
test_cohort <- refit_external(refit_model_D,df_dana)
test_cohort <- test_cohort[[2]]
temp <- calc_effects_table(test_cohort)
resp <-test_cohort$y * ((test_cohort$wts)/2)
TX <-test_cohort$TX
recom <-test_cohort$recom

test_cohort <-data.frame(resp,TX,recom)
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col = "#ebebeb") #  gray background
grid(nx = NULL, ny = NULL,col = "gray", lwd = 2) # grid lines
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2,add=TRUE)

anov <- aov(resp ~ TX * recom, data = test_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Testing p-int : %f',pvalue)
rm(pvalue)

#-----------------------------------------
#            Loss E
#-----------------------------------------
#--Discovery---
loss_type <- "sq_loss_xgboost"
plot(refit_model_E, type = "interaction")
trt <-refit_model_E$trts
y <-refit_model_E$y
pi.x <- refit_model_E$pi.x
wts <- 1 / (pi.x * (trt[1]) + (1 - pi.x) * (trt[2]))
resp <-y*wts
recom<-refit_model_E$recommended.trts
TX<-refit_model_E$trt.received

train_cohort <-data.frame(resp,TX,recom)
anov <- aov(resp ~ TX * recom, data = train_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Discovery p-int : %f',pvalue)
rm(pvalue)


#--Testing----
subset <- Freq_E[Freq_E$Freq>=freq_num,]
rownames(subset) <-NULL
sel_varnames <- as.character(subset$features)
test_cohort <- refit_external(refit_model_E,df_dana)
test_cohort <- test_cohort[[2]]
temp <- calc_effects_table(test_cohort)
resp <-test_cohort$y * ((test_cohort$wts)/2)
TX <-test_cohort$TX
recom <-test_cohort$recom

test_cohort <-data.frame(resp,TX,recom)
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col = "#ebebeb") #  gray background
grid(nx = NULL, ny = NULL,col = "gray", lwd = 2) # grid lines
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2,add=TRUE)

anov <- aov(resp ~ TX * recom, data = test_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Testing p-int : %f',pvalue)
rm(pvalue)


#-----------------------------------------
#            Composite
#-----------------------------------------
#--Discovery---
att <- att_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_disc,'disc',x)
train_cohort <-att[[2]]
temp <- calc_effects_table(train_cohort)
resp <-train_cohort$y * ((train_cohort$wts)/2)
TX <-train_cohort$TX
recom <-train_cohort$recom

train_cohort <-data.frame(resp,TX,recom)
interaction.plot(x.factor = train_cohort$recom, trace.factor = train_cohort$TX,
                 response = train_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col = "#ebebeb") #  gray background
grid(nx = NULL, ny = NULL,col = "gray", lwd = 2) # grid lines
interaction.plot(x.factor = train_cohort$recom, trace.factor = train_cohort$TX,
                 response = train_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2,add=TRUE)

anov <- aov(resp ~ TX * recom, data = train_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Discovery p-int : %f',pvalue)
rm(pvalue)


#--Testing---
att <- att_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_dana,'test',x)
test_cohort <-att[[2]]
temp <- calc_effects_table(test_cohort)
resp <-test_cohort$y * ((test_cohort$wts)/2)
TX <-test_cohort$TX
recom <-test_cohort$recom

test_cohort <-data.frame(resp,TX,recom)
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col = "#ebebeb") #  gray background
grid(nx = NULL, ny = NULL,col = "gray", lwd = 2) # grid lines
interaction.plot(x.factor = test_cohort$recom, trace.factor = test_cohort$TX,
                 response = test_cohort$resp,trace.label="TX",col =  c("#F8766D","#00BFC4"),lty = 1,lwd = 2,add=TRUE)

anov <- aov(resp ~ TX * recom, data = test_cohort)
p<-summary(anov)
pvalue <- p[[1]]$`Pr(>F)`
pvalue <- pvalue[1:3]
pvalue<-stats::p.adjust(pvalue,'BH')
pvalue<-pvalue[3]
sprintf('Testing p-int : %f',pvalue)
rm(pvalue)