rm(list=ls(all=TRUE))
setwd("C:\\Users\\mbsaad\\Desktop\\Recovered Files\\Projects\\QASEM_NC\\Finalized_Qasem\\Modeling")
load('trained_models.RData')
#load('GWO_iter_30_tp53.RData')

library("survival")
library("ggplot2")
library("survminer")

source("fn_biomarkers\\fn_recom_plots.R")
source("fn_biomarkers\\fn_train_model.R")
source("fn_biomarkers\\fn_eval_model.R")
source("fn_biomarkers\\fn_ensemble.R")
source("fn_biomarkers\\summarize_subgroups.R")
source("fn_biomarkers\\fn_wes.R")


#-----------------------------------------
#            Prep Survival data
#-----------------------------------------
outcome <- c('OS','OS_events','PFS','PFS_events')

df_temp <- read.csv("Matched_MDA.csv")
mda_outcome <-df_temp[outcome]
df_temp <- read.csv("Matched_MAYO.csv")
mayo_outcome <-df_temp[outcome]
df_temp <-read.csv("Matched_NATGEN.csv")
natgen_outcome <-df_temp[outcome]
df_temp <-read.csv("Matched_DANA.csv")
dana_outcome <-df_temp[outcome]
rm(df_temp)

disc_outcome <-rbind(mda_outcome,mayo_outcome,natgen_outcome)

#-----------------------------------------
#            Censoring
# #-----------------------------------------
t <-24
for (i in 1:dim(disc_outcome)[1])
{
  time <-disc_outcome[i,]['PFS']
  event <-disc_outcome[i,]['PFS_events']
  if (time > t)
  {
    disc_outcome[i,]['PFS'] <- t
    disc_outcome[i,]['PFS_events'] <-0

  }
}


for (i in 1:dim(dana_outcome)[1])
{
  time <-dana_outcome[i,]['PFS']
  event <-dana_outcome[i,]['PFS_events']
  if (time > t)
  {
    dana_outcome[i,]['PFS'] <- t
    dana_outcome[i,]['PFS_events'] <-0

  }
}


#---Custom theme for plotting------------------------------
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(size = 14, color = "black",hjust=0.5,face = "bold"),
      axis.text.x = element_text(size = 14, color = "black", face = "bold"),
      legend.text = element_text(size = 14, color = "black", face = "bold"),
      legend.title = element_text(size = 14, color = "black", face = "bold"),
      axis.text.y = element_text(size = 14, color = "black", face = "bold"),
      axis.title.x = element_text(size = 14, color = "black", face = "bold"),
      axis.title.y = element_text(size = 14, color = "black", face = "bold") , #angle=(90))
    )
}


#-----------------------------------------
#            Loss A
#-----------------------------------------
recom <- refit_model_A$recommended.trts
actual <- df_disc$treatment
pid <- df_disc$patient_id
recomA <- cbind(pid,recom,actual)
recomA <- data.frame(recomA,disc_outcome)
fit <-recom_plot(recomA,custom_theme)
#debug(real_plot)
#fit <-real_plot(recomA,custom_theme)

test_cate_A <-refit_external(refit_model_A,df_dana)
recomA <-test_cate_A[[2]]
recomA <- recomA[c("TX","recom")]
colnames(recomA)[1]<-'actual'
recomA <- cbind(recomA,dana_outcome)
fit <-recom_plot(recomA,custom_theme)
#fit <-real_plot(recomA,custom_theme)

#-----------------------------------------
#            Loss B
#-----------------------------------------
recom <- refit_model_B$recommended.trts
actual <- df_disc$treatment
pid <- df_disc$patient_id
recomB <- cbind(pid,recom,actual)
recomB <- data.frame(recomB,disc_outcome)
fit<-recom_plot(recomB,custom_theme)
#fit<-real_plot(recomB,custom_theme)


test_cate_B <-refit_external(refit_model_B,df_dana)
recomB <-test_cate_B[[2]]
recomB <- recomB[c("TX","recom")]
colnames(recomB)[1]<-'actual'
recomB <- cbind(recomB,dana_outcome)
fit <-recom_plot(recomB,custom_theme)
#fit <-real_plot(recomB,custom_theme)


#-----------------------------------------
#            Loss C
#-----------------------------------------
recom <- refit_model_C$recommended.trts
actual <- df_disc$treatment
pid <- df_disc$patient_id
recomC <- cbind(pid,recom,actual)
recomC <- data.frame(recomC,disc_outcome)
fit<-recom_plot(recomC,custom_theme)
#fit<-real_plot(recomC,custom_theme)


test_cate_C <-refit_external(refit_model_C,df_dana)
recomC <-test_cate_C[[2]]
recomC <- recomC[c("TX","recom")]
colnames(recomC)[1]<-'actual'
recomC <- cbind(recomC,dana_outcome)
fit <-recom_plot(recomC,custom_theme)
#fit <-real_plot(recomC,custom_theme)


#-----------------------------------------
#            Loss D
#-----------------------------------------
recom <- refit_model_D$recommended.trts
actual <- df_disc$treatment
pid <- df_disc$patient_id
recomD <- cbind(pid,recom,actual)
recomD <- data.frame(recomD,disc_outcome)
fit<-recom_plot(recomD,custom_theme)
#fit<-real_plot(recomD,custom_theme)


test_cate_D <-refit_external(refit_model_D,df_dana)
recomD <-test_cate_D[[2]]
recomD <- recomD[c("TX","recom")]
colnames(recomD)[1]<-'actual'
recomD <- cbind(recomD,dana_outcome)
fit<-recom_plot(recomD,custom_theme)
#fit <-real_plot(recomD,custom_theme)


#-----------------------------------------
#            Loss E
#-----------------------------------------
recom <- refit_model_E$recommended.trts
actual <- df_disc$treatment
pid <- df_disc$patient_id
recomE <- cbind(pid,recom,actual)
recomE <- data.frame(recomE,disc_outcome)
fit<-recom_plot(recomE,custom_theme)
#fit<-real_plot(recomE,custom_theme)


test_cate_E <-refit_external(refit_model_E,df_dana)
recomE <-test_cate_E[[2]]
recomE <- recomE[c("TX","recom")]
colnames(recomE)[1]<-'actual'
recomE <- cbind(recomE,dana_outcome)
fit<-recom_plot(recomE,custom_theme)
#fit <-real_plot(recomE,custom_theme)



#-----------------------------------------
#            Composite
#-----------------------------------------
#x <- data.frame(0.10000000, 0.40000000, 0.10000000, 0.1129036, 0.48948990)
#x <- t(x)
#colnames(x) <- "w"
composite <- att_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_disc,'disc',x)
recomF <-composite[[2]]
drops <- c('y','wts')
recomF <- recomF[ , !(names(recomF) %in% drops)]
colnames(recomF)[3] <- "actual"
recomF <- data.frame(recomF,disc_outcome)
#debug(recom_plot)
fit<-recom_plot(recomF,custom_theme)
#fit<-real_plot(recomF,custom_theme)


composite <- att_scheme(refit_model_A,refit_model_B,refit_model_C,refit_model_D,refit_model_E,df_dana,'test',x)
recomF <-composite[[2]]
drops <- c('y','wts')
recomF <- recomF[ , !(names(recomF) %in% drops)]
colnames(recomF)[3] <- "actual"
recomF <- data.frame(recomF,dana_outcome)
fit<-recom_plot(recomF,custom_theme)
#fit <-real_plot(recomF,custom_theme)





