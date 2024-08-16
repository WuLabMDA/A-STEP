rm(list=ls(all=TRUE))
setwd("C:\\Users\\MBSaad\\Desktop\\Projects\\QASEM_NC\\Version_3(Repeat_CV)\\Modeling\\CV_lock_features")
load('trained_models.RData')

#---Data reading-----
#df_mda <- read.csv("Matched_MDA.csv")
#df_mayo <- read.csv("Matched_MAYO.csv")
#drops <- c('Liver.met','Brain.met','Met.status')
#df_mda <- df_mda[ , !(names(df_mda) %in% drops)]
#df_mayo <- df_mayo[ , !(names(df_mayo) %in% drops)]

#df_natgen <-read.csv("Matched_NATGEN.csv")
#drops <- c('OS','OS_events','PFS','PFS_events')
#df_natgen <- df_natgen[ , !(names(df_natgen) %in% drops)]
#df_disc <-rbind(df_mda,df_mayo,df_natgen)


model <- refit_cate_C[[2]]
b.scores <- model$benefit.scores
subset <-model$var.names
features <- df_disc[subset]
shap_data <- cbind(features,b.scores)

write.csv(shap_data,'SHAP_C.csv')