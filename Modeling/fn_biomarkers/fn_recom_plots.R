#------------------------------------------
#        recommended by arms
#-----------------------------------------
recom_plot <- function(df,custom_theme)
{
  #--------------------------------
  # Recommended ICI-mono
  #--------------------------------
  follow_T0 <- df[df$recom==0 & df$actual==0,]
  follow_T0$Group <-'Follow'
  anti_T0 <- df[df$recom==0 & df$actual==1,]
  anti_T0$Group <-'Anti'
  df_surv <- rbind(follow_T0,anti_T0)
  rownames(df_surv) <-NULL
  fit <- survfit(Surv(PFS,PFS_events)~ Group, data = df_surv)
  med1<-surv_median(fit)
  
  print(ggsurvplot(fit, data = df_surv,title = "Recom ICI-mono",ggtheme=custom_theme(),
             conf.int = FALSE,
             pval = TRUE,
             fun = "pct",
             risk.table = TRUE,
             risk.table.fontsize =5,
             size = 1,
             xlab = "Time (months)",
             ylab = "PFS (%)",
             #xlim = c(0, 5),
             linetype = "strata",
             palette = c("#6666FF",
                         "#FF3366"),
             risk.table.col = "strata",
             #legend = "bottom",
             legend.title = "Risk",
             legend.labs = c("Anti",
                             "Follow")))
  p1<-summary(coxph(Surv(df_surv$PFS, df_surv$PFS_events)~relevel(as.factor(df_surv$Group),ref='Anti')))
  #pi1 <-summary(coxph(Surv(df_surv$PFS, df_surv$PFS_events)~df_surv$Group*df_surv$actual))
  #pi1 <-pi1$coefficients[15]
  
  
  #--------------------------------
  # Recommended ICI-chemo
  #--------------------------------
  follow_T1 <- df[df$recom==1 & df$actual==1,]
  follow_T1$Group <-'Follow'
  anti_T1 <- df[df$recom==1 & df$actual==0,]
  anti_T1$Group <-'Anti'
  df_surv <- rbind(follow_T1,anti_T1)
  rownames(df_surv) <-NULL
  fit <- survfit(Surv(PFS,PFS_events)~ Group, data = df_surv)
  med2<-surv_median(fit)
  
  print(ggsurvplot(fit, data = df_surv,title = "Recom ICI-chemo",ggtheme=custom_theme(),
             conf.int = FALSE,
             pval = TRUE,
             fun = "pct",
             risk.table = TRUE,
             risk.table.fontsize =5,
             size = 1,
             xlab = "Time (months)",
             ylab = "PFS (%)",
             #xlim = c(0, 5),
             linetype = "strata",
             palette = c("#FF3366",
                         "#6666FF"),
             risk.table.col = "strata",
             #legend = "bottom",
             legend.title = "Risk",
             legend.labs = c("Anti",
                             "Follow")))
  p2<-summary(coxph(Surv(df_surv$PFS, df_surv$PFS_events)~relevel(as.factor(df_surv$Group),ref='Anti')))
  #pi2 <-summary(coxph(Surv(df_surv$PFS, df_surv$PFS_events)~df_surv$Group*df_surv$actual))
  #pi2 <-pi2$coefficients[15]
  
  
  
  return(list(p1,p2,med1,med2))
}



#------------------------------------------
#        Real TX by arms
#-----------------------------------------
real_plot <- function(df,custom_theme)
{
  #--------------------------------
  # Received ICI-mono
  #--------------------------------
  follow_T0 <- df[df$actual==0 & df$recom==0,]
  follow_T0$Group <-'Follow'
  anti_T0 <- df[df$actual==0 & df$recom==1,]
  anti_T0$Group <-'Anti'
  df_surv <- rbind(follow_T0,anti_T0)
  rownames(df_surv) <-NULL
  fit <- survfit(Surv(PFS,PFS_events)~ Group, data = df_surv)
  med1<-surv_median(fit)
  print(ggsurvplot(fit, data = df_surv,title = "Received ICI-mono",ggtheme=custom_theme(),
                   conf.int = FALSE,
                   pval = TRUE,
                   fun = "pct",
                   risk.table = TRUE,
                   risk.table.fontsize =5,
                   size = 1,
                   xlab = "Time (months)",
                   ylab = "PFS (%)",
                   #xlim = c(0, 5),
                   linetype = "strata",
                   palette = c("Blue",
                               "Red3"),
                   risk.table.col = "strata",
                   #legend = "bottom",
                   legend.title = "Risk",
                   legend.labs = c("ICI-Chemo",
                                   "ICI-Mono")))
  p1<-summary(coxph(Surv(df_surv$PFS, df_surv$PFS_events)~relevel(as.factor(df_surv$Group),ref='Anti')))

  
  
  #--------------------------------
  # Received ICI-chemo
  #--------------------------------
  follow_T1 <- df[df$actual==1 & df$recom==1,]
  follow_T1$Group <-'Follow'
  anti_T1 <- df[df$actual==1 & df$recom==0,]
  anti_T1$Group <-'Anti'
  df_surv <- rbind(follow_T1,anti_T1)
  rownames(df_surv) <-NULL
  fit <- survfit(Surv(PFS,PFS_events)~ Group, data = df_surv)
  med2<-surv_median(fit)
  
  print(ggsurvplot(fit, data = df_surv,title = "Received ICI-chemo",ggtheme=custom_theme(),
                   conf.int = FALSE,
                   pval = TRUE,
                   fun = "pct",
                   risk.table = TRUE,
                   risk.table.fontsize =5,
                   size = 1,
                   xlab = "Time (months)",
                   ylab = "PFS (%)",
                   #xlim = c(0, 5),
                   linetype = "strata",
                   palette = c("Red3",
                               "Blue"),
                   risk.table.col = "strata",
                   #legend = "bottom",
                   legend.title = "Risk",
                   legend.labs = c("ICI-Mono",
                                   "ICI-Chemo")))
  p2<-summary(coxph(Surv(df_surv$PFS, df_surv$PFS_events)~relevel(as.factor(df_surv$Group),ref='Anti')))

  
  
  return(list(p1,p2,med1,med2))
}

#------------------------------------------
#        P-interaction
#-----------------------------------------
p_interaction <- function(df)
{
  p <-summary(coxph(Surv(df$PFS, df$PFS_events)~df$recom*df$actual))
  return(p$coefficients[15])
    
}
  