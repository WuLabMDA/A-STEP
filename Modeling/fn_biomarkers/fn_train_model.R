
#------------------------------------------
#        with feature selection
#-----------------------------------------
cross_validation <- function(df_disc,matched_ids,x.varnames,style,numFold,loss_type)
{
  train_loss <- data.frame(matrix(ncol = 2, nrow = numFold))
  colnames(train_loss) <- c('CATE0', 'CATE1')
  train_sample <- data.frame(matrix(ncol = 2, nrow = numFold))
  colnames(train_sample) <- c('CATE0', 'CATE1')
  
  valid_loss <- data.frame(matrix(ncol = 2, nrow = numFold))
  colnames(valid_loss) <- c('CATE0', 'CATE1')
  valid_sample <- data.frame(matrix(ncol = 2, nrow = numFold))
  colnames(valid_sample) <- c('CATE0', 'CATE1')
  
  trained_models <- list("m1","m2","m3","m4","m5")
  recursive_select <- list("m1","m2","m3","m4","m5")
  
  train_idx <- list("m1","m2","m3","m4","m5")
  valid_idx <- list("m1","m2","m3","m4","m5")
  
  
  
  folds <- cut(seq(1,nrow(matched_ids)),breaks=numFold,labels=FALSE)
  
  #---K Fold CV-----
  for (i in 1:numFold)
    
    {
     print(sprintf("Fold : %d",i))
     cate0 <- 0
     cate1 <- 0
     testIndexes <-which(folds==i,arr.ind=TRUE) 
     valid_ids <- matched_ids[testIndexes, ]
     train_ids <- matched_ids[-testIndexes, ]
     
     #---train set----
     control_id <- train_ids['Control_PID']
     treated_id <- train_ids['Treated_PID']
     control_set <- filter(df_disc,patient_id %in% control_id$Control_PID)
     treated_set <- filter(df_disc,patient_id %in% treated_id$Treated_PID)
     train_set <- rbind(control_set,treated_set)
     
     #---valid set----
     control_id <- valid_ids['Control_PID']
     treated_id <- valid_ids['Treated_PID']
     control_set <- filter(df_disc,patient_id %in% control_id$Control_PID)
     treated_set <- filter(df_disc,patient_id %in% treated_id$Treated_PID)
     valid_set <- rbind(control_set,treated_set)
     
     y <- train_set$prog_3_mo 
     trt <- train_set$treatment 
     features <- train_set[, x.varnames]
     rand.select <-sample(x.varnames,1)
     new.varnames <- c(rand.select)
     
     
     for (k in 1:length(features))
     {
       if(x.varnames[k] != rand.select)
       {
         sel_features <-data.frame(features[,k])
         colnames(sel_features) <- x.varnames[k]
         interest_features <-cbind(features[new.varnames],sel_features)
         x <- model.matrix(~ -1 + ., data = interest_features)
         
         #----fitting------
         set.seed(123)
         trained_model <-fit.subgroup(
         x=x,
         y= y,
         trt=trt,
         propensity.func = NULL,
         loss = loss_type,
         method = "weighting",
         match.id = NULL,
         augment.func = NULL,
         fit.custom.loss = NULL,
         cutpoint = 0,
         larger.outcome.better = FALSE,
         reference.trt = NULL,
         retcall = TRUE,
         verbose = FALSE)
       
       #---calculate CATE on trained model---
       train_cate <- adjusted_cate(trained_model,"train")
       valid <-eval_model(trained_model,valid_set)
       temp <- calc_effects_table(valid)
       valid_cate <- adjusted_cate(temp,"valid")
       if (train_cate[[1]] < cate0 & train_cate[[2]] < cate1)
        {
         new.varnames <-c(new.varnames,x.varnames[k])
          cate0 <-valid_cate[[1]]
          cate1 <-valid_cate[[2]]
        }
      } 
     }
     
     
     refit_cate <- refit_features(train_set,valid_set,new.varnames,style,loss_type,x.varnames)
     train_cate <-refit_cate[[1]]
     train_loss[i,1] <- round(train_cate[[1]],digits=2)
     train_loss[i,2] <- round(train_cate[[2]],digits=2)
     train_sample[i,1] <- train_cate[[3]]
     train_sample[i,2] <- train_cate[[4]]
     
     valid_cate <-refit_cate[[2]]
     valid_loss[i,1] <- round(valid_cate[[1]],digits=2)
     valid_loss[i,2] <- round(valid_cate[[2]],digits=2)
     valid_sample[i,1] <- valid_cate[[3]]
     valid_sample[i,2] <- valid_cate[[4]]
     
     trained_models[[i]]<-refit_cate[[3]]
     train_idx[[i]] <- train_ids
     valid_idx[[i]] <- valid_ids
     recursive_select[[i]]<-new.varnames

     
  }
  return(list(train_loss,valid_loss,trained_models,train_idx,valid_idx,recursive_select,train_sample,valid_sample))
}

#------------------------------------------
#        Refit K-fold model
#-----------------------------------------
refit_features <-function(train_set,valid_set,new.varnames,style,loss_type,x.varnames)
  {
  y <- train_set$prog_3_mo 
  trt <- train_set$treatment 
  features_refit <- train_set[, new.varnames]
  features_refit <- data.frame(features_refit)
  x <- model.matrix(~ -1 + ., data = features_refit)
  set.seed(123)
  refit_model <-fit.subgroup(
    x=x,
    y= y,
    trt=trt,
    propensity.func = NULL,
    loss = loss_type,
    method = "weighting",
    match.id = NULL,
    augment.func = NULL,
    fit.custom.loss = NULL,
    cutpoint = 0,
    larger.outcome.better = FALSE,
    reference.trt = NULL,
    retcall = TRUE,
    verbose = FALSE
  )
  
  
  if (is_empty(valid_set))
  {
    train_cate <- adjusted_cate(refit_model,"train")
    return (list(train_cate,refit_model))
  }
  else
  {
    train_cate <- adjusted_cate(refit_model,"train")
    valid <-eval_model(refit_model,valid_set)
    temp <- calc_effects_table(valid)
    valid_cate <- adjusted_cate(temp,"valid")
    return (list(train_cate,valid_cate,refit_model))
  }

  
  
 }




#------------------------------------------
#         adjusted CATE
#-----------------------------------------

adjusted_catexx <- function(model,set)
{
  if (set == 'train')
  {
    E00 = model$subgroup.trt.effects$avg.outcomes[1]
    if (E00 == 'NaN'){E00 <- 0.001}
    E01 = model$subgroup.trt.effects$avg.outcomes[2]
    if (E01 == 'NaN'){E01 <- 0.001}
    
    E10 = model$subgroup.trt.effects$avg.outcomes[3]
    if (E10 == 'NaN'){E10 <- 0.001}
    E11 = model$subgroup.trt.effects$avg.outcomes[4]
    if (E11 == 'NaN'){E11 <- 0.001}
    
    conditional_effect_G0 <- (E00 - E01) 
    conditional_effect_G1 <- (E11 - E10) 
    
    sum_g0 = model$subgroup.trt.effects$sample.sizes[1] + model$subgroup.trt.effects$sample.sizes[2]
    sum_g1 = model$subgroup.trt.effects$sample.sizes[3] + model$subgroup.trt.effects$sample.sizes[4]
   
    
    n00 <- model$subgroup.trt.effects$sample.sizes[1]/sum_g0
    n01 <- model$subgroup.trt.effects$sample.sizes[2]/sum_g0
    
    n10 <- model$subgroup.trt.effects$sample.sizes[3]/sum_g1
    n11 <- model$subgroup.trt.effects$sample.sizes[4]/sum_g1
    
    conditional_effect_G0 <- conditional_effect_G0/(n00/n01)
    if (conditional_effect_G0 == 'NaN'){conditional_effect_G0 <- 0.001}
    conditional_effect_G1 <- conditional_effect_G1/(n11/n10)
    if (conditional_effect_G1 == 'NaN'){conditional_effect_G1 <- 0.001}
    
    
    conditional_effect_G0 <- conditional_effect_G0*100
    conditional_effect_G1 <- conditional_effect_G1*100
    
    return(list(conditional_effect_G0,conditional_effect_G1))
  }
  else
  {
    E00 = model$avg.outcomes[1]
    if (E00 == 'NaN'){E00 <- 0.001}
    E01 = model$avg.outcomes[2]
    if (E01 == 'NaN'){E01 <- 0.001}
    
    E10 = model$avg.outcomes[3]
    if (E10 == 'NaN'){E10 <- 0.001}
    E11 = model$avg.outcomes[4]
    if (E11 == 'NaN'){E11 <- 0.001}
    
    conditional_effect_G0 <- (E00 - E01) 
    conditional_effect_G1 <- (E11 - E10) 
    
    sum_g0 = model$sample.sizes[1] + model$sample.sizes[2]
    sum_g1 = model$sample.sizes[3] + model$sample.sizes[4]
    
    n00 <- model$sample.sizes[1]/sum_g0
    n01 <- model$sample.sizes[2]/sum_g0
    
    n10 <- model$sample.sizes[3]/sum_g1
    n11 <- model$sample.sizes[4]/sum_g1
    
    conditional_effect_G0 <- conditional_effect_G0/(n00/n01)
    if (conditional_effect_G0 == 'NaN'){conditional_effect_G0 <- 0.001}
    conditional_effect_G1 <- conditional_effect_G1/(n11/n10)
    if (conditional_effect_G1 == 'NaN'){conditional_effect_G1 <- 0.001}
    
    conditional_effect_G0 <- conditional_effect_G0*100
    conditional_effect_G1 <- conditional_effect_G1*100
    
    return(list(conditional_effect_G0,conditional_effect_G1))
    
  }
  
}


#------------------------------------------
#         Original CATE
#-----------------------------------------

adjusted_cate <- function(model,set)
{
  if (set == 'train')
  {
    E00 = model$subgroup.trt.effects$avg.outcomes[1] 
    if (E00 == 'NaN'){E00 <- 0.001}
    E01 = model$subgroup.trt.effects$avg.outcomes[2] 
    if (E01 == 'NaN'){E01 <- 0.001}
    
    E10 = model$subgroup.trt.effects$avg.outcomes[3] 
    if (E10 == 'NaN'){E10 <- 0.001}
    E11 = model$subgroup.trt.effects$avg.outcomes[4] 
    if (E11 == 'NaN'){E11 <- 0.001}
    
    conditional_effect_G0 <- (E00 - E01) 
    conditional_effect_G1 <- (E11 - E10) 
    
    sum_g0 = model$subgroup.trt.effects$sample.sizes[1] + model$subgroup.trt.effects$sample.sizes[2]
    sum_g1 = model$subgroup.trt.effects$sample.sizes[3] + model$subgroup.trt.effects$sample.sizes[4]
    
    
    conditional_effect_G0 <- conditional_effect_G0*100
    conditional_effect_G1 <- conditional_effect_G1*100
    
    return(list(conditional_effect_G0,conditional_effect_G1,sum_g0,sum_g1))
  }
  else
  {
    E00 = model$avg.outcomes[1] 
    if (E00 == 'NaN'){E00 <- 0.001}
    E01 = model$avg.outcomes[2] 
    if (E01 == 'NaN'){E01 <- 0.001}
    
    E10 = model$avg.outcomes[3] 
    if (E10 == 'NaN'){E10 <- 0.001}
    E11 = model$avg.outcomes[4] 
    if (E11 == 'NaN'){E11 <- 0.001}
    
    conditional_effect_G0 <- (E00 - E01) 
    conditional_effect_G1 <- (E11 - E10) 
    
    sum_g0 = model$sample.sizes[1] + model$sample.sizes[2]
    sum_g1 = model$sample.sizes[3] + model$sample.sizes[4]
    
    conditional_effect_G0 <- conditional_effect_G0*100
    conditional_effect_G1 <- conditional_effect_G1*100
    
    return(list(conditional_effect_G0,conditional_effect_G1,sum_g0,sum_g1))
    
  }
  
}


#------------------------------------------
#           Print subgroup effects
#-----------------------------------------
print_effects_table <- function(model)
{
  
  digits = max(getOption('digits')-3, 3)
  Cf <- matrix(paste0(round(model$avg.outcomes, digits),
                      " (n = ", model$sample.sizes, ")"), ncol = ncol(model$avg.outcomes))
  
  
  dimnames(Cf) <- dimnames(model$avg.outcomes)
  
  print('=============Average Outcomes===========')
  #cat("Average Outcomes:\n")
  print(Cf, quote = FALSE, right = TRUE, na.print = "NA")
  cat("\n")
  print('=======Treatment effects conditional on subgroups=========')
  
  E00 = model$avg.outcomes[1]
  E01 = model$avg.outcomes[2]
  
  E10 = model$avg.outcomes[3]
  E11 = model$avg.outcomes[4]
  
  conditional_effect_G0 <- (E00 - E01) 
  conditional_effect_G1 <- (E11 - E10) 
  
  sum_g0 = model$sample.sizes[1] + model$sample.sizes[2]
  sum_g1 = model$sample.sizes[3] + model$sample.sizes[4]
  
  conditional_effect_G0 <- conditional_effect_G0/(model$sample.sizes[1]/model$sample.sizes[2])
  conditional_effect_G1 <- conditional_effect_G1/(model$sample.sizes[4]/model$sample.sizes[3])
  
  conditional_effect_G0 <- conditional_effect_G0*100
  conditional_effect_G1 <- conditional_effect_G1*100
  
  
  Cf2 <- paste0(round(conditional_effect_G0, digits)," (n = ", sum_g0, ")", "      ",   
                round(conditional_effect_G1, digits)," (n = ", sum_g1, ")")
  
  Cf2 <- data.frame(Cf2)
  names(Cf2) <- paste0("Recommended_0","       ",",Recommended_1")
  print(Cf2, quote = FALSE, right = TRUE, na.print = "NA")
  
}




#------------------------------------------
#           Create weights (original package)
#-----------------------------------------

create.weights.binary.trt <- function(pi.x, trt, method)
{
  
  if (is.factor(trt))
  {
    # drop any unused levels of trt
    trt         <- droplevels(trt)
    unique.trts <- levels(trt)
    n.trts      <- length(unique.trts)
  } else
  {
    unique.trts <- sort(unique(trt))
    n.trts      <- length(unique.trts)
  }
  
  if (n.trts != 2) stop("two trtment levels only for binary trt weighting function")
  
  if (method == "weighting")
  {
    wts     <- 1 / (pi.x * (trt == unique.trts[2L]) + (1 - pi.x) * (trt == unique.trts[1L]))
  } else
  {   # A-learning method
    wts     <- rep(1, length(pi.x))
  }
  wts
}





#------------------------------------------
#           composite 4x4 effects
#-----------------------------------------
calc_effects_table <- function(df)
{
  n.trts <- length(unique(df$TX))
  recommended.trt <- df$recom
  unique.trts <- sort(unique(df$TX))
  wts <- df$wts
  trt <-df$T
  y <- df$y
    
  idx.list <- rep(list(vector(mode = "list", length = n.trts)), n.trts)
  
  for (t.recom in 1:n.trts)
  {
    for (t.receiv in 1:n.trts)
    {
      idx.list[[t.recom]][[t.receiv]] <- (recommended.trt == unique.trts[t.recom]) &
        (trt == unique.trts[t.receiv])
    }
  }
  
  res.mat <- matrix(0, ncol = n.trts, nrow = n.trts)
  colnames(res.mat) <- paste("Recommended", unique.trts)
  rownames(res.mat) <- paste("Received", unique.trts)
  
  sample.size.mat <- res.mat
  subgroup.effects <- numeric(n.trts)
  
  for (t.recom in 1:n.trts)
  {
    for (t.receiv in 1:n.trts)
    {
      idx.cur <- idx.list[[t.recom]][[t.receiv]]
      res.mat[t.receiv, t.recom] <- weighted.mean(y[idx.cur], w = wts[idx.cur])
      sample.size.mat[t.receiv, t.recom] <- sum(idx.cur)
      
      if (t.recom == t.receiv)
      {
        #idx.disagree <- Reduce("|", idx.list[[t.recom]][-t.receiv])
        idx.disagree <- (recommended.trt == unique.trts[t.recom]) &
          (trt != unique.trts[t.recom])
        subgroup.effects[t.recom] <- res.mat[t.receiv, t.recom] - weighted.mean(y[idx.disagree], w = wts[idx.disagree])
      }
      
    }
  }
  
  idx.agree <- recommended.trt == trt
  overall.subgroup.effect <- weighted.mean(y[idx.agree], w = wts[idx.agree]) - weighted.mean(y[!idx.agree], w = wts[!idx.agree])
  
  fitted.model <- list(avg.outcomes=res.mat,sample.sizes=sample.size.mat,overall.subgroup.effect=overall.subgroup.effect)
  return(fitted.model)
  
  
  
}


