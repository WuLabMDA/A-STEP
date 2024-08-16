#------------------------------------------
#        External test
#-----------------------------------------
refit_external <- function(trained_model,test_set)
{
  model <- trained_model
  test <-eval_model(model,test_set)
  temp <- calc_effects_table(test)
  temp_cate <- adjusted_cate(temp,"test")
  return(list(temp_cate,test))
}


#------------------------------------------
#           Evaluate on external data
#-----------------------------------------
eval_model <- function(model,df)
{

  sel_varnames <-model$var.names
  if (length(sel_varnames)==1)
  {
    sel_varnames <-sample(x.varnames,1)
  }
  
  y <- df$prog_3_mo
  TX <- df$treatment 
  sel_features <- df[, sel_varnames]
  sel_features <- data.frame(sel_features)
  x <- model.matrix(~ -1 + ., data = sel_features)
    
  bene_score <- predict(model, newx = x, type = "benefit.score")
  recom <- predict(model, newx = x, type = "trt.group")
  pi.x <-rep(c(0.5),times=nrow(df))
  wts <-rep(c(2),times=nrow(df))
  delta <- wts*bene_score
    
  results <- data.frame(TX,recom,y,delta,bene_score,pi.x,wts)
  return(results)
  

}

#------------------------------------------
#           Avegraging CATE
#-----------------------------------------
avg_cate <- function(train_cate,valid_cate,test_cate,name)
{ 
  train_cate0 <- round(train_cate[[1]],2)
  train_cate1 <- round(train_cate[[2]],2)
  test_cate0 <- round(test_cate[[1]],2)
  test_cate1 <- round(test_cate[[2]],2)
  
  train_sample0 <- train_cate[[3]]
  train_sample1 <- train_cate[[4]]
  test_sample0 <- test_cate[[3]]
  test_sample1 <- test_cate[[4]]
  
  if (is_empty(valid_cate))
  {
    temp <- data.frame(set = c('Disc','DFCI'),
                       CATE_T0 = c(paste0(train_cate0," ", "(n=",train_sample0,")"),paste0(test_cate0," ", "(n=",test_sample0,")")),
                       CATE_T1 = c(paste0(train_cate1," ", "(n=",train_sample1,")"),paste0(test_cate1," ", "(n=",test_sample1,")")))
    colnames(temp)[1] <-name
    return(temp)
  }
  else
  {
    valid_cate0 <- mean(valid_cate[[1]])
    valid_cate1 <- mean(valid_cate[[2]])
    
    temp <- data.frame(set = c('Train','valid','DFCI'),
                       CATE_T0 = c(train_cate0,valid_cate0,test_cate0),
                       CATE_T1 = c(train_cate1,valid_cate1,test_cate1))
    colnames(temp)[1] <-name
    return(temp)
  }
  
  
}


#-----------------------------------------
#           Averaging Kfold results
#-----------------------------------------
avg_iter_results <- function(train_cate,valid_cate,train_sample,valid_sample)
{
  train_cate0 <- round(mean(train_cate[[1]]),1)
  train_cate1 <- round(mean(train_cate[[2]]),1)
  valid_cate0 <- round(mean(valid_cate[[1]]),1)
  valid_cate1 <- round(mean(valid_cate[[2]]),1)
  
  train_sample0 <- round(mean(train_sample[[1]],0))
  train_sample1 <- round(mean(train_sample[[2]],0))
  valid_sample0 <- round(mean(valid_sample[[1]],0))
  valid_sample1 <- round(mean(valid_sample[[2]],0))
  
  cf1 <- paste0(train_cate0," (n = ", train_sample0, ")", "      ",train_cate1," (n = ", train_sample1, ")")
  cf2 <- paste0(valid_cate0," (n = ", valid_sample0, ")", "      ",valid_cate1," (n = ", valid_sample1, ")")
  
  cf1 <- data.frame(cf1)
  names(cf1) <- paste0("Recommended_0","       ",",Recommended_1")
  
  cf2 <- data.frame(cf2)
  names(cf2) <- paste0("Recommended_0","       ",",Recommended_1")
  
  cf <- rbind(cf1,cf2)
  rownames(cf) <-c('Train','valid')

  print(cf, quote = FALSE, right = TRUE)
  
}