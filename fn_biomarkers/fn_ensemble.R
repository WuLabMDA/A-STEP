#-----------------------------------------
#            Composite : Voting
#-----------------------------------------
voting_scheme <- function(model_A,model_B,model_C,model_D,model_E,df_data,cohort,T)
{
  if (cohort == 'disc')
    {
      df <- data.frame(model_A$recommended.trts,model_B$recommended.trts,model_C$recommended.trts,model_D$recommended.trts,model_E$recommended.trts)
      df$recom <- rowSums(df)
      df$recom <- df$recom>T
      df$recom <- as.integer(as.logical(df$recom))
      df$pid <- df_data$patient_id
      df <- df[c("pid","recom")]
      
      df$TX <- model_B$trt.received
      df$y <- model_B$y
      df$wts <- model_B$model$call$weights
      temp <-calc_effects_table(df)
      data_cate <- adjusted_cate(temp,"refit")

    }
  else
    {
      valid_A <-eval_model(model_A,df_data)
      valid_B <-eval_model(model_B,df_data)
      valid_C <-eval_model(model_C,df_data)
      valid_D <-eval_model(model_D,df_data)
      valid_E <-eval_model(model_E,df_data)
      
      df <- data.frame(valid_A$recom,valid_B$recom,valid_C$recom,valid_D$recom,valid_E$recom)
      df$recom <- rowSums(df)
      df$recom <- df$recom>T
      df$recom <- as.integer(as.logical(df$recom))
      df$pid <- df_data$patient_id
      df <- df[c("pid","recom")]
      
      df$TX <- df_data$treatment 
      df$y <- df_data$prog_3_mo
      df$wts <-rep(c(2),times=nrow(df_data))
      temp <-calc_effects_table(df)
      data_cate <- adjusted_cate(temp,"refit")
    }
  
  return (data_cate)
}

#-----------------------------------------
#            Composite : AVERAGING
#-----------------------------------------
avg_scheme <- function(model_A,model_B,model_C,model_D,model_E,df_data,cohort)
{ cutpoint <-0
  if (cohort == 'disc')
  { 
    df <- data.frame(model_A$benefit.scores,model_B$benefit.scores,model_C$benefit.scores,model_D$benefit.scores,model_E$benefit.scores)
    df$recom <- rowMeans(df)
    df$recom <- df$recom<cutpoint
    df$recom <- as.integer(as.logical(df$recom))
    df$pid <- df_data$patient_id
    df <- df[c("pid","recom")]
    
    df$TX <- model_B$trt.received
    df$y <- model_B$y
    df$wts <- model_B$model$call$weights
    temp <-calc_effects_table(df)
    data_cate <- adjusted_cate(temp,"refit")
    
  }
  else
  {
    valid_A <-eval_model(model_A,df_data)
    valid_B <-eval_model(model_B,df_data)
    valid_C <-eval_model(model_C,df_data)
    valid_D <-eval_model(model_D,df_data)
    valid_E <-eval_model(model_E,df_data)
    
    df <- data.frame(valid_A$bene_score,valid_B$bene_score,valid_C$bene_score,valid_D$bene_score,valid_E$bene_score)
    df$recom <- rowMeans(df)
    df$recom <- df$recom<cutpoint
    df$recom <- as.integer(as.logical(df$recom))
    df$pid <- df_data$patient_id
    df <- df[c("pid","recom")]
    
    df$TX <- df_data$treatment 
    df$y <- df_data$prog_3_mo
    df$wts <-rep(c(2),times=nrow(df_data))
    temp <-calc_effects_table(df)
    data_cate <- adjusted_cate(temp,"refit")
  }
  
  return (data_cate)
}


#-----------------------------------------
#            Composite : ATTENTION
#-----------------------------------------
att_scheme <- function(model_A,model_B,model_C,model_D,model_E,df_data,cohort,x)
{ cutpoint <-0
  if (cohort == 'disc')
  {
    df <- data.frame(model_A$benefit.scores*x[1],model_B$benefit.scores*x[2],model_C$benefit.scores*x[3],
                     model_D$benefit.scores*x[4],model_E$benefit.scores*x[5])
    sum_wts <-sum(x)
    df$recom <- rowSums(df)
    df$recom <- df$recom/sum_wts
    df$recom <- df$recom<cutpoint
    df$recom <- as.integer(as.logical(df$recom))
    df$pid <- df_data$patient_id
    df <- df[c("pid","recom")]
    
    df$TX <- model_B$trt.received
    df$y <- model_B$y
    df$wts <- model_B$model$call$weights

    temp <-calc_effects_table(df)
    data_cate <- adjusted_cate(temp,"refit")
  }

  else
  {
    valid_A <-eval_model(model_A,df_data)
    valid_B <-eval_model(model_B,df_data)
    valid_C <-eval_model(model_C,df_data)
    valid_D <-eval_model(model_D,df_data)
    valid_E <-eval_model(model_E,df_data)
  
    df <- data.frame(valid_A$bene_score*x[1],valid_B$bene_score*x[2],valid_C$bene_score*x[3],valid_D$bene_score*x[4],valid_E$bene_score*x[5])
    sum_wts <-sum(x)
    df$recom <- rowSums(df)
    df$recom <- df$recom/sum_wts
    df$recom <- df$recom<cutpoint
    df$recom <- as.integer(as.logical(df$recom))
    df$pid <- df_data$patient_id
    df <- df[c("pid","recom")]
  
    df$TX <- df_data$treatment 
    df$y <- df_data$prog_3_mo
    df$wts <-rep(c(2),times=nrow(df_data))
    temp <-calc_effects_table(df)
    data_cate <- adjusted_cate(temp,"refit")
  }

return (list(data_cate,df))
}

