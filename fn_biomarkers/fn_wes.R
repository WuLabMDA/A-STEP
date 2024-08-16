WeightedEnsemble<-function(model_A,model_B,model_C,model_D,model_E,df_data,Method=optimizer,cohort,test_data=NULL, forecast=NULL)
{ 
  Optimization<-Method
  cutpoint <-0
  input_data <- data.frame(model_A$benefit.scores,model_B$benefit.scores,model_C$benefit.scores,
                   model_D$benefit.scores,model_E$benefit.scores)
  numVar <- ncol(input_data)
  rangeVar <- matrix(c(0,1), nrow=2)

  
  if (cohort == 'disc')
  {
    sphere <- function(x)
    {
      df <- data.frame(input_data[1]*x[1],input_data[2]*x[2],input_data[3]*x[3],
                       input_data[4]*x[4],input_data[5]*x[5])
      
      
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
      #debug(adjusted_cate)
      data_cate <- adjusted_cate(temp,"refit")
      err0 <- data_cate[[1]]
      err1 <- data_cate[[2]]
      #error_fn <-err1
      #error_fn <-sum(err0,err1)
      error_fn <- (0.1*err0) + (0.9*err1)
      return(error_fn)
    }
  }
  
 if (Optimization=="PSO")
   {
  #debug(PSO)
  opt <- PSO(sphere, optimType="MIN", numVar, numPopulation=30,
             maxIter=1000,rangeVar)
  optimum_value <- sphere(opt)  
    }
  
  else if (Optimization=="GWO")
    {
    #debug(GWO)
    opt <- GWO(sphere, optimType="MIN", numVar, numPopulation=30,
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
    }
  
  else if (Optimization=="GA")
  {
    opt <- GA(sphere, optimType="MIN", numVar, numPopulation=30,
               maxIter=1000,rangeVar)
    optimum_value <- sphere(opt)
  }
  
    else {
      message("Please select a valid optimization technique")
    
       }


  return(list(Weights=opt,optimum_value))
}

