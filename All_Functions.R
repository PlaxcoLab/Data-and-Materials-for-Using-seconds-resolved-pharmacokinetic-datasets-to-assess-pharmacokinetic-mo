#================================================================================================================================================================
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# PK_NLS_Estimation
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#================================================================================================================================================================

{
  
  # This function performs the NLS modeling for all the considered 
  # one-compartment models and the constant elimination rate 
  # two-compartment model. NOTE: the modified function "nlsLM" is used 
  # in place of base function "nls" as it has been more stable in 
  # practice. 
  #
  # Function inputs: input_data = the observed data containing at least a 
  #                               time column (labeled "time"), a plasma 
  #                               concentration column (labeled "conc.C_P"), 
  #                               and a column containing the initial dose 
  #                               (labeled "D_0")
  #                  model = specification for which model to fit to the 
  #                          data (i.e. C_1CM, L_1CM, E_1CM, R_1CM, C_2CM)
  #                  optim_method = All time-varying models use optim in 
  #                                 the fitting process, the allowable 
  #                                 method inputs are "Nelder-Mead" and 
  #                                 "SANN"
  #                  max_iter = the number of maximum iterations to allow
  #                             for in the NLS fitting
  #                  converge_crit = the convergence criteria for the NLS 
  #                                  fitting
  # 
  # Function outputs: "NLS Model" = the NLS model fit
  
  PK_NLS_Estimation <- function(input_data,
                                model,
                                optim_method,
                                max_iter = 50,
                                converge_crit = 1e-05){
    
    if (model=='C_1CM'){ # Constant-rate One-compartment model
      
      # grab initial dose from input_data
      
      D_0 <- min(input_data$D_0) # nmol
      
      #===========================================
      # Find initial estimates with lm by taking
      # the logarithm of the C_1CM solution
      #===========================================
      
      # center independent variable
      input_data$centered_time <- input_data$time - mean(input_data$time)
      
      # if observed data is negative, shift data to be positive
      if (min(input_data$conc.C_P) <= 0){
        log.C_P.lm <- lm(log(conc.C_P + 0.0001 - min(conc.C_P)) ~ centered_time,
                         data = input_data)
      }else{
        log.C_P.lm <- lm(log(conc.C_P) ~ centered_time, data = input_data)
      }
      
      # save initial estimates
      V_P <- unname(D_0/exp(coef(log.C_P.lm)[1])) # mL
      k_E <- unname(       -coef(log.C_P.lm)[2] ) # 1 / hr
      
      #===========================================
      # Define function C_P,
      # the solution to the
      # Constant-rate One-compartment model
      #===========================================
      
      C_P <- function(t,V_P,k_E){
        D_0/V_P * exp(-k_E * t)
      }
      
      #===========================================
      # Apply NLS Model using `nls`
      #===========================================
      
      NLS_Cp_Model <- nlsLM(conc.C_P ~ C_P(time,V_P,k_E),
                            data = input_data,
                            start = list(V_P = V_P,
                                         k_E = k_E),
                            control = list(maxiter = max_iter,
                                           tol = converge_crit))
      
    }
    else if (model=='L_1CM'){ # Linear-rate One-compartment model
      
      # grab initial dose from input_data
      
      D_0 <- min(input_data$D_0) # nmol
      
      #===========================================
      # Find initial estimates with lm by taking
      # the logarithm of the L_1CM solution
      #===========================================
      
      # center independent variable
      input_data$centered_time <- input_data$time - mean(input_data$time)
      
      # if observed data is negative, shift data to be positive
      if (min(input_data$conc.C_P) <= 0){
        log.C_P.lm <- lm(log(conc.C_P + 0.0001 - min(conc.C_P)) ~ poly(centered_time,2),
                         data = input_data)
      }else{
        log.C_P.lm <- lm(log(conc.C_P) ~ poly(centered_time,2), data = input_data)
      }
      
      # save initial estimates
      V_P  <- unname( D_0/exp( coef(log.C_P.lm)[1] ) ) # mL
      k_E0 <- unname(         -coef(log.C_P.lm)[2]   ) # 1 / hr
      b_E  <- unname(       -2*coef(log.C_P.lm)[3]   ) # 1 / hr^2
      
      #===========================================
      # Define function C_P,
      # the solution to the
      # Linear-rate One-compartment model
      #===========================================
      
      C_P <- function(t,V_P,k_E0,b_E){
        D_0/V_P * exp( -( b_E/2*t + k_E0 ) * t )
      }
      
      #===========================================
      # The function `optim` is used to find
      # secondary estimates prior to the NLS
      # estimation
      # First define penalty function (here, SSE)
      # to use in `optim`
      #===========================================
      
      C_P.optim <- function(func_inputs){
        C_P(func_inputs[1],func_inputs[2],func_inputs[3],func_inputs[4])
      }
      
      SSE <- function(theta)
      {
        sum( (input_data$conc.C_P - apply(X = cbind(input_data[,1],
                                                    theta[1],
                                                    theta[2],
                                                    theta[3]),
                                          MARGIN = 1,
                                          FUN = C_P.optim)
        )^2 )
      }
      
      #===========================================
      # Apply `optim` and grab estimates
      #===========================================
      
      optim.C_P <- optim(par = c(V_P = V_P, k_E0 = k_E0, b_E  = b_E),
                         fn = SSE,
                         method = optim_method)
      
      V_P  <- unname(optim.C_P$par[1]) # mL
      k_E0 <- unname(optim.C_P$par[2]) # 1 / hr
      b_E  <- unname(optim.C_P$par[3]) # 1 / hr^2
      
      #===========================================
      # Apply NLS Model using `nlsLM`
      #===========================================
      
      NLS_Cp_Model <- nlsLM(conc.C_P ~ C_P(time,V_P,k_E0,b_E),
                            data = input_data,
                            start = list(V_P = V_P,
                                         k_E0 = k_E0,
                                         b_E  = b_E),
                            control = list(maxiter = max_iter,
                                           tol = converge_crit))
      
    }
    else if (model=='E_1CM'){ # Exponential-rate One-compartment model
      
      # grab initial dose from input_data
      
      D_0 <- min(input_data$D_0) # nmol
      
      #===========================================
      # Find initial estimates with lm by taking
      # the logarithm of the E_1CM solution
      #===========================================
      
      # center independent variable
      input_data$centered_time <- input_data$time - mean(input_data$time)
      
      # if observed data is negative, shift data to be positive
      if (min(input_data$conc.C_P) <= 0){
        log.C_P.lm <- lm(log(conc.C_P + 0.0001 - min(conc.C_P)) ~ poly(centered_time,2),
                         data = input_data)
      }else{
        log.C_P.lm <- lm(log(conc.C_P) ~ poly(centered_time,2), data = input_data)
      }
      
      # save initial estimates
      V_P  <- unname(D_0/exp(coef(log.C_P.lm)[1])) # mL
      k_E0 <- unname(        coef(log.C_P.lm)[2] ) # 1 / hr
      b_E  <- unname(        coef(log.C_P.lm)[3] ) # 1 / hr
      
      #===========================================
      # Define function C_P,
      # the solution to the
      # Exponential-rate One-compartment model
      #===========================================
      
      C_P <- function(t,V_P,k_E0,b_E){
        D_0/V_P * exp(  (k_E0/b_E) * (  exp(-b_E*t)  -  1  )  )
      }
      
      #===========================================
      # The function `optim` is used to find
      # secondary estimates prior to the NLS
      # estimation
      # First define penalty function (here, SSE)
      # to use in `optim`
      #===========================================
      
      C_P.optim <- function(func_inputs){
        C_P(func_inputs[1],func_inputs[2],func_inputs[3],func_inputs[4])
      }
      
      SSE <- function(theta)
      {
        sum( (input_data$conc.C_P - apply(X = cbind(input_data[,1],
                                                    theta[1],
                                                    theta[2],
                                                    theta[3]),
                                          MARGIN = 1,
                                          FUN = C_P.optim)
        )^2 )
      }
      
      #===========================================
      # Apply `optim` and grab estimates
      #===========================================
      
      optim.C_P <- optim(par = c(V_P = V_P, k_E0 = k_E0, b_E  = b_E),
                         fn = SSE,
                         method = optim_method)
      
      V_P  <- unname(optim.C_P$par[1]) # mL
      k_E0 <- unname(optim.C_P$par[2]) # 1 / hr
      b_E  <- unname(optim.C_P$par[3]) # 1 / hr
      
      #===========================================
      # Apply NLS Model using `nlsLM`
      #===========================================
      
      NLS_Cp_Model <- nlsLM(conc.C_P ~ C_P(time,V_P,k_E0,b_E),
                            data = input_data,
                            start = list(V_P = V_P,
                                         k_E0 = k_E0,
                                         b_E  = b_E),
                            control = list(maxiter = max_iter,
                                           tol = converge_crit))
      
    }
    else if (model=='R_1CM'){ # Reciprocal-rate One-compartment model
      
      # grab initial dose from input_data
      
      D_0 <- min(input_data$D_0) # nmol
      
      #===========================================
      # Find initial estimates with lm by taking
      # the logarithm of the C_1CM solution
      #===========================================
      
      # center independent variable
      input_data$centered_time <- input_data$time - mean(input_data$time)
      
      # if observed data is negative, shift data to be positive
      if (min(input_data$conc.C_P) <= 0){
        log.C_P.lm <- lm(log(conc.C_P + 0.0001 - min(conc.C_P)) ~ poly(centered_time,2),
                         data = input_data)
      }else{
        log.C_P.lm <- lm(log(conc.C_P) ~ poly(centered_time,2), data = input_data)
      }
      
      # save initial estimates
      V_P  <- unname(D_0/exp(coef(log.C_P.lm)[1])) # mL
      k_E0 <- unname(        coef(log.C_P.lm)[2] ) # 1 / hr
      b_E  <- unname(    abs(coef(log.C_P.lm)[3])) # 1 / hr
      
      #===========================================
      # Define function C_P,
      # the solution to the
      # Reciprocal-rate One-compartment model
      #===========================================
      
      C_P <- function(t,V_P,k_E0,b_E){
        D_0/V_P * exp(  -(k_E0/b_E) * log(  b_E*t  +  1  )  )
      }
      
      C_P.cobyla <- function(func_inputs){
        
        C_P(func_inputs[1],func_inputs[2],func_inputs[3],func_inputs[4])
        
      }
      
      constraint <- function(time, y, theta, ...) {
        con    <- numeric(2)
        con[1] <- theta[3] * min(time) + 1 + 0.0001 # the 0.0001 ensures we do reach 0
        con[2] <- theta[3] * max(time) + 1 + 0.0001 # the 0.0001 ensures we do reach 0
        
      }
      
      SSE <- function(time,y,theta)
      {
        sum((y-apply(cbind(time,theta[1],theta[2],theta[3]),1,C_P.cobyla))^2)
      }
      
      #===========================================
      # Apply `cobyla` and grab estimates
      #===========================================
      
      # messages are suppressed as this function always outputs a comment about the package
      constrained_est.C_P <- suppressMessages(
        cobyla(x0      = c(V_P, k_E0, b_E),
               fn      = SSE,
               hin     = constraint,
               time    = input_data$time,
               y       = input_data$conc.C_P,
               control = list(xtol_rel = converge_crit, maxeval = max_iter))
      )
      
      V_P  <- constrained_est.C_P$par[1] # mL
      k_E0 <- constrained_est.C_P$par[2] # 1 / hr
      b_E  <- constrained_est.C_P$par[3] # 1 / hr
      
      #===========================================
      # Apply NLS Model using `nlsLM`
      #===========================================
      
      NLS_Cp_Model <- nlsLM(conc.C_P ~ C_P(time,V_P,k_E0,b_E),
                            data = input_data,
                            start = list(V_P = V_P,
                                         k_E0 = k_E0,
                                         b_E  = b_E),
                            upper = c(Inf, Inf, Inf),
                            lower = c(-Inf, -Inf, -1/max(input_data$time)),
                            control = list(maxiter = max_iter,
                                           tol = converge_crit))
      
    }
    else if (model=='C_2CM'){ # Constant-rate Two-compartment model
      
      # grab initial dose from input_data
      
      D_0 <- min(input_data$D_0) # nmol
      
      #===========================================
      # Find initial estimates with the function
      # `PK_NLS_Estimation.Seq` by fitting the 
      # C_1CM model
      #===========================================
      
      C_1CM_fit <- PK_NLS_Estimation.Seq(input_data = input_data,
                                     model = 'C_1CM',
                                     max_iter = max_iter,
                                     converge_crit = converge_crit)
      
      
      Approx.CI.df <- Approx.CI(model = 'C_1CM',
                            D_0 = D_0.temp,
                            NLS_Cp_Model = C_1CM_fit)

      A <- D_0/Approx.CI.df$`Estimates and Approximate Confidence Bounds`$Estimates[1]
      a <- Approx.CI.df$`Estimates and Approximate Confidence Bounds`$Estimates[2]
      B <- A
      b <- a
      
      #===========================================
      # Define function C_P,
      # the solution to the
      # Constant-rate Two-compartment model
      #===========================================
      
      C_P <- function(t,A,a,B,b){
        (A*exp(-(t*a))+B*exp(-(t*b)))
      }
      
      #===========================================
      # The function `optim` is used to find
      # secondary estimates prior to the NLS
      # estimation
      # First define penalty function (here, SSE)
      # to use in `optim`
      #===========================================
      
      C_P.optim <- function(func_inputs){
        C_P(func_inputs[1],func_inputs[2],func_inputs[3],func_inputs[4],func_inputs[5])
      }
      
      SSE <- function(theta)
      {
        sum( (input_data$conc.C_P - apply(X = cbind(input_data[,1],
                                                    theta[1],
                                                    theta[2],
                                                    theta[3],
                                                    theta[4]),
                                          MARGIN = 1,
                                          FUN = C_P.optim)
        )^2 )
      }
      
      #===========================================
      # Apply `optim` and grab estimates
      #===========================================
      
      optim.C_P <- optim(par = c(A = A, a = a, B  = B, b = b),
                         fn = SSE,
                         method = optim_method)
      
      A <- unname(optim.C_P$par[1])
      a <- unname(optim.C_P$par[2])
      B <- unname(optim.C_P$par[3])
      b <- unname(optim.C_P$par[4])
      
      #===========================================
      # Apply NLS Model using `nlsLM`
      #===========================================
      
      NLS_Cp_Model <- nlsLM(conc.C_P ~ C_P(time,A,a,B,b),
                            data = input_data,
                            start = list(A = A,
                                         a = a,
                                         B = B,
                                         b = b),
                            control = list(maxiter = max_iter,
                                           tol = converge_crit))
      
    }
    
      return('NLS Model' = NLS_Cp_Model)

  }  
  
}

#================================================================================================================================================================
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# PK_NLS_Estimation.try
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#================================================================================================================================================================

{
  
  # This function pulls sepecific warnings or errors from the function
  # 'PK_NLS_Estimation'. The inputs are the same as those in 'PK_NLS_Estimation'
  # and returns the same output as 'PK_NLS_Estimation' if no errors or warnigns 
  # occur. Otherwise, "sigular.error", "max.iter.reached", or "other.error" is 
  # returned.
  
  PK_NLS_Estimation.try <- function(input_data, model, optim_method, max_iter, converge_crit){
    
    nls.fit <- tryCatch( 
      
      PK_NLS_Estimation(input_data = input_data,
                        model = model,
                        optim_method = optim_method,
                        max_iter = max_iter,
                        converge_crit = converge_crit)
      
      , 
      error   = function(e) e, 
      warning = function(w) w) 
    
    value <- toupper(as.character(nls.fit))
    
    chars.maxit <- "NUMBER OF ITERATIONS HAS REACHED `MAXITER'"
    chars.siglr <- "SINGULAR GRADIENT MATRIX AT INITIAL PARAMETER ESTIMATES"
    
    output <-      if(NROW(value) == 1 & grepl(chars.siglr, value[1], fixed = TRUE) == TRUE  & grepl(chars.maxit, value[1], fixed = TRUE) == FALSE){"sigular.error"}
              else if(NROW(value) == 1 & grepl(chars.siglr, value[1], fixed = TRUE) == FALSE & grepl(chars.maxit, value[1], fixed = TRUE) == TRUE){"max.iter.reached"}
              else if(NROW(value) == 1 & grepl(chars.siglr, value[1], fixed = TRUE) == FALSE & grepl(chars.maxit, value[1], fixed = TRUE) == FALSE){"other.error"}
              else   {nls.fit}
    
    return(output)
    
  }
  
}

#================================================================================================================================================================
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# PK_NLS_Estimation_Uniformitive
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#================================================================================================================================================================

{
  
  # This function is almost identical to 'PK_NLS_Estimation' with the exception 
  # that initial values for the NLS fitting are uninformative, i.e., "1" for all
  # parameters.
  #
  # Function inputs: input_data = the observed data containing at least a 
  #                               time column (labeled "time"), a plasma 
  #                               concentration column (labeled "conc.C_P"), 
  #                               and a column containing the initial dose 
  #                               (labeled "D_0")
  #                  model = specification for which model to fit to the 
  #                          data (i.e. C_1CM, L_1CM, E_1CM, R_1CM, C_2CM)
  #                  max_iter = the number of maximum iterations to allow
  #                             for in the NLS fitting
  #                  converge_crit = the convergence criteria for the NLS 
  #                                  fitting
  # 
  # Function outputs: "NLS Model" = the NLS model fit
  
  PK_NLS_Estimation_Uniformitive <- function(input_data,
                                             model,
                                             max_iter = 50,
                                             converge_crit = 1e-05){
    
    if (model=='C_1CM'){ # Constant-rate One-compartment model
      
      # grab initial dose from input_data
      
      D_0 <- min(input_data$D_0) # nmol
      
      #===========================================
      # Define function C_P,
      # the solution to the
      # Constant-rate One-compartment model
      #===========================================
      
      C_P <- function(t,V_P,k_E){
        D_0/V_P * exp(-k_E * t)
      }
      
      #===========================================
      # Apply NLS Model using `nlsLM`
      #===========================================
      
      NLS_Cp_Model <- nlsLM(conc.C_P ~ C_P(time,V_P,k_E),
                            data = input_data,
                            start = list(V_P = 1,
                                         k_E = 1),
                            control = list(maxiter = max_iter,
                                           tol = converge_crit))
      
    }
    else if (model=='L_1CM'){ # mLinear-rate One-compartment model
      
      # grab initial dose from input_data
      
      D_0 <- min(input_data$D_0) # nmol
      
      #===========================================
      # Define function C_P,
      # the solution to the
      # Linear-rate One-compartment model
      #===========================================
      
      C_P <- function(t,V_P,k_E0,b_E){
        D_0/V_P * exp( -( b_E/2*t + k_E0 ) * t )
      }
      
      #===========================================
      # Apply NLS Model using `nlsLM`
      #===========================================
      
      NLS_Cp_Model <- nlsLM(conc.C_P ~ C_P(time,V_P,k_E0,b_E),
                            data = input_data,
                            start = list(V_P = 1,
                                         k_E0 = 1,
                                         b_E  = 1),
                            control = list(maxiter = max_iter,
                                           tol = converge_crit))
      
      }
    else if (model=='E_1CM'){ # Exponential-rate One-compartment model
      
      # grab initial dose from input_data
      
      D_0 <- min(input_data$D_0) # nmol
      
      #===========================================
      # Define function C_P,
      # the solution to the
      # Exponential-rate One-compartment model
      #===========================================
      
      C_P <- function(t,V_P,k_E0,b_E){
        D_0/V_P * exp(  (k_E0/b_E) * (  exp(-b_E*t)  -  1  )  )
      }
      
      #===========================================
      # Apply NLS Model using `nlsLM`
      #===========================================
      
      NLS_Cp_Model <- nlsLM(conc.C_P ~ C_P(time,V_P,k_E0,b_E),
                            data = input_data,
                            start = list(V_P = 1,
                                         k_E0 = 1,
                                         b_E  = 1),
                            control = list(maxiter = max_iter,
                                           tol = converge_crit))
      
    }
    else if (model=='R_1CM'){ # Reciprocal-rate One-compartment model
      
      # grab initial dose from input_data
      
      D_0 <- min(input_data$D_0) # nmol
      
      #===========================================
      # Define function C_P,
      # the solution to the
      # Reciprocal-rate One-compartment model
      #===========================================
      
      C_P <- function(t,V_P,k_E0,b_E){
        D_0/V_P * exp(  -(k_E0/b_E) * log(  b_E*t  +  1  )  )
      }
      
      #===========================================
      # Apply NLS Model using `nlsLM`
      #===========================================
      
      NLS_Cp_Model <- nlsLM(conc.C_P ~ C_P(time,V_P,k_E0,b_E),
                            data = input_data,
                            start = list(V_P = 1,
                                         k_E0 = 1,
                                         b_E  = 1),
                            upper = c(Inf, Inf, Inf),
                            lower = c(-Inf, -Inf, -1/max(input_data$time)),
                            control = list(maxiter = max_iter,
                                           tol = converge_crit))
      
    }
    else if (model=='C_2CM'){ # Constant-rate Two-compartment model
      
      # grab initial dose from input_data
      
      D_0 <- min(input_data$D_0) # nmol
      
      #===========================================
      # Define function C_P,
      # the solution to the
      # Constant-rate Two-compartment model
      #===========================================
      
      C_P <- function(t,A,a,B,b){
        (A*exp(-(t*a))+B*exp(-(t*b)))
      }
      
      #===========================================
      # Apply NLS Model using `nlsLM`
      #===========================================
      
      NLS_Cp_Model <- nlsLM(conc.C_P ~ C_P(time,A,a,B,b),
                            data = input_data,
                            start = list(A = 1,
                                         a = 1,
                                         B = 1,
                                         b = 1),
                            control = list(maxiter = max_iter,
                                           tol = converge_crit))
      
    }
    
    
    return('NLS Model' = NLS_Cp_Model)
    
  }  
  
}

#================================================================================================================================================================
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# PK_NLS_Estimation_Uniformitive.try
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#================================================================================================================================================================

{
  
  # This function pulls sepecific warnings or errors from the function
  # 'PK_NLS_Estimation_Uniformitive'. The inputs are the same as those in 'PK_NLS_Estimation_Uniformitive'
  # and returns the same output as 'PK_NLS_Estimation_Uniformitive' if no errors or warnigns 
  # occur. Otherwise, "sigular.error", "max.iter.reached", or "other.error" is 
  # returned.
  
  PK_NLS_Estimation_Uniformitive.try <- function(input_data, model, max_iter, converge_crit){
    
    nls.fit <- tryCatch( 
      
      PK_NLS_Estimation_Uniformitive(input_data = input_data,
                                     model = model,
                                     max_iter = max_iter,
                                     converge_crit = converge_crit)
      
      , 
      error   = function(e) e, 
      warning = function(w) w) 
    
    value <- toupper(as.character(nls.fit))
    
    chars.maxit <- "NUMBER OF ITERATIONS HAS REACHED `MAXITER'"
    chars.siglr <- "SINGULAR GRADIENT MATRIX AT INITIAL PARAMETER ESTIMATES"
    
    output <-      if(NROW(value) == 1 & grepl(chars.siglr, value[1], fixed = TRUE) == TRUE  & grepl(chars.maxit, value[1], fixed = TRUE) == FALSE){"sigular.error"}
              else if(NROW(value) == 1 & grepl(chars.siglr, value[1], fixed = TRUE) == FALSE & grepl(chars.maxit, value[1], fixed = TRUE) == TRUE){"max.iter.reached"}
              else if(NROW(value) == 1 & grepl(chars.siglr, value[1], fixed = TRUE) == FALSE & grepl(chars.maxit, value[1], fixed = TRUE) == FALSE){"other.error"}
              else   {nls.fit}
    
    return(output)
    
  }
  
}

#================================================================================================================================================================
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# PK_NLS_Estimation.Seq
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#================================================================================================================================================================

{
  
  # This function applies the NLS estimation functions in a sequential order
  # until convergence (or until failure across all functions). First, NLS
  # estimation is conducted with 'Nelder-Mead' in optim. If this fails, the
  # estimation is run again but with 'SANN' instead. If this second function 
  # fails, a estimation is last attempted with uninformative initial 
  # parameters.
  #
  # The inputs for this functions are the same as PK_NLS_Estimation_Uniformitive.
  # The output returns the nls model fit (if converged) and a status on which 
  # estimation method convergence was attained.
  
  
  PK_NLS_Estimation.Seq <- function(input_data, model, max_iter, converge_crit){
    
    Status <- data.frame(c(NA,NA,NA))
    names(Status) <- NULL
    row.names(Status) <- c('LS/NM', 'SANN', 'Uninformitive')
    Status <- as.data.frame(Status)
    
    nls.fit <- PK_NLS_Estimation.try(input_data = input_data,
                                     model = model,
                                     optim_method = 'Nelder-Mead',
                                     max_iter = max_iter,
                                     converge_crit = converge_crit)
    
    nls.fit_1.NM <- nls.fit
    
    if(!(nls.fit[1] == "max.iter.reached" | nls.fit[1] == "sigular.error" | nls.fit[1] == "other.error")){
      
      Status[1,1] <- 'Converged'
      
    }else{
      
      Status[1,1] <- nls.fit[1]
      
    }
    # 
    # Status <- c()
    
    if(nls.fit[1] == "max.iter.reached" | nls.fit[1] == "sigular.error" | nls.fit[1] == "other.error"){
      
      # Status <- c(Status,nls.fit[1])
      
      nls.fit <- PK_NLS_Estimation.try(input_data = input_data,
                                       model = model,
                                       optim_method = 'SANN',
                                       max_iter = max_iter,
                                       converge_crit = converge_crit)
      
      nls.fit_2.SANN <- nls.fit
      
      if(nls.fit[1] == "max.iter.reached" | nls.fit[1] == "sigular.error" | nls.fit[1] == "other.error"){
        
        Status[2,1] <- nls.fit[1]
        
      }else{
        
        Status[2,1] <- 'Converged'
        
      }
      
    }
    
    
    
    
    if(nls.fit[1] == "max.iter.reached" | nls.fit[1] == "sigular.error" | nls.fit[1] == "other.error"){
      
      nls.fit <- PK_NLS_Estimation_Uniformitive.try(input_data = input_data,
                                                    model = model,
                                                    max_iter = max_iter,
                                                    converge_crit = converge_crit)
      
      nls.fit_3.Uniformitive <- nls.fit
      
      if(nls.fit[1] == "max.iter.reached" | nls.fit[1] == "sigular.error" | nls.fit[1] == "other.error"){
        
        Status[3,1] <- nls.fit[1]
        
      }else{
        
        Status[3,1] <- 'Converged'
        
      }
      
    }
    
    
    
    
    
    
    if(NROW(nls.fit) != 1){
      nls.fit$Status <- Status
      return(nls.fit)
    }else{
      names(Status) <- 'Status'
      return(Status)
    }
    
  }
  
}

#================================================================================================================================================================
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# C_P.fit
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#================================================================================================================================================================

{
  
  # This function calculates the estimated fitted function for each model.
  #
  # Function inputs: model = specification for which model to fit to the 
  #                          data (i.e. C_1CM, L_1CM, E_1CM, R_1CM, C_2CM)
  #                  time  = the observed time data (labeled "time"), a plasma 
  #                  D_0   = the initial dose
  #                  Other Variables = model parameters for related to the specified 'model'
  #
  # Function outputs: C_P = vector of fitted curve values.
  
  C_P.fit <- function(model, time, D_0 = NULL, V_P = NULL,
                      k_E = NULL, k_E0 = NULL, b_E = NULL,
                      A = NULL, a = NULL, B = NULL, b = NULL){
    
    if(model == 'C_1CM'){
      
      C_P <- D_0/V_P * exp(-k_E * time)
      
    }else if(model == 'L_1CM'){
      
      C_P <- D_0/V_P * exp( -( b_E/2*time + k_E0 ) * time )
      
    }else if(model == 'E_1CM'){
      
      C_P <- D_0/V_P * exp(  (k_E0/b_E) * (  exp(-b_E*time)  -  1  )  )
      
    }else if(model == 'R_1CM'){
      
      C_P <- D_0/V_P * exp(  -(k_E0/b_E) * log(  b_E*time  +  1  )  )
      
    }else if(model == 'C_2CM'){
      
      C_P <- A * exp(-time*a) + B * exp(-time*b)
      
    }
    
    return(C_P)
    
  }
  
}

#================================================================================================================================================================
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Approx.CI
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#================================================================================================================================================================

{
  
  # This function returns the estimates from the provided NLS model and calculates
  # returns approximate 95% confidence interval bound for each parameter in the
  # model.
  #
  # Function inputs: model = specification for which model to fit to the 
  #                          data (i.e. C_1CM, L_1CM, E_1CM, R_1CM, C_2CM)
  #                  D_0   = the initial dose
  #                  NLS_Cp_Model = NLS model fit object
  #
  # Function outputs: "Estimates and Approximate Confidence Bounds" 
  #                           = a data frame contain the NLS model estimates and
  #                             corresponding 95% confidence bounds.
  #                   "Estimates and Approximate Confidence Bounds - Untransformed
  #                           = (for C_2CM) a data frame contain the NLS model 
  #                             estimates and corresponding 95% confidence bounds
  #                             for the untranformed NLS model parameters.
  
  Approx.CI <- function(model, D_0 = NULL, NLS_Cp_Model){
    
    if(model == 'C_1CM'){
      
      Asym_CI <- nlConfint(NLS_Cp_Model, "b[1];b[2]")
      
      Asym_CI <- as.data.frame(round(Asym_CI[,c(2,1,3)],4))
      
      names(Asym_CI) <- c('95% Approx. Lower Bound','Estimates','95% Approx. Upper Bound')
      row.names(Asym_CI) <- c('V_P','k_E')
      
    }else if(model == 'C_2CM'){
      
      Asym_CI <- nlConfint(NLS_Cp_Model, "1/(b[1]+b[3]);
                                        (b[2])*(b[4])*(b[1]+b[3])/(b[1]*(b[4])+b[3]*(b[2]));
                                        b[1]*b[3]*((b[4])-(b[2]))^2/((b[1]+b[3])*(b[1]*(b[4])+b[3]*(b[2])));
                                        (b[1]*(b[4])+b[3]*(b[2]))/(b[1]+b[3])")
      
      Asym_CI <- as.data.frame(round(Asym_CI[,c(2,1,3)],4))
      Asym_CI[1,] <- Asym_CI[1,]*D_0
      
      names(Asym_CI) <- c('95% Approx. Lower Bound','Estimates','95% Approx. Upper Bound')
      row.names(Asym_CI) <- c('V_P','k_E','k_PD','k_DP')
      
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
      Asym_CI.Untransformed <- nlConfint(NLS_Cp_Model, "b[1]; b[2]; b[3]; b[4]")
      
      Asym_CI.Untransformed <- as.data.frame(round(Asym_CI.Untransformed[,c(2,1,3)],4))
      
      names(Asym_CI.Untransformed) <- c('95% Approx. Lower Bound','Estimates','95% Approx. Upper Bound')
      row.names(Asym_CI.Untransformed) <- c('A','a','B','b')
      
    }else{
      
      Asym_CI <- nlConfint(NLS_Cp_Model, "b[1];b[2];b[3]")
      
      Asym_CI <- as.data.frame(round(Asym_CI[,c(2,1,3)],4))
      
      names(Asym_CI) <- c('95% Approx. Lower Bound','Estimates','95% Approx. Upper Bound')
      row.names(Asym_CI) <- c('V_P','k_E0','b_E')
      
    }
    
    if(model=='C_2CM'){
      return(list(
        "Estimates and Approximate Confidence Bounds" = Asym_CI,
        "Estimates and Approximate Confidence Bounds - Untransformed" = Asym_CI.Untransformed
      ))
    }else{
      return(list(
        "Estimates and Approximate Confidence Bounds" = Asym_CI
      ))
    }
    
  }
  
}

#================================================================================================================================================================
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rat_Data_Load
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#================================================================================================================================================================

{
  
  # This script loads in the observed intravenous (IV) rat data into a 
  # list called 'rat_list'. This script truncates all observations prior
  # to the peak concentration value for each rat dataset and rescales all
  # time points such that the first datapoint is observed at time 0. 
  # This script requires "Rat_Data_Long_Format.csv" to live in the same
  # folder as this file. NOTE: One rat shares its peak concentration 
  # value at two consecutive time-points, resulting in a warning message.
  # This warning does not affect the intention of the truncation.
  
  # Function inputs: path.to.files = the file path to where Rat_Data_Long_Format.csv lives 
  #
  # Function outputs: A list containing the follwoing:
  #                           rat_list_full = a list of the data for all rats
  #                           rat_list      = a list of the data for the rats used in paper
  
  Rat_Data_Load <- function(path.to.files){
    
    # load in data
    rat.data <- read.csv(here(path.to.files,'Rat_Data_Long_Format.csv'))
    
    # identify unique rats
    unique.rats <- unique(paste0(rat.data$Sex,'.',rat.data$Rat))
    
    # create a data frame of the unique rats
    unique.rats.df <- data.frame(do.call('rbind', strsplit(unique.rats,'.',fixed=TRUE)))
    names(unique.rats.df) <- c('Sex', "Rat")
    
    # initialize the full rat list
    rat_list_full <- c()
    
    for(i in 1:NROW(unique.rats.df)){
      
      # grab specific rat
      temp.df <- rat.data[which(rat.data$Sex == unique.rats.df$Sex[i] & rat.data$Rat == unique.rats.df$Rat[i]),]
      # name columns
      names(temp.df) <- c("Time.hr", "conc.C_Pentration.muM", "Rat", "Sex", "Molecular.Weight", "Weight", "Dose.Amount")
      
      # put relevent information into temporary data frame
      temp.df <- data.frame(temp.df$Time.hr, 
                            temp.df$conc.C_Pentration.muM, 
                            temp.df$Weight*temp.df$Dose.Amount/temp.df$Molecular.Weight*1e9)
      # rename columns
      names(temp.df) <- c("time", "conc.C_P", "D_0")
      
      # removs data prior to peak concentration
      temp.df <- temp.df[-1:-(which(temp.df$conc.C_P==max(temp.df$conc.C_P))-1),]
      
      # reset row numbering
      rownames(temp.df) <- seq(length=nrow(temp.df))
      
      # shift time to start at 0
      temp.df$time <- temp.df$time - min(temp.df$time)
      
      # save data frame for specific rat
      assign(paste0(unique.rats.df$Sex[i],'.IV.',unique.rats.df$Rat[i]),
             temp.df)
      
      # put data frame into list
      temp.list <- list(get(paste0(unique.rats.df$Sex[i],'.IV.',unique.rats.df$Rat[i])))
      
      # name said list
      names(temp.list) <- paste0(unique.rats.df$Sex[i],'.IV.',unique.rats.df$Rat[i])
      
      # add lists together
      rat_list_full <- c(rat_list_full, temp.list)
      
    }
    
    # define function to truncate data after 2.2 hrs
    truncation.2.2hrs <- function(x){
      
      x <- x[x$time <= 2.2,]
      
      return(list(x))
      
    }
    
    # apply function to 'rat_list_full'
    rat_list <- sapply(rat_list_full, truncation.2.2hrs)
    
    # remove rats in which are not considered in the analsyis due to
    # have less than 1 hr of obserbations
    rat_list <- rat_list[!(names(rat_list) %in% c('Female.IV.Rat5',
                                                  'Male.IV.Rat6',
                                                  'Male.IV.Rat7',
                                                  'Male.IV.Rat9'))]
    
    return(list(rat_list_full = rat_list_full,
                rat_list = rat_list))
    
  }
  
}




