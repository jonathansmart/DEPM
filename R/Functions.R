#' Calculate Daily Egg Production (P0)
#' @description Daily egg production is estimated using mean egg density-at-age for each sample which is combined with a pre-defined
#'     estimate of egg mortality (Z). If multple years and/or regions are included in the data then a separate estimate of P0 can be
#'     provided by specifying the columns that correspond to these variables
#' @param data A data.frame that contains at least the following columns for Age,  density, and hatching time.
#'     The function can determine these columns based on similarly named variables
#' @param site A string containing the column name to idenify the site or sample
#' @param Z A numeric vector or single numeric value of egg mortality (Z) value to be used in P0 calculation
#' @param Region A string containing the column name for a region if desired as a grouping variable
#' @param Time A string containing the column name for a timestep if desired as a grouping variable
#'
#' @importFrom plyr ldply
#' @return A data.frame of P0 estimates and standard errors for each included Region and Time as well as the values of Z provided
#' @export

Estimate_P0 <- function(data, site, Z =  NULL, Region = NULL, Time = NULL){

  ## Error checks
  if(!is.numeric(Z) | is.null(Z)) stop("Egg mortality (Z) has not been correctly spcified")
  if(!is.null(Region) & !any(names(data) %in% Region)) stop("Region column could not be determined")
  if(!is.null(Time) & !any(names(data) %in% Time)) stop("Time column could not be determined")
  if(!any(names(data) %in% site)) stop("site column could not be determined")
  if(!any(names(data) %in% "Age")) stop("Age column could not be determined")
  if(any(is.na(data$Age))) stop("Age column contains NA's")

  # HAtch and density columns could be called different things. Therefore, use grep
  # to figure out which column name corresponds to these variables
  Hatch_col <- names(data)[grep("hatch", substr(tolower(names(data)),1,30))]
  Density_col <- names(data)[grep("dens", substr(tolower(names(data)),1,30))]

  if(length(Hatch_col) !=1) stop("Column for hatching time could not be distinguished")
  if(length(Density_col) !=1) stop("Column for Egg density could not be distinguished")

  # create a data subset to be used in the functions based on function arguments
  if(is.null(Time) & !is.null(Region)){
    processed_data <- dplyr::select(data,
                                    site = paste(site),
                                    Density = paste(Density_col),
                                    Age ,
                                    Hatch = paste(Hatch_col),
                                    Region = paste(Region))
  } else if(!is.null(Time) & is.null(Region)) {
    processed_data <- dplyr::select(data,
                                    site = paste(site),
                                    Density = paste(Density_col),
                                    Age ,
                                    Hatch = paste(Hatch_col),
                                    Time = paste(Time))
  }else if(!is.null(Time) & !is.null(Region)) {
    processed_data <- dplyr::select(data,
                                    site = paste(site),
                                    Density = paste(Density_col),
                                    Age ,
                                    Hatch = paste(Hatch_col),
                                    Region = paste(Region),
                                    Time = paste(Time))
  } else {
    processed_data <- dplyr::select(data,
                                    site = paste(site),
                                    Density = paste(Density_col),
                                    Age ,
                                    Hatch = paste(Hatch_col))
  }

  # P0 calculation function. This takes each of the users arguments and calculates P0
  # for a specified Z value. This is looped over before finishing Estimate_P0 to produce
  # estimates for multiple Z's
  calc_P0 <- function(data, site, Z =  NULL, Region = NULL, Time = NULL){

    Density_estimates <- data
    # Calculate egg density for each site based on age and hatch time
    Density_estimates$P <- Density_estimates$Density * exp(Z*Density_estimates$Age)/(Density_estimates$Hatch)

    ## sum egg density at each site
    # if 'Region' is specified, group_by Region
    if(is.null(Time) & !is.null(Region)){
      # if 'Region' is specified, group_by Region as well as site
      Density_estimates <- dplyr::group_by(Density_estimates, Region, site)
      # sum over egg densitities
      Density_estimates <- dplyr::summarise(Density_estimates, Pt = sum(P, na.rm = TRUE))
      Density_estimates <- dplyr::group_by(Density_estimates, Region)

    } else if(!is.null(Time) & is.null(Region)) {
      # if 'Time' is specified, group_by Time as well as site
      # rename the grouping variable for 'Time' so that it is recognised by dplyr::group_by
      Density_estimates <- dplyr::group_by(Density_estimates, Time, site)
      Density_estimates <- dplyr::summarise(Density_estimates, Pt = sum(P, na.rm = TRUE))
      Density_estimates <- dplyr::group_by(Density_estimates, Time)

    }else if(!is.null(Time) & !is.null(Region)) {
      # if 'Region' and 'Time are specified, group_by both as well as site
      Density_estimates <- dplyr::group_by(Density_estimates, Time, Region,  site)
      # sum over egg densitities
      Density_estimates <- dplyr::summarise(Density_estimates, Pt = sum(P, na.rm = TRUE))
      Density_estimates <- dplyr::group_by(Density_estimates, Time, Region)

    } else {
      # if neither 'Region' or 'Time are specified, only group by site
      Density_estimates <- dplyr::group_by(Density_estimates, site)
      # sum over egg densitities
      Density_estimates <- dplyr::summarise(Density_estimates, Pt = sum(P, na.rm = TRUE))
    }

    # Calculate the mean and se for P0 based on grouping variables.
    Density_estimates <- dplyr::summarise(Density_estimates, P0 = mean(Pt, na.rm = TRUE),
                                          P0_se = sqrt(var(Pt,na.rm = TRUE)/length(Pt)),
                                          Z = Z)
    return(Density_estimates)
  }

  results_list <- list()
  for(i in seq_along(Z)){
    results_list[[i]] <- calc_P0(processed_data, site = site, Z = Z[i], Region = Region, Time = Time)
  }

  results <-  plyr::ldply(results_list)
  return(results)
}

#' Generic ratio estimator
#' The internal function used to estimate spawning fraction and sex ratio
#' @param data A data.frame with numeric columns named "Top" and "Bottom" which represent the numerator and denominator, respectively
#' @param incl.extended.ratio.est A logical statement that determines whether the extended ratio method is also returned
#'
#' @return A data.frame that includes the ratio estimate, variance, standard error and coefficient of variation (CV). The
#'     extended ratio test is also included if requested.
#'
Estimate_Mean_Ratio <- function(data, incl.extended.ratio.est = FALSE){

  colnames(data) <- c("Top","Bottom")
  data <- na.omit(data)

  Top_mean <- mean(data$Top)
  Top_sd <- sd(data$Top)
  Bot_mean <- mean(data$Bottom)
  Bot_sd <- sd(data$Bottom)
  n <- length(data$Top)
  covar <- cov(data$Top, data$Bottom)

  rho <- covar/(Top_sd*Bot_sd)

  ratio_est <- Top_mean/Bot_mean
  extended_ratio_est <- ratio_est + (1/12)*(1/(Bot_mean^2))*
    (ratio_est*Bot_sd^2 - rho*Top_sd*Bot_sd)
  rel_diff <- ((extended_ratio_est-ratio_est)/ratio_est)*100
  ratio_var <- (1/n)*(1/Bot_mean^2)*(ratio_est^2*Bot_sd^2+Top_sd^2-2*ratio_est*covar)
  ratio_se <- sqrt(ratio_var)
  CV <- ratio_se/ratio_est


  if(incl.extended.ratio.est == TRUE){
    results <- data.frame(ratio_est,
                          extended_ratio_est,
                          rel_diff,
                          ratio_var,
                          ratio_se,
                          CV)

    colnames(results) <- c("Simple ratio estimate","Extended ratio estimate","Relative diff of extended and simple",
                           "Simple ratio var", "Simple ratio SE", "CV")
  } else{
    results <- data.frame(ratio_est,
                          ratio_var,
                          ratio_se,
                          CV)

    colnames(results) <- c("Ratio estimate", "Variance", "SE", "CV")

  }

  return(results)
}

#' Sex ratio estimator
#' @description Estimate sex ratio (R) of females for each sample in a survey. If multple years and/or regions are included in the data then a separate estimate
#'     of R can be provided by specifying the columns that correspond to these variables
#' @param data A data.frame that contains at least the following columns for number of females,  number of males, and total sample size.
#'     The function can determine these columns based on similarly named variables. Each line of the data.frame is assumed to be an indvivudal sample
#' @param Region A string containing the column name for a region if desired as a grouping variable
#' @param Time A string containing the column name for a timestep if desired as a grouping variable
#'
#' @return A data.frame that includes the female sex ratio estimate, variance, standard error and coefficient of variation (CV).
#' @export
#'
Estimate_sex_ratio <- function(data, Region = NULL, Time = NULL){
  if(!is.null(Region) & !any(names(data) %in% Region)) stop("Region column could not be determined")
  if(!is.null(Time) & !any(names(data) %in% Time)) stop("Time column could not be determined")

  data <- dplyr::ungroup(data)

  F_col <-  names(data)[grep("f", substr(tolower(names(data)),1,1))]
  Tot_col <-  names(data)[grep("tot", substr(tolower(names(data)),1,30))]

  if(is.null(Time) & !is.null(Region)){
    processed_data <- dplyr::select(data,
                                    Region = paste(Region),
                                    Top = paste(F_col),
                                    Bottom = paste(Tot_col))

    results <- expand.grid(Region = unique(processed_data$Region),
                           `Ratio estimate` = NA,
                           Variance = NA,
                           SE= NA,
                           CV= NA)

    for(i in unique(processed_data$Region)){
      tmp <- dplyr::filter(processed_data, Region == i)
      tmp <- dplyr::select(tmp,Top, Bottom)
      results[which(results$Region == i),c(2:5)]  <- Estimate_Mean_Ratio(tmp)
    }

  } else if(!is.null(Time) & is.null(Region)) {
    processed_data <- dplyr::select(data,
                                    Time = paste(Time),
                                    Top = paste(F_col),
                                    Bottom = paste(Tot_col))
    results <- expand.grid(Time = unique(processed_data$Time),
                           `Ratio estimate` = NA,
                           Variance = NA,
                           SE= NA,
                           CV= NA)

    for(i in unique(processed_data$Time)){
      tmp <- dplyr::filter(processed_data, Time == i)
      tmp <- dplyr::select(tmp,Top, Bottom)
      results[which(results$Time == i),c(2:5)]  <- Estimate_Mean_Ratio(tmp)
    }

  }else if(!is.null(Time) & !is.null(Region)) {

    processed_data <- dplyr::select(data,
                                    Region = paste(Region),
                                    Time = paste(Time),
                                    Top = paste(F_col),
                                    Bottom = paste(Tot_col))
    results <- expand.grid(Time = unique(processed_data$Time),
                           Region = unique(processed_data$Region),
                           `Ratio estimate` = NA,
                           Variance = NA,
                           SE= NA,
                           CV= NA)

    for(i in unique(processed_data$Time)){
      for(j in unique(processed_data$Region)){
        tmp <- dplyr::filter(processed_data, Time == i, Region == j)
        tmp <- dplyr::select(tmp,Top, Bottom)
        results[which(results$Time == i & results$Region == j),c(3:6)]  <- Estimate_Mean_Ratio(tmp)
      }
    }

  } else {
    processed_data <- dplyr::select(data,
                                    Top = paste(F_col),
                                    Bottom = paste(Tot_col))


    results <- Estimate_Mean_Ratio(processed_data)
  }

  return(results)

}

#' Spawning fraction estimator
#' @description Estimate Spawning fraction (S) of females for each sample in a survey. If multple years and/or regions are included in the data then a separate estimate
#'     of s can be provided by specifying the columns that correspond to these variables
#' @param data A data.frame that contains at least the following columns: "No","Yes" and "Total". These correspond to the numbers of individuals in each
#'     that were in spawning condition.
#'     The function can determine these columns based on similarly named variables. Each line of the data.frame is assumed to be an indvivudal sample
#' @param Region A string containing the column name for a region if desired as a grouping variable
#' @param Time A string containing the column name for a timestep if desired as a grouping variable
#'
#' @return A data.frame that includes the spawning fraction, variance, standard error and coefficient of variation (CV).
#' @export
#'
#'
Estimate_Spawning_fraction <- function(data, Region = NULL, Time = NULL){
  if(!is.null(Region) & !any(names(data) %in% Region)) stop("Region column could not be determined")
  if(!is.null(Time) & !any(names(data) %in% Time)) stop("Time column could not be determined")

  data <- dplyr::ungroup(data)

  Y_col <-  names(data)[grep("y", substr(tolower(names(data)),1,1))]

  if(length(Y_col)>1)
    Y_col <- Y_col[grep("yes",substr(tolower(Y_col),1,3))]

  if(length(Y_col)!=1)
    Y_col <-  names(data)[grep("spawn", substr(tolower(names(data)),1,30))]

  if(length(Y_col)>1) stop("Multiple columns for spawning fraction identified")

  Tot_col <-  names(data)[grep("tot", substr(tolower(names(data)),1,30))]

  if(is.null(Time) & !is.null(Region)){
    processed_data <- dplyr::select(data,
                                    Region = paste(Region),
                                    Top = paste(Y_col),
                                    Bottom = paste(Tot_col))

    results <- expand.grid(Region = unique(processed_data$Region),
                           `Ratio estimate` = NA,
                           Variance = NA,
                           SE= NA,
                           CV= NA)

    for(i in unique(processed_data$Region)){
      tmp <- dplyr::filter(processed_data, Region == i)
      tmp <- dplyr::select(tmp,Top, Bottom)
      results[which(results$Region == i),c(2:5)]  <- Estimate_Mean_Ratio(tmp)
    }


  } else if(!is.null(Time) & is.null(Region)) {
    processed_data <- dplyr::select(data,
                                    Time = paste(Time),
                                    Top = paste(Y_col),
                                    Bottom = paste(Tot_col))
    results <- expand.grid(Time = unique(processed_data$Time),
                           `Ratio estimate` = NA,
                           Variance = NA,
                           SE= NA,
                           CV= NA)

    for(i in unique(processed_data$Time)){
      tmp <- dplyr::filter(processed_data, Time == i)
      tmp <- dplyr::select(tmp,Top, Bottom)
      results[which(results$Time == i),c(2:5)]  <- Estimate_Mean_Ratio(tmp)
    }

  }else if(!is.null(Time) & !is.null(Region)) {

    processed_data <- dplyr::select(data,
                                    Region = paste(Region),
                                    Time = paste(Time),
                                    Top = paste(Y_col),
                                    Bottom = paste(Tot_col))
    results <- expand.grid(Time = unique(processed_data$Time),
                           Region = unique(processed_data$Region),
                           `Ratio estimate` = NA,
                           Variance = NA,
                           SE= NA,
                           CV= NA)

    for(i in unique(processed_data$Time)){
      for(j in unique(processed_data$Region)){
        tmp <- dplyr::filter(processed_data, Time == i, Region == j)
        tmp <- dplyr::select(tmp,Top, Bottom)
        results[which(results$Time == i & results$Region == j),c(3:6)]  <- Estimate_Mean_Ratio(tmp)
      }
    }

  } else {
    processed_data <- dplyr::select(data,
                                    Top = paste(Y_col),
                                    Bottom = paste(Tot_col))


    results <- Estimate_Mean_Ratio(processed_data)
  }

  return(results)

}



#' Estimate batch fecundity
#' @description Batch fecundity (number of eggs by female fish weight) is estimated with an allometric error structure allowing
#'     greater variance with increasing weight. This is estimated as a 4 parameter model where specified parameters can be fixed.
#' @param data A dataframe with 2 numeric variables: Fish weight and fecundity (number of eggs). Names and order are not important
#'    as the larger of the two variables is assumed to be number of eggs and is automatically assigned as this.
#' @param start_pars A list of 4 start_pars tha must include: "alpha", "beta", "Sigma0" and "Sigma1"
#' @param prediction.int A numeric vector of weights to predict batch fecundity over. Must be on the same scale as the data
#' @param fixed.pars A character vector of any start_pars that require fixing.
#' @param verbose If TRUE, parameter estimates are printed to the screen
#' @param return.parameters If TRUE, parameter estimates are returned instead of estimates
#' @useDynLib DEPM
#' @return Parameter estimates are automatically printed to the screen. If a prediction interval is provided then predicted fecundity-at-weight
#'     is provided at those intervals. If no prediction interval is provided than predictions for the raw data are returned. Variance (Var) is
#'     returned for both options. if `return.parameters == TRUE`, parameter estimates are returned instead of estimates
#' @import dplyr tidyr purrr stats
#' @export
#'
Estimate_Batch_Fecundity <- function(data, start_pars, prediction.int = NULL, return.parameters = FALSE,  fixed.pars= NULL, verbose = FALSE){

  if(any(prediction.int < 1) ) stop("Prediction intervals must be in grams not kilos")

  # suppressWarnings(dyn.load(dynlib('C:/UserData/Documents/SARDI/Modelling and programming/R/Packages/DEPM package/allometric_FvsW')))

  # function to get predictions from TMB
  create_TMB_sd_report_data.frame <- function(x){
    y <- unique(rownames(x))
    z <- NULL
    for(i in 1:length(y)){
      par <- data.frame(Parameter = y[i],
                        Val = x[rownames(x)==as.character(y[i]),1],
                        var = x[rownames(x)==as.character(y[i]),2])

      z <- rbind(z, par)
    }
    return(z)
  }

  if(ncol(data) != 2) stop("Two columns needed in data")

  if(!is.null(fixed.pars) &  any(!fixed.pars %in% names(start_pars)))
    stop("fixed.pars must be a vector of parameter names")

  if(any(!names(start_pars) %in% c("alpha","beta","Sigma0", "Sigma1") ))
    stop("start_pars must be a list containing starting values for the following parameters:
         alpha, beta, Sigma0 and Sigma1.\n Starting values must be provided even if a parameter is being fixed")

  if(mean(data[[1]]) < mean(data[[2]])){
    data <- list(x=data[[1]], y=data[[2]])
  } else{
    data <- list(x=data[[2]], y=data[[1]])
  }

  # alter TMB model if parameters need fixing.
  if(!is.null(fixed.pars)){
    # if there are parametes to be fixed.
    fixed_pars <- list(fixed.pars)

    fixed_pars <- vector("list", length = length(fixed.pars))
    names(fixed_pars) <- fixed.pars

    for(i in 1:length(fixed_pars)){
      fixed_pars[[i]] <- factor(NA)
    }
    model <- TMB::MakeADFun(data,start_pars,DLL="DEPM",
                       silent = TRUE,
                       checkParameterOrder=FALSE,
                       map = fixed_pars)
    alpha <- ifelse(!is.na(as.numeric(fit$par["alpha"])),as.numeric(fit$par["alpha"]),start_pars$alpha)
    beta <- ifelse(!is.na(as.numeric(fit$par["beta"])),as.numeric(fit$par["beta"]),start_pars$beta)
    Sigma0 <- ifelse(!is.na(as.numeric(fit$par["Sigma0"])),as.numeric(fit$par["Sigma0"]),start_pars$Sigma0)
    Sigma1 <- ifelse(!is.na(as.numeric(fit$par["Sigma1"])),as.numeric(fit$par["Sigma1"]),start_pars$Sigma1)

    if(return.parameters == TRUE)
      return(list(alpha = alpha, beta = beta, Sigma0 = Sigma0, Sigma1 = Sigma1))

  }else{
    model <- TMB::MakeADFun(data,start_pars,DLL="DEPM",checkParameterOrder=FALSE,silent = TRUE)
  }

  #fit model
  fit   <- nlminb(model$par, model$fn, model$gr)

  # results and covariance matrix
  rep   <- TMB::sdreport(model, getReportCovariance = T)



  # get the predictions and parameters with standard errors from TMB
  Derived_Quants <- create_TMB_sd_report_data.frame(summary(rep))
  if(return.parameters == TRUE){
    final_pars <- head(Derived_Quants, 4)
    return(final_pars)
  }

  alpha <- ifelse(!is.na(as.numeric(fit$par["alpha"])),as.numeric(fit$par["alpha"]),start_pars$alpha)
  beta <- ifelse(!is.na(as.numeric(fit$par["beta"])),as.numeric(fit$par["beta"]),start_pars$beta)
  Sigma0 <- ifelse(!is.na(as.numeric(fit$par["Sigma0"])),as.numeric(fit$par["Sigma0"]),start_pars$Sigma0)
  Sigma1 <- ifelse(!is.na(as.numeric(fit$par["Sigma1"])),as.numeric(fit$par["Sigma1"]),start_pars$Sigma1)



  if(verbose == TRUE)
    cat("alpha = ",alpha,"\n beta =", beta, "\n Sigma0 =", Sigma0, "\n Sigma1 =", Sigma1,"\n" )


  # If prediction data is provided then return that with the associated estimates
  if(!is.null(prediction.int)){


    results <- data.frame(Wt = prediction.int,
                          Fecundity = alpha * prediction.int^beta)
    results <-  dplyr::mutate(results, Var = (Sigma0*Fecundity^Sigma1)^2)

  } else{


    # Produce 95% CI's using the sigma parameters so that the variance is for the data rather than
    # the model
    results <- dplyr::filter(Derived_Quants,Parameter == "F_pred")
    results <- dplyr::bind_cols(Wt = data$x, Observed = data$y,results)
    results <- dplyr::select(results, -Parameter)
    results <- purrr::set_names(results, c("Wt", "Fecundity", "Predicted", "SD"))
    results <- dplyr::mutate(results, upp = Predicted +((Sigma0*Predicted^Sigma1)*1.96),
             low = Predicted -((Sigma0*Predicted^Sigma1)*1.96))
  }

  return(results)
}

#' Estimate_proportion_female
#' @description Determine the proportion of of the population in each weight bin. This is performed for females
#'     only and is an input parameter for the DEPMwt model.
#' @param data A dataframe that conatins female weight in grams. No males should be included and must be removed before use
#' @param Weight A character string of the variable in `data` that specifies weight in grams
#' @param max.weight An integer of the upper weight bin boundary in grams
#' @param bin.width An integer for the bin width in grams
#' @param Time A string containing the column name for a timestep if desired as a grouping variable
#' @param Region A string containing the column name for a region if desired as a grouping variable
#'
#' @return A dataframe with the columns Wt_bin, n, Prop and Prop_var. Timestep and region columns are returned if used
#' @export
#'
Estimate_proportion_female <- function(data, Weight, max.weight, bin.width, Time = NULL, Region = NULL){
  if(!is.null(Region) & !any(names(data) %in% Region)) stop("Region column could not be determined")
  if(!is.null(Time) & !any(names(data) %in% Time)) stop("Time column could not be determined")
  if(bin.width > max.weight) stop("Bin widths are larger than maximum weight")

  # Create weight bins based on the max weight for largest bin and their width
  wt_bins <- seq(0,max.weight, bin.width)

  # Group the data and perform analyses by Time, Region, Neither or both depending on
  # whether the variables were provided
  if(is.null(Time) & !is.null(Region)){
    # Subset necessary variables
    processed_data <- dplyr::select(data,
                                    Weight = paste(Weight),
                                    Region = paste(Region))
    # create weight bin integers
    processed_data <-  dplyr::mutate(processed_data,
                                     Wt_bin = cut(Weight,wt_bins,labels = FALSE))
    # Assign grouping variables
    processed_data <- dplyr::group_by(processed_data, Region, Wt_bin)
    # Get total number of females for each weight bin within each group
    processed_data <-  dplyr::summarise(processed_data, n = n())
    # calculate proportion females and the multionomial variance in each weigth bin
    processed_data <- dplyr::mutate(processed_data, Prop = n /sum(n),
                                    Prop_var = (Prop*(1-Prop))/sum(n))

    # At this point the data contains no zeroes and will be missing these weight bins.
    # create dummy data with all bins for each group and join together
    # This creates NA's which can be converted to zeroes
    processed_data <- suppressWarnings(
      suppressMessages(dplyr::right_join(processed_data,
                                         expand.grid(Region = unique(processed_data$Region),
                                                     Wt_bin = seq_along(wt_bins)))
      )
    )

  } else if(!is.null(Time) & is.null(Region)) {
    processed_data <- dplyr::select(data,
                                    Weight = paste(Weight),
                                    Time = paste(Time))

    processed_data <-  dplyr::mutate(processed_data,
                                     Wt_bin = cut(Weight,wt_bins,labels = FALSE))

    processed_data <- dplyr::group_by(processed_data, Time, Wt_bin)

    processed_data <-  dplyr::summarise(processed_data, n = n())

    processed_data <- dplyr::mutate(processed_data, Prop = n /sum(n),
                                    Prop_var = (Prop*(1-Prop))/sum(n))

    processed_data <-  suppressWarnings(
      suppressMessages(
        dplyr::right_join(processed_data,
                          expand.grid(Time = unique(processed_data$Time),
                                      Wt_bin = seq_along(wt_bins)))))


  }else if(!is.null(Time) & !is.null(Region)) {
    processed_data <- dplyr::select(data,
                                    Weight = paste(Weight),
                                    Region = paste(Region),
                                    Time = paste(Time))

    processed_data <-  dplyr::mutate(processed_data,
                                     Wt_bin = cut(Weight,wt_bins,labels = FALSE))

    processed_data <- dplyr::group_by(processed_data,Time, Region, Wt_bin)

    processed_data <-  dplyr::summarise(processed_data, n = n())

    processed_data <- dplyr::mutate(processed_data, Prop = n /sum(n),
                                    Prop_var = (Prop*(1-Prop))/sum(n))

    Bins <-  tidyr::expand(processed_data, tidyr::nesting(Time = Time, Region = Region), Wt_bin = seq_along(wt_bins))

    processed_data <-  suppressWarnings(
      suppressMessages(
        dplyr::right_join(processed_data,Bins)
      )
    )

  } else {
    processed_data <- dplyr::select(data,
                                    Weight = paste(Weight))

    processed_data <-  dplyr::mutate(processed_data,
                                     Wt_bin = cut(Weight,wt_bins,labels = FALSE))

    processed_data <- dplyr::group_by(processed_data, Wt_bin)

    processed_data <-  dplyr::summarise(processed_data, n = n())

    processed_data <- dplyr::mutate(processed_data, Prop = n /sum(n),
                                    Prop_var = (Prop*(1-Prop))/sum(n))

    processed_data <-  suppressWarnings(
      suppressMessages(
        dplyr::right_join(processed_data,
                          expand.grid( Wt_bin = seq_along(wt_bins)))
      )
    )


  }



  # Convert NA's to zeroes
  processed_data <-  dplyr::mutate_at(processed_data,vars(-group_cols()), funs(ifelse(is.na(.), 0, .)))

  return(processed_data)
}


#' Create a dataset of DEPM parameters to be passed to Estimate_biomass()
#'
#' @param P0 The results of Estimate_P0() with a single value of Z. Multiple values are not yet supported.
#'     Can be specified for each Time/Region combination if necessary
#' @param A A data.frame of spawning area in meters^2. Can be specified for each Time/Region combination if necessary
#' @param R The results of Estimate_sex_ratio(). Can be specified for each Time/Region combination if necessary
#' @param S The results of Estimate_Spawning_fraction(). Can be specified for each Time/Region combination if necessary
#' @param W_F The results of  Estimate_mean_W_F(). This is optional and only required if using the standard DEPM methods.
#'     Can be specified for each Time/Region combination if necessary
#'
#' @description Combine all of the estimates determined from DEPM parameter functions in this package. The structure of this object
#'    is compatible with the Estimate_biomass() function to facilitate its use. If multiple time steps/Regions are included then
#'    parameters specific to those surveys will be combined. Parameters that are used across years (for example a constant sex ratio)
#'    will automatically be included for each survey even if other parameters are being specified for specific steps/Regions.
#'    Mean Weight and mean fecundity are only required if using the standard DEPM approach that does not require weight bins.
#' @return A data.frame of DEPM parameters compatible with the Estimate_biomass() function
#' @export
#'
#'
combine_estimates <- function(P0, A, R, S, W_F = NULL){

  P0 <- dplyr::select(P0, -P0_se, -Z)
  results <- P0
  R <- dplyr::select(R, -Variance, -SE,-CV  )
  S <- dplyr::select(S, -Variance, -SE,-CV   )

  if(length(R)== 1) {
    results$R <- as.numeric(R)
  } else{
    results <- tryCatch(
      dplyr::left_join(results, dplyr::rename(R, R = "Ratio estimate")),
      error = function(e) stop("Inconsistent use of time and region")
    )
  }

  if(length(A)== 1) {
    results$A <- as.numeric(A)
  } else{
    results <- tryCatch(
      dplyr::left_join(results, A),
      error = function(e) stop("Inconsistent use of time and region")
    )
  }

  if(length(S)== 1) {
    results$S <- as.numeric(S)
  } else{
    results <- tryCatch(
      dplyr::left_join(results, dplyr::rename(S, S = "Ratio estimate")),
      error = function(e) stop("Inconsistent use of time and region")
    )
  }

  # mean weight and mean fecundity are optional.
  if(!is.null(W_F)){
    if(nrow(W_F)== 1) {
      results$W <- as.numeric(W_F$Mean_W)
      results$.F <- as.numeric(W_F$Mean_F)
    } else{
      W_F <- dplyr::select(W_F, -var_W, -var_F)
      results <- tryCatch(
        dplyr::left_join(results, dplyr::rename(W_F, W = "Mean_W", .F = "Mean_F")),
        error = function(e) stop("Inconsistent use of time and region")
      )
    }
  }


  return(results)

}



#' Create a dataset of DEPM parameter variances to be passed to Estimate_biomass()
#'
#' @description Combine all of the variances determined from DEPM parameter functions in this package. The structure of this object
#'    is compatible with the Estimate_biomass() function to facilitate its use. If multiple time steps/Regions are included then
#'    parameter variances specific to those surveys will be combined. Parameters that are used across years (for example a constant sex ratio)
#'    will automatically be included for each survey even if other parameters are being specified for specific steps/Regions.
#' @param P0 The results of Estimate_P0() with a single value of Z. Multiple values are not yet supported.
#' @param R A data.frame of sex ratios. Can be specified for each Time/Region combination if necessary
#' @param S A data.frame of spawning fractions. Can be specified for each Time/Region combination if necessary
#' @param W_F The results of  Estimate_mean_W_F(). This is optional and only required if using the standard DEPM methods.
#'     Can be specified for each Time/Region combination if necessary
#'
#' @return A data.frame of DEPM parameters compatible with the Estimate_biomass() function
#' @export
#'
#'
#'
combine_variances <- function(P0,R, S, W_F = NULL){

  P0 <- dplyr::select(P0, -P0, -Z)
  results <-  dplyr::rename(P0, P0 = P0_se)
  R <- dplyr::select(R, -`Ratio estimate`, -SE,-CV  )
  S <- dplyr::select(S, -`Ratio estimate`, -SE,-CV   )

  if(length(R)== 1) {
    results$R <- as.numeric(R)
  } else{
    results <- tryCatch(
      dplyr::left_join(results, dplyr::rename(R, R = "Variance")),
      error = function(e) stop("Inconsistent use of time and region")
    )
  }

  if(length(S)== 1) {
    results$S <- as.numeric(S)
  } else{
    results <- tryCatch(
      dplyr::left_join(results, dplyr::rename(S, S = "Variance")),
      error = function(e) stop("Inconsistent use of time and region")
    )
  }

  # mean weight and mean fecundity are optional.
  if(!is.null(W_F)){
    if(nrow(W_F)== 1) {
      results$W <- as.numeric(W_F$var_W)
      results$.F <- as.numeric(W_F$var_F)
    } else{
      W_F <- dplyr::select(W_F, -Mean_W, -Mean_F)
      results <- tryCatch(
        dplyr::left_join(results, dplyr::rename(W_F, W = "var_W", .F = "var_F")),
        error = function(e) stop("Inconsistent use of time and region")
      )
    }
  }


  return(results)

}

#' Create a dataset of DEPM parameter weight class parameter estimates and variances to be passed to Estimate_biomass()
#' @description Combine the estimates and variances determined from the Estimate_proportion_female() and Estimate_Batch_Fecundity()
#'     functions in this package. The structure of this object is compatible with the Estimate_biomass() function to facilitate
#'     its use. If multiple time steps/Regions are included then parameter estimates and variances specific to those surveys will
#'     be combined.
#' @param prop.fem.data A dataframe produced from the the Estimate_proportion_female() function
#' @param fecundity.data A dataframe produced from the the Estimate_Batch_Fecundity()  function
#'
#' @return A dataframe with the columns Wt_bin, n, Prop, Prop_var, Wt, Fecundity and Fec_var.
#'     Time and Region are included if included in either dataset
#' @export
#'
combine_wt_class_estimates <- function(prop.fem.data, fecundity.data){

  # Check that the proportion females data have columns with "prop" in them
  prop_fem_check <- length(names(prop.fem.data)[grep("prop", substr(tolower(names(prop.fem.data)),1,30))])
  if(prop_fem_check== 0) stop("prop.fem.data does not appear to have been created using `Estimate_proportion_female`")

  # Check that the proportion females data have columns with "fecundity" in them
  fec_check <- length(names(fecundity.data)[grep("fecundity", substr(tolower(names(fecundity.data)),1,30))])
  if(fec_check== 0) stop("fecundity.data does not appear to have been created using `Estimate_Batch_Fecundity`")

  ## Process the fecundity data by
  # convert the var name to fec_var
  process_fec_data <- dplyr::rename(fecundity.data, Fec_var = Var)
  # make sure the Wt variable is in ascending order
  process_fec_data <- dplyr::arrange(process_fec_data, Wt)
  # Use the row number to get an integer weight bin
  process_fec_data <-  dplyr:: mutate(process_fec_data, Wt_bin = row_number())

  # There should be one more weight bin in the female prop data than the fecundity data
  if(max(process_fec_data$Wt_bin) != max(prop.fem.data$Wt_bin)-1) stop("Weight bins in data sets do not align")

  # create a vector of joining variables depending on whether Time or Region were specified
  if(all((any(names(prop.fem.data) %in% "Region") &  any(names(process_fec_data) %in% "Region")),
         (any(names(prop.fem.data) %in% "Time") &  any(names(process_fec_data) %in% "Time")))){
    join <- c("Time","Region", "Wt_bin")
  }else if(any(names(prop.fem.data) %in% "Time") &  any(names(process_fec_data) %in% "Time")){
    join <- c("Time", "Wt_bin")
  } else if(any(names(prop.fem.data) %in% "Region") &  any(names(process_fec_data) %in% "Region")){
    join <- c("Region", "Wt_bin")
  } else{
    join <- "Wt_bin"
  }

  # join data sets togther
  combined_data <-  dplyr::left_join(prop.fem.data, process_fec_data, by = join)
  # NA's are produced when weight bins don't correspond. get rid of these
  combined_data <- dplyr::filter(combined_data, !is.na(Wt))

  return(combined_data)
}

#' Combined DEPM parameters to calculate spawning biomass, total number of females and the number of females per weight class using the DEPMwt approach
#'
#' @param adult.pars The output of `combine_estimates`. This included each of the DEPM parameters estimated using this package.
#'     These can be grouped by Timestep/Region and will be used to estimate results for each combination of these
#' @param adult.vars The output of `combine_variances`. This included each of the DEPM parameters estimated using this package.
#'     These can be grouped by Timestep/Region and will be used to estimate results for each combination of these
#' @param weight.pars.vars  A dataframe that must include the columns: "Wt_bin", "Prop", "Prop_var", "Wt", "Fecundity", "Fec_var". These correspond
#'     to weight bin number (must start at 1), proportion of females in each weight bin, the multinomial variance of the proportion of females
#'     in each weight bin, weight of each weight bin, fecundity-at-weight of that weight bin and the variance of fecundity-at-weight for that
#'     weight bin, respectively. Wt, fecundity and Fec_var can be calculated using the `Estimate_batch_fecundity` function. These can be
#'     grouped by Timestep/Region and will be used to estimate results for each combination of these
#'
#' @return A list of results for spawning biomass, total number of females and the number of females per weight class. These list elements
#'     will be grouped by Timestep/Region if these were provided with the parameters. The standard deviations of these estimates are also returned
#' @export
#'
Estimate_DEPMWt_biomass<- function(adult.pars, adult.vars, weight.pars.vars){

  if(!any(names(adult.pars) %in% names(adult.vars))) stop("Different Time/Regions provided between adult parameter and variance objects")

  if(any(names(adult.pars)== "Time") & !any(names(weight.pars.vars) == "Time"))
    stop("Different Timesteps provided between adult parameter objects and weight parameter objects")
  if(any(names(adult.pars)== "Region") & !any(names(weight.pars.vars) == "Region"))
    stop("Different Regions provided between adult parameter objects and weight parameter objects")



  ## Underlying Biomass and variance estimation functions by Rick ------------------------

  BspNfemNfemWt = function(P0, A, R, S, Fw, PNw, Ww, nw){


    # if(nw==1){ #For this case of a single wt class, traditional DEPM, uisng mean F and W, applies. (R McGarvey 1Aug19)
    #   print("You have entered nw=1.  The Bsp and Nfem estimates will use the standard Parker equation where it is assumed that the function inputs for Fw=batch fecundity and Ww=weight of females are the means")
    #   print("The sole value for PNw must equal 1.")
    #   if(nw==1 & PNw !=1) stop("The sole value for PNw must equal 1.")
    #   #If these conditions (means for F and W, PNw=1)  apply, the equations for Nfem and Bsp as coded below will work as written in the traditional Parker equation.
    # }
      #Check inputs
    if(length(Ww) != nw) stop("length(Ww) != nw")
    if(length(Fw) != nw) stop("length(Fw) != nw")
    if(length(PNw) != nw) stop("length(PNw) != nw")
    if(sum(PNw) != 1) stop("sum(PNw) != 1")

    #Compute Nfem, Bsp, and NfemWt
    Nfem = P0 * A / ( S * sum(Fw * PNw) )
    Bsp = Nfem * sum(PNw * Ww) / R
    NfemWt = Nfem * PNw
    return( list(Nfem=Nfem, Bsp=Bsp, NfemWt=NfemWt ) )
    #was:  return( data.frame(Nfem=Nfem, Bsp=Bsp, NfemWt1=NfemWt[1],NfemWt2=NfemWt[2], NfemWt3=NfemWt[3],NfemWt4=NfemWt[4], NfemWt5=NfemWt[5],NfemWt6=NfemWt[6], NfemWt7=NfemWt[7],NfemWt8=NfemWt[8], NfemWt9=NfemWt[9], NfemWt10=NfemWt[10], NfemWt11=NfemWt[11],NfemWt12=NfemWt[12], NfemWt13=NfemWt[13],NfemWt14=NfemWt[14], NfemWt15=NfemWt[15],NfemWt16=NfemWt[16], NfemWt17=NfemWt[17],NfemWt18=NfemWt[18], NfemWt19=NfemWt[19],NfemWt20=NfemWt[20], NfemWt21=NfemWt[21],NfemWt22=NfemWt[22], NfemWt23=NfemWt[23],NfemWt24=NfemWt[24], NfemWt25=NfemWt[25],NfemWt26=NfemWt[26],NfemWt27=NfemWt[27]) )
  }


  VarBsp = function(P0, A, R, S, Fw, PNw, VP0, VA, VR, VS, VFw, VPNw, Ww, nw ){
    #Check that the vectors are all of the correct length:
    if(length(Fw) != nw) stop("length(Fw) does not = nw")
    if(length(PNw) != nw) stop("length(PNw) does not = nw")
    if(length(VFw) != nw) stop("length(VFw) does not = nw")
    if(length(VPNw) != nw) stop("length(VPNw) does not = nw")
    if(length(Ww) != nw) stop("length(Ww) does not = nw")

    F1=1/(R^4 * S^4 * (sum(Fw*PNw))^4 )
    S1 = P0^2 * A^2 * R^2 * VS +  P0^2 * A^2 * VR * S^2 +  P0^2 * VA * R^2 * S^2 +  VP0 * A^2 * R^2 * S^2
    S21sq = ( sum(Fw*PNw) )^2
    S22sq = ( sum(PNw*Ww) )^2
    F3 = P0^2 * A^2 * R^2 * S^2
    S31 = sum( VFw* (PNw^2) )
    S32sq = S22sq
    F4=F3
    S4Outer=0.
    for (i in 1:nw){
      wvec=1:nw
      wNoti = wvec[-i]
      S41=0.
      for (w in wNoti) {
        S41 = S41 + Fw[w]*PNw[w]
      }
      S42=0.
      for (w in wNoti) {
        S42 = S42 + PNw[w]*Ww[w]
      }
      S4Innersqi = ( Ww[i]*S41 - Fw[i]*S42 )^2
      S4Outer = S4Outer + VPNw[i] * S4Innersqi
    }  #end loop over i

    VarBsp = F1 * (S1*S21sq*S22sq + F3*S31*S32sq + F4*S4Outer)
    return(VarBsp)

  }


  VarNfem = function(P0, A, S, Fw, PNw, VP0, VA, VS, VFw, VPNw, nw ){
    #Check that the vectors are all of the correct length:
    if(length(Fw) != nw) stop("length(Fw) does not = nw")
    if(length(PNw) != nw) stop("length(PNw) does not = nw")
    if(length(VFw) != nw) stop("length(VFw) does not = nw")
    if(length(VPNw) != nw) stop("length(VPNw) does not = nw")

    F1=1/(S^4 * (sum(Fw*PNw))^4 )
    S2 = P0^2 * A^2 * VS  +  P0^2 * VA * S^2 +  VP0 * A^2 *  S^2
    S21sq = ( sum(Fw*PNw) )^2
    F3 = P0^2 * A^2 *S^2
    S3 = sum( VFw* (PNw^2) )
    S4 = sum( VPNw* (Fw^2) )
    F4=F3

    VarNfem = F1 * (S2*S21sq + F3*S3 + F4*S4)
    return(VarNfem)

  }

  VarNfemWt = function(P0, A, S, Fw, PNw, VP0, VA, VS, VFw, VPNw, nw ){
    #Check that the vectors are all of the correct length:
    if(length(Fw) != nw) stop("length(Fw) does not = nw")
    if(length(PNw) != nw) stop("length(PNw) does not = nw")
    if(length(VFw) != nw) stop("length(VFw) does not = nw")
    if(length(VPNw) != nw) stop("length(VPNw) does not = nw")

    #Some quantities vary with w (the specific wt class var being computed, and some will be the same for all w.
    #We compute these 5 quantities prior to entering the sum over w
    F1=1/(S^4 * (sum(Fw*PNw))^4 )
    S21 = P0^2 * A^2 * VS  +  P0^2 * VA * S^2 +  VP0 * A^2 *  S^2
    S22sq = ( sum(Fw*PNw) )^2
    F3 = P0^2 * A^2 * S^2
    F4=F3

    #Define vector into which we will assign the final outputted VarNfemWt vector
    VarNfemWt = rep(NA, nw)  #nw=9
    for (w in 1:nw){
      F21 = PNw[w]^2
      F31 = VFw[w] *  (PNw[w])^4
      F32 = VPNw[w]
      F41 = F21

      #Now we are looping over all wprime's except w in the outer loop
      wvec=1:nw
      wpNotw = wvec[-w]
      S33=0.
      for (wp in wpNotw) {
        S33 = S33 + Fw[wp]*PNw[wp]
      }
      S33sq = S33^2
      S42=0.
      for (wp in wpNotw) {
        S42 = S42 + VFw[wp]*(PNw[wp])^2
      }
      S43=0.
      for (wp in wpNotw) {
        S43 = S43 + VPNw[wp]*(Fw[wp])^2
      }


      VarNfemWt[w] = F1 * (F21*S21*S22sq + F3* (F31 + F32 * S33sq) + F4 * F41 * (S42+S43) )

    } #end loop over w

    return(VarNfemWt)

  }


  if(any(names(adult.pars) %in% "Region" ) & any(names(adult.pars) %in% "Time")){
    # List of results to be filled: Number of females, Biomass and number of females at weight class
    resultlist <- list(
      Nfem =  tidyr::expand(adult.pars, tidyr::nesting(Time = Time, Region = Region),  Nfem = NA, SD = NA ),
      Biomass = tidyr::expand(adult.pars, tidyr::nesting(Time = Time, Region = Region),  Biomass = NA, SD = NA ),
      NfemWt = tidyr::expand(weight.pars.vars, tidyr::nesting(Time = Time, Region = Region, Nwt = Wt_bin), NfemWt = NA , SD = NA)
      )

    # Loop over Time and Region and save estimates and variances to correct position in resultlist
    for(i in unique(adult.pars$Time)){
      for(j in unique(adult.pars$Region)){

        tmp <-  dplyr::filter(adult.pars, Time == i, Region == j)

        if(nrow(tmp) == 0)next

        tmp_var <-  dplyr::filter(adult.vars, Time == i, Region == j)

        wt_tmp <-  dplyr::filter(weight.pars.vars,Time == i, Region == j, !is.na(Wt))

        res <- BspNfemNfemWt(P0 = tmp$P0, A = tmp$A, R = tmp$R, S = tmp$S,
                             Ww = wt_tmp$Wt, Fw = wt_tmp$Fecundity, PNw = wt_tmp$Prop, nw = max(wt_tmp$Wt_bin))

        bspvar <- sqrt(VarBsp(P0 = tmp$P0, VP0 = tmp_var$P0^2,
                              A = tmp$A, VA = 0,
                              R = tmp$R,  VR = tmp_var$R,
                              S = tmp$S,  VS = tmp_var$S,
                              Ww = wt_tmp$Wt,
                              nw = max(wt_tmp$Wt_bin),
                              Fw = wt_tmp$Fecundity, VFw = wt_tmp$Fec_var,
                              PNw = wt_tmp$Prop, VPNw = wt_tmp$Prop_var ))/1000


        NfemVar <- sqrt(VarNfem(P0 = tmp$P0, VP0 = tmp_var$P0^2,
                                A = tmp$A, VA = 0,
                                S = tmp$S,  VS = tmp_var$S,
                                nw = max(wt_tmp$Wt_bin),
                                Fw = wt_tmp$Fecundity, VFw = wt_tmp$Fec_var,
                                PNw = wt_tmp$Prop, VPNw = wt_tmp$Prop_var ))


        NfemWtVar <- sqrt(VarNfemWt(P0 = tmp$P0, VP0 = tmp_var$P0^2,
                                    A = tmp$A, VA = 0,
                                    S = tmp$S,  VS = tmp_var$S,
                                    nw = max(wt_tmp$Wt_bin),
                                    Fw = wt_tmp$Fecundity, VFw = wt_tmp$Fec_var,
                                    PNw = wt_tmp$Prop, VPNw = wt_tmp$Prop_var ))


        resultlist[["Nfem"]][resultlist[["Nfem"]]$Time == i & resultlist[["Nfem"]]$Region == j,
                             "Nfem"] <- res$Nfem

        resultlist[["Biomass"]][resultlist[["Biomass"]]$Time == i & resultlist[["Biomass"]]$Region == j,
                                "Biomass"] <- res$Bsp/1000

        resultlist[["NfemWt"]][resultlist[["NfemWt"]]$Time == i & resultlist[["NfemWt"]]$Region == j,
                               "NfemWt" ] <- res$NfemWt

        resultlist[["Nfem"]][resultlist[["Nfem"]]$Time == i & resultlist[["Nfem"]]$Region == j,
                             "SD"] <- NfemVar

        resultlist[["Biomass"]][resultlist[["Biomass"]]$Time == i & resultlist[["Biomass"]]$Region == j,
                                "SD"] <- bspvar

        resultlist[["NfemWt"]][resultlist[["NfemWt"]]$Time == i & resultlist[["NfemWt"]]$Region == j,
                               "SD" ] <- NfemWtVar

      }
    }

  } else if(any(names(adult.pars) %in% "Time")){
    resultlist <- list(
      Nfem = expand.grid(Time = unique(adult.pars$Time), Nfem = NA, SD = NA ),
      Biomass = expand.grid(Time = unique(adult.pars$Time),  Biomass = NA , SD = NA),
      NfemWt = expand.grid(Time = unique(adult.pars$Time),
                           Nwt = seq(1, length(unique(weight.pars.vars$Wt_bin)),1), NfemWt = NA , SD = NA)
    )

    for(i in unique(adult.pars$Time)){


      tmp <-  dplyr::filter(adult.pars,Time == i)

      if(nrow(tmp) == 0)next

      tmp_var <- dplyr::filter(adult.vars,Time == i)

      wt_tmp <- dplyr::filter(weight.pars.vars,Time == i , !is.na(Wt))

      res <- BspNfemNfemWt(P0 = tmp$P0, A = tmp$A, R = tmp$R, S = tmp$S,
                           Ww = wt_tmp$Wt, Fw = wt_tmp$Fecundity, PNw = wt_tmp$Prop, nw = max(wt_tmp$Wt_bin))

      bspvar <- sqrt(VarBsp(P0 = tmp$P0, VP0 = tmp_var$P0^2,
                            A = tmp$A, VA = 0,
                            R = tmp$R,  VR = tmp_var$R,
                            S = tmp$S,  VS = tmp_var$S,
                            Ww = wt_tmp$Wt,
                            nw = max(wt_tmp$Wt_bin),
                            Fw = wt_tmp$Fecundity, VFw = wt_tmp$Fec_var,
                            PNw = wt_tmp$Prop, VPNw = wt_tmp$Prop_var ))/1000


      NfemVar <- sqrt(VarNfem(P0 = tmp$P0, VP0 = tmp_var$P0^2,
                              A = tmp$A, VA = 0,
                              S = tmp$S,  VS = tmp_var$S,
                              nw = max(wt_tmp$Wt_bin),
                              Fw = wt_tmp$Fecundity, VFw = wt_tmp$Fec_var,
                              PNw = wt_tmp$Prop, VPNw = wt_tmp$Prop_var ))


      NfemWtVar <- sqrt(VarNfemWt(P0 = tmp$P0, VP0 = tmp_var$P0^2,
                                  A = tmp$A, VA = 0,
                                  S = tmp$S,  VS = tmp_var$S,
                                  nw = max(wt_tmp$Wt_bin),
                                  Fw = wt_tmp$Fecundity, VFw = wt_tmp$Fec_var,
                                  PNw = wt_tmp$Prop, VPNw = wt_tmp$Prop_var ))


      resultlist[["Nfem"]][resultlist[["Nfem"]]$Time == i, "Nfem"] <- res$Nfem

      resultlist[["Biomass"]][resultlist[["Biomass"]]$Time == i ,"Biomass"] <- res$Bsp/1000

      resultlist[["NfemWt"]][resultlist[["NfemWt"]]$Time == i,"NfemWt" ] <- res$NfemWt

      resultlist[["Nfem"]][resultlist[["Nfem"]]$Time == i ,"SD"] <- NfemVar

      resultlist[["Biomass"]][resultlist[["Biomass"]]$Time == i,"SD"] <- bspvar

      resultlist[["NfemWt"]][resultlist[["NfemWt"]]$Time == i, "SD" ] <- NfemWtVar

    }


  } else if(any(names(adult.pars) %in% "Region")){
    resultlist <- list(
      Nfem = expand.grid(Region = unique(adult.pars$Region), Nfem = NA, SD = NA ),
      Biomass = expand.grid(Region = unique(adult.pars$Region),  Biomass = NA , SD = NA),
      NfemWt = expand.grid(Region = unique(adult.pars$Region),
                           Nwt = seq(1, length(unique(weight.pars.vars$Wt_bin)),1), NfemWt = NA , SD = NA)
    )


    for(j in unique(adult.pars$Region)){

      tmp<- dplyr::filter(adult.pars, Region == j)

      if(nrow(tmp) == 0)next

      tmp_var <- dplyr::filter(adult.vars, Region == j)

      wt_tmp <-  dplyr::filter(weight.pars.vars, Region == j, !is.na(Wt))

      res <- BspNfemNfemWt(P0 = tmp$P0, A = tmp$A, R = tmp$R, S = tmp$S,
                           Ww = wt_tmp$Wt, Fw = wt_tmp$Fecundity, PNw = wt_tmp$Prop, nw = max(wt_tmp$Wt_bin))

      bspvar <- sqrt(VarBsp(P0 = tmp$P0, VP0 = tmp_var$P0^2,
                            A = tmp$A, VA = 0,
                            R = tmp$R,  VR = tmp_var$R,
                            S = tmp$S,  VS = tmp_var$S,
                            Ww = wt_tmp$Wt,
                            nw = max(wt_tmp$Wt_bin),
                            Fw = wt_tmp$Fecundity, VFw = wt_tmp$Fec_var,
                            PNw = wt_tmp$Prop, VPNw = wt_tmp$Prop_var ))/1000


      NfemVar <- sqrt(VarNfem(P0 = tmp$P0, VP0 = tmp_var$P0^2,
                              A = tmp$A, VA = 0,
                              S = tmp$S,  VS = tmp_var$S,
                              nw = max(wt_tmp$Wt_bin),
                              Fw = wt_tmp$Fecundity, VFw = wt_tmp$Fec_var,
                              PNw = wt_tmp$Prop, VPNw = wt_tmp$Prop_var ))


      NfemWtVar <- sqrt(VarNfemWt(P0 = tmp$P0, VP0 = tmp_var$P0^2,
                                  A = tmp$A, VA = 0,
                                  S = tmp$S,  VS = tmp_var$S,
                                  nw = max(wt_tmp$Wt_bin),
                                  Fw = wt_tmp$Fecundity, VFw = wt_tmp$Fec_var,
                                  PNw = wt_tmp$Prop, VPNw = wt_tmp$Prop_var ))


      resultlist[["Nfem"]][resultlist[["Nfem"]]$Region == j,"Nfem"] <- res$Nfem

      resultlist[["Biomass"]][resultlist[["Biomass"]]$Region == j,"Biomass"] <- res$Bsp/1000

      resultlist[["NfemWt"]][ resultlist[["NfemWt"]]$Region == j, "NfemWt" ] <- res$NfemWt

      resultlist[["Nfem"]][resultlist[["Nfem"]]$Region == j,"SD"] <- NfemVar

      resultlist[["Biomass"]][ resultlist[["Biomass"]]$Region == j,"SD"] <- bspvar

      resultlist[["NfemWt"]][resultlist[["NfemWt"]]$Region == j, "SD" ] <- NfemWtVar

    }


  } else {
    stop("Could not determine Time and Region columns")
  }

  return(resultlist)
}


#' Combined DEPM parameters to calculate spawning biomass, total number of females and the number of females per weight class using the standard
#' DEPM approach based off the Parker equation.
#'
#' @param adult.pars The output of `combine_estimates` with mean Weight and mean F included. This included each of the DEPM parameters estimated using this package.
#'     These can be grouped by Timestep/Region and will be used to estimate results for each combination of these
#' @param adult.vars The output of `combine_variances`with var Weight and var F inlcuded. This included each of the DEPM parameters estimated using this package.
#'     These can be grouped by Timestep/Region and will be used to estimate results for each combination of these
#'
#' @return A list of results for spawning biomass and total number of females. These list elements
#'     will be grouped by Timestep/Region if these were provided with the parameters. The standard deviations of these estimates are also returned
#' @export
#'
Estimate_DEPM_Biomass <- function(adult.pars, adult.vars){

  if(!any(names(adult.pars) %in% names(adult.vars))) stop("Different Time/Regions provided between adult parameter and variance objects")

    ## Underlying Biomass and variance estimation functions by Rick ------------------------

  BspNfem = function(P0, A, R, S, Fw, Ww){

    #Compute Nfem, Bsp, and NfemWt
    Nfem = P0 * A / ( S * Fw)
    Bsp = Nfem * Ww / R


    #Return outputs as a list
    return( list(Nfem=Nfem, Bsp=Bsp ) )
  }



  VarBspNfemnw1 = function(P0, A, R, S, F1, W1, VP0, VA, VR, VS, VF1, VW1 ){


    #Compute Nfem, Bsp (but not NfemWt, which is only computed for the nw>1 case).
    Nfem = P0 * A / ( S * F1 )
    Bsp = Nfem * W1 / R

    #Compute all the CV's squared:
    CV2P0 = VP0 / (P0^2)
    CV2A = VA / (A^2)
    CV2R = VR / (R^2)
    CV2S = VS / (S^2)
    CV2F1 = VF1 / (F1^2)
    CV2W1 = VW1 / (W1^2)


    #Compute variances for Bsp and Nfem
    VarNfem = Nfem^2 * (CV2P0 + CV2A + CV2S +CV2F1)
    VarBsp = Bsp^2 * (CV2P0 + CV2A + CV2R + CV2S + CV2F1 + CV2W1)

    #Return outputs as a list
    return(list(VarNfem=VarNfem, VarBsp=VarBsp))

  }






  if(any(names(adult.pars) %in% "Region" ) & any(names(adult.pars) %in% "Time")){
    # List of results to be filled: Number of females, Biomass and number of females at weight class
    resultlist <- list(
      Nfem =  tidyr::expand(adult.pars, tidyr::nesting(Time = Time, Region = Region),  Nfem = NA, SD = NA ),
      Biomass = tidyr::expand(adult.pars, tidyr::nesting(Time = Time, Region = Region),  Biomass = NA, SD = NA )
    )



    # Loop over Time and Region and save estimates and variances to correct position in resultlist
    for(i in unique(adult.pars$Time)){
      for(j in unique(adult.pars$Region)){

        tmp <-  dplyr::filter(adult.pars, Time == i, Region == j)

        if(nrow(tmp) == 0)next

        tmp_var <-  dplyr::filter(adult.vars, Time == i, Region == j)


        res <- BspNfem(P0 = tmp$P0, A = tmp$A, R = tmp$R, S = tmp$S,
                       Ww = tmp$W, Fw = tmp$.F)



        vars <-  VarBspNfemnw1(P0 = tmp$P0, A = tmp$A, R = tmp$R, S = tmp$S, F1 = tmp$.F, W1= tmp$W,
                               VP0 = tmp_var$P0, VA= tmp_var$A, VR= tmp_var$R, VS= tmp_var$S,
                               VF1= tmp_var$.F, VW1 = tmp_var$W)

        resultlist[["Nfem"]][resultlist[["Nfem"]]$Time == i & resultlist[["Nfem"]]$Region == j,
                             "Nfem"] <- res$Nfem

        resultlist[["Biomass"]][resultlist[["Biomass"]]$Time == i & resultlist[["Biomass"]]$Region == j,
                                "Biomass"] <- res$Bsp/1000



        resultlist[["Nfem"]][resultlist[["Nfem"]]$Time == i & resultlist[["Nfem"]]$Region == j,
                             "SD"] <- sqrt(vars$VarNfem)

        resultlist[["Biomass"]][resultlist[["Biomass"]]$Time == i & resultlist[["Biomass"]]$Region == j,
                                "SD"] <- sqrt(vars$VarBsp)/1000


      }
    }

  } else if(any(names(adult.pars) %in% "Time")){
    resultlist <- list(
      Nfem = expand.grid(Time = unique(adult.pars$Time), Nfem = NA, SD = NA ),
      Biomass = expand.grid(Time = unique(adult.pars$Time),  Biomass = NA , SD = NA)
    )

    for(i in unique(adult.pars$Time)){


      tmp <-  dplyr::filter(adult.pars,Time == i)

      if(nrow(tmp) == 0)next

      tmp_var <- dplyr::filter(adult.vars,Time == i)


      res <- BspNfem(P0 = tmp$P0, A = tmp$A, R = tmp$R, S = tmp$S,
                     Ww = tmp$W, Fw = tmp$.F)



      vars <-  VarBspNfemnw1(P0 = tmp$P0, A = tmp$A, R = tmp$R, S = tmp$S, F1 = tmp$.F, W1= tmp$W,
                             VP0 = tmp_var$P0, VA= tmp_var$A, VR= tmp_var$R, VS= tmp_var$S,
                             VF1= tmp_var$.F, VW1 = tmp_var$W)


      resultlist[["Nfem"]][resultlist[["Nfem"]]$Time == i, "Nfem"] <- res$Nfem

      resultlist[["Biomass"]][resultlist[["Biomass"]]$Time == i ,"Biomass"] <- res$Bsp/1000



      resultlist[["Nfem"]][resultlist[["Nfem"]]$Time == i ,"SD"] <- sqrt(vars$VarNfem)

      resultlist[["Biomass"]][resultlist[["Biomass"]]$Time == i,"SD"] <- sqrt(vars$VarBsp)/1000



    }


  } else if(any(names(adult.pars) %in% "Region")){
    resultlist <- list(
      Nfem = expand.grid(Time = unique(adult.pars$Region), Nfem = NA, SD = NA ),
      Biomass = expand.grid(Time = unique(adult.pars$Region),  Biomass = NA , SD = NA)
    )


    for(j in unique(adult.pars$Region)){

      tmp<- dplyr::filter(adult.pars, Region == j)

      if(nrow(tmp) == 0)next

      tmp_var <- dplyr::filter(adult.vars, Region == j)


      res <- BspNfem(P0 = tmp$P0, A = tmp$A, R = tmp$R, S = tmp$S,
                     Ww = tmp$W, Fw = tmp$.F)



      vars <-  VarBspNfemnw1(P0 = tmp$P0, A = tmp$A, R = tmp$R, S = tmp$S, F1 = tmp$.F, W1= tmp$W,
                             VP0 = tmp_var$P0, VA= tmp_var$A, VR= tmp_var$R, VS= tmp_var$S,
                             VF1= tmp_var$.F, VW1 = tmp_var$W)


      resultlist[["Nfem"]][resultlist[["Nfem"]]$Region == j,"Nfem"] <- res$Nfem

      resultlist[["Biomass"]][resultlist[["Biomass"]]$Region == j,"Biomass"] <- res$Bsp/1000



      resultlist[["Nfem"]][resultlist[["Nfem"]]$Region == j,"SD"] <- sqrt(vars$VarNfem)

      resultlist[["Biomass"]][ resultlist[["Biomass"]]$Region == j,"SD"] <- sqrt(vars$VarBsp)/1000



    }


  } else {

    res <- BspNfem(P0 = adult.pars$P0, A = adult.pars$A, R = adult.pars$R, S = adult.pars$S,
                   Ww = adult.pars$W, Fw = adult.pars$.F)

    vars <-  VarBspNfemnw1(P0 = adult.pars$P0, A = adult.pars$A, R = adult.pars$R, S = adult.pars$S, F1 = adult.pars$.F, W1= adult.pars$W,
                           VP0 = adult.vars$P0, VA= adult.vars$A, VR= adult.vars$R, VS= adult.vars$S, VF1= adult.vars$.F, VW1 = adult.vars$W)

    resultlist <- list(Nfem = data.frame(Nfem = res$Nfem,
                                         SD = sqrt(vars$VarNfem)),
                       Biomass = data.frame(Biomass = res$Bsp/1000,
                                            SD = sqrt(vars$VarBsp)/1000))

  }

  return(resultlist)
}




#' Estimate mean weight (W) and mean (F) for the standard DEPM approach.
#' @description The standard DEPM approach does not use weight bins and does not require a set of weight class parameters.
#'     However, because of this it requires additional parameters for mean female weight (W) and mean batch fecundity (F).
#'     These are produced using this function and can be combined with the other weight invariant parameters (P0, A, S and R)
#'     using the combine_estimates() function. Two data.frames are required which contain weight data and batch fecundity data
#'     which can be used to estimated each of these parameters.
#' @param weight.data Adult data used to calculate the mean total female weight (W). This data set should ONLY include
#'     females and should include variables for Total weight and/or gonad free weight (Total weight - Gonad weight). These
#'     are used to calculate mean W and provide estimates the intended weight estimates to the fecundity estimator.
#' @param TotalWt The column heading in weight.data that represents Total Weight in grams
#' @param GonadFrWt The column heading in weight.data that represents Gonad Free Weight in grams. If Gonad weight is not being used
#'     in the fecundity estimator (as total weight is used instead), then this parameter does not need to be specified.
#' @param Time A string containing the column name for a timestep if desired as a grouping variable
#' @param Region A string containing the column name for a region if desired as a grouping variable
#'
#' @param parameters  A list of 4 starting parameters for the batch fecundity relationship that
#'     must include: "alpha", "beta", "Sigma0" and "Sigma1"
#'
#' @return A data.frame that includes the Mean W, variance of W, Mean F, variance of F, the covariance of the fecundity-at-weight data and Time and Region if grouping variables were provided
#' @export

Estimate_mean_W_F <- function(weight.data, TotalWt, GonadFrWt= NULL, Time = NULL, Region = NULL, parameters=NULL){
  if(!is.null(Region) & !any(names(weight.data) %in% Region)) stop("Region column could not be determined")
  if(!is.null(Time) & !any(names(weight.data) %in% Time)) stop("Time column could not be determined")
  if(all(!names(weight.data) %in% TotalWt)) stop("Total Weight column incorrectly specified")
  if(!is.null(GonadFrWt) & all(!names(weight.data) %in% GonadFrWt)) stop("Gonad Free Weight column incorrectly specified")
  if(is.null(GonadFrWt)) {
    message("Total weight used to estimate fecundity as gonad free weight was not provided")
  }else{
    if(GonadFrWt == TotalWt) {
      GonadFrWt <- NULL # Total weight will be used automatically in this situation so make explicit to avoid bugs
      message("Total weight used to estimate fecundity as gonad free weight was not provided")
    }
  }
  if(is.null(parameters)) stop("Batch fecundity relationship parameters must be provided from the Estimate_batch_fecundity() function")
  if(!is.data.frame(parameters)) stop("Batch fecundity relationship parameters should be list obtained from the Estimate_batch_fecundity() function with 'return.parameters == TRUE'")


  # extract individual parameters from the parameter list
  alpha <- parameters[parameters$Parameter== "alpha","Val"]
  beta <- parameters[parameters$Parameter== "beta","Val"]
  Sigma0 <- parameters[parameters$Parameter== "Sigma0","Val"]
  Sigma1 <- parameters[parameters$Parameter== "Sigma1","Val"]


  if(is.null(Time) & !is.null(Region)){
    processed_data <- dplyr::select(weight.data,
                                    Region = paste(Region),
                                    GonadFrWt =  paste(GonadFrWt),
                                    TotalWt = paste(TotalWt))

    processed_data <- na.omit(processed_data)

    results <- expand.grid(Region = unique(processed_data$Region),
                           Mean_W = NA,
                           var_W = NA,
                           Mean_F= NA,
                           var_F= NA)

    for(i in unique(processed_data$Region)){
      tmp <- dplyr::filter(processed_data, Region == i)
      TotalWt <- tmp$TotalWt
      results[which(results$Region == i),"Mean_W"]  <- mean(TotalWt)
      results[which(results$Region == i),"var_W"]  <- var(TotalWt)/length(TotalWt) #variance of the mean

      if(is.null(GonadFrWt)){
        prediction_wt <- TotalWt

      }else{
        prediction_wt <- tmp$GonadFrWt
      }

      # remove data that may not be in grams
      prediction_wt <- prediction_wt[prediction_wt >1]

      # The fecundity relationship is not time invariant but the mean weight used is.
      Fec_results <- data.frame(Wt = prediction_wt,
                                Fecundity = alpha * prediction_wt^beta)
      Fec_results <-  dplyr::mutate(Fec_results, Var = (Sigma0*Fecundity^Sigma1)^2)

      results[which(results$Region == i),"Mean_F"] <- mean(Fec_results$Fecundity)
      results[which(results$Region == i),"var_F"] <- (1/length(Fec_results$Fecundity)^2)*(sum(Fec_results$var))
      #var(Fec_results$Fecundity)/length(Fec_results$Fecundity)

    }


  } else if(!is.null(Time) & is.null(Region)) {

    processed_data <- dplyr::select(weight.data,
                                    Time = paste(Time),
                                    GonadFrWt =  paste(GonadFrWt),
                                    TotalWt = paste(TotalWt))
    processed_data <- na.omit(processed_data)

    results <- expand.grid(Time = unique(processed_data$Time),
                           Mean_W = NA,
                           var_W = NA,
                           Mean_F= NA,
                           var_F= NA)

    for(i in unique(processed_data$Time)){
      tmp <- dplyr::filter(processed_data, Time == i)
      TotalWt <- tmp$TotalWt
      results[which(results$Time == i),"Mean_W"]  <- mean(TotalWt)
      results[which(results$Time == i),"var_W"]  <- var(TotalWt)/length(TotalWt)

      if(is.null(GonadFrWt)){
        prediction_wt <- TotalWt
      }else{
        prediction_wt <- tmp$GonadFrWt
      }

      # remove data that may not be in grams
      prediction_wt <- prediction_wt[prediction_wt >1]

      Fec_results <- data.frame(Wt = prediction_wt,
                            Fecundity = alpha * prediction_wt^beta)
      Fec_results <-  dplyr::mutate(Fec_results, Var = (Sigma0*Fecundity^Sigma1)^2)

      results[which(results$Time == i),"Mean_F"] <- mean(Fec_results$Fecundity)
      results[which(results$Time == i),"var_F"] <- (1/length(Fec_results$Fecundity)^2)*(sum(Fec_results$var))
      #var(Fec_results$Fecundity)/length(Fec_results$Fecundity)
    }

  } else if(!is.null(Time) & !is.null(Region)) {

    processed_data <- dplyr::select(weight.data,
                                    Time = paste(Time),
                                    Region = paste(Region),
                                    GonadFrWt =  paste(GonadFrWt),
                                    TotalWt = paste(TotalWt))
    processed_data <- na.omit(processed_data)

    results <-  tidyr::expand(processed_data, tidyr::nesting(Time = Time, Region = Region),
                              Mean_W = NA,
                              var_W = NA,
                              Mean_F= NA,
                              var_F= NA)

    for(i in unique(processed_data$Time)){
      for(j in unique(processed_data$Region)){
        tmp <- dplyr::filter(processed_data, Time == i, Region == j)
        TotalWt <- tmp$TotalWt
        results[which(results$Time == i & results$Region == j),"Mean_W"]  <- mean(TotalWt)
        results[which(results$Time == i & results$Region == j),"var_W"]  <- var(TotalWt)/length(TotalWt)

        if(is.null(GonadFrWt)){
          prediction_wt <- TotalWt
        }else{
          prediction_wt <- tmp$GonadFrWt
        }

        # remove data that may not be in grams
        prediction_wt <- prediction_wt[prediction_wt >1]

        Fec_results <- data.frame(Wt = prediction_wt,
                                  Fecundity = alpha * prediction_wt^beta)
        Fec_results <-  dplyr::mutate(Fec_results, Var = (Sigma0*Fecundity^Sigma1)^2)

        results[which(results$Time == i & results$Region == j),"Mean_F"] <- mean(Fec_results$Fecundity)
        results[which(results$Time == i & results$Region == j),"var_F"] <- (1/length(Fec_results$Fecundity)^2)*(sum(Fec_results$var))
        #var(Fec_results$Fecundity)/length(Fec_results$Fecundity)
      }
    }

  } else {
    processed_data <- dplyr::select(weight.data,
                                    GonadFrWt =  paste(GonadFrWt),
                                    TotalWt = paste(TotalWt))

    processed_data <- na.omit(processed_data)

    TotalWt <- processed_data$TotalWt
    Mean_W  <- mean(TotalWt)
    var_W  <- var(TotalWt)/length(TotalWt)

    if(is.null(GonadFrWt)){
      prediction_wt <- TotalWt
    }else{
      prediction_wt <- tmp$GonadFrWt
    }

    # remove data that may not be in grams
    prediction_wt <- prediction_wt[prediction_wt >1]

    Fec_results <- data.frame(Wt = prediction_wt,
                              Fecundity = alpha * prediction_wt^beta)
    Fec_results <-  dplyr::mutate(Fec_results, Var = (Sigma0*Fecundity^Sigma1)^2)

    Mean_F <- mean(Fec_results$Fecundity)
    var_F <- (1/length(Fec_results$Fecundity)^2)*(sum(Fec_results$var))
    #var(Fec_results$Fecundity)/length(Fec_results$Fecundity)

    results <- data.frame(Mean_W, var_W, Mean_F,var_F )
  }

  return(results)

}



#' Estimate mean weight (W) and mean (F) for the standard DEPM approach.
#' @description The standard DEPM approach does not use weight bins and does not require a set of weight class parameters.
#'     However, because of this it requires additional parameters for mean female weight (W) and mean batch fecundity (F).
#'     These are produced using this function and can be combined with the other weight invariant parameters (P0, A, S and R)
#'     using the combine_estimates() function. Two data.frames are required which contain weight data and batch fecundity data
#'     which can be used to estimated each of these parameters.
#' @param weight.data Adult data used to calculate the mean total female weight (W). This data set should ONLY include
#'     females and should include variables for Total weight and/or gonad free weight (Total weight - Gonad weight). These
#'     are used to calculate mean W and provide estimates the intended weight estimates to the fecundity estimator.
#' @param TotalWt The column heading in weight.data that represents Total Weight in grams
#' @param GonadFrWt The column heading in weight.data that represents Gonad Free Weight in grams. If Gonad weight is not being used
#'     in the fecundity estimator (as total weight is used instead), then this parameter does not need to be specified.
#' @param Time A string containing the column name for a timestep if desired as a grouping variable
#' @param Region A string containing the column name for a region if desired as a grouping variable
#' @param fecundity.data A data frame of weights and batch fecundity used to determine the relationship
#'     which informs mean F. The weight data can be in Gonad free weight or total weight. If Gonad free weight is used,
#'     then `GonadFrWt` needs to be specified as otherwise total weight will be incorrectly used to estimate fecundity.
#'
#' @param parameters  A list of 4 starting parameters for the batch fecundity relationship that
#'     must include: "alpha", "beta", "Sigma0" and "Sigma1"
#'
#' @return A data.frame that includes the Mean W, variance of W, Mean F, variance of F, the covariance of the fecundity-at-weight data and Time and Region if grouping variables were provided
#' @export

#
# Estimate_mean_W_F <- function(weight.data, TotalWt, GonadFrWt= NULL, Time = NULL, Region = NULL, fecundity.data, parameters){
#   if(!is.null(Region) & !any(names(weight.data) %in% Region)) stop("Region column could not be determined")
#   if(!is.null(Time) & !any(names(weight.data) %in% Time)) stop("Time column could not be determined")
#   if(all(!names(weight.data) %in% TotalWt)) stop("Total Weight column incorrectly specified")
#   if(!is.null(GonadFrWt) & all(!names(weight.data) %in% GonadFrWt)) stop("Gonad Free Weight column incorrectly specified")
#   if(is.null(GonadFrWt)) {
#     message("Total weight used to estimate fecundity as gonad free weight was not provided")
#   }else{
#     if(GonadFrWt == TotalWt) {
#       GonadFrWt <- NULL # Total weight will be used automatically in this situation so make explicit to avoid bugs
#       message("Total weight used to estimate fecundity as gonad free weight was not provided")
#     }
#   }
#
#
#
#   if(is.null(Time) & !is.null(Region)){
#     processed_data <- dplyr::select(weight.data,
#                                     Region = paste(Region),
#                                     GonadFrWt =  paste(GonadFrWt),
#                                     TotalWt = paste(TotalWt))
#
#     processed_data <- na.omit(processed_data)
#
#     results <- expand.grid(Region = unique(processed_data$Region),
#                            Mean_W = NA,
#                            var_W = NA,
#                            Mean_F= NA,
#                            var_F= NA)
#
#     for(i in unique(processed_data$Region)){
#       tmp <- dplyr::filter(processed_data, Region == i)
#       TotalWt <- tmp$TotalWt
#       results[which(results$Region == i),"Mean_W"]  <- mean(TotalWt)
#       results[which(results$Region == i),"var_W"]  <- var(TotalWt)
#
#       if(is.null(GonadFrWt)){
#         prediction_wt <- TotalWt
#
#       }else{
#         prediction_wt <- tmp$GonadFrWt
#       }
#
#       # remove data that may not be in grams
#       prediction_wt <- prediction_wt[prediction_wt >1]
#
#       # The fecundity relationship is not time invariant but the mean weight used is.
#       Fec_results <- Estimate_Batch_Fecundity(fecundity.data, parameters,
#                                               prediction.int = prediction_wt,
#                                               verbose = FALSE)
#       results[which(results$Region == i),"Mean_F"] <- mean(Fec_results$Fecundity)
#       results[which(results$Region == i),"var_F"] <- var(Fec_results$Fecundity)
#
#     }
#
#
#   } else if(!is.null(Time) & is.null(Region)) {
#
#     processed_data <- dplyr::select(weight.data,
#                                     Time = paste(Time),
#                                     GonadFrWt =  paste(GonadFrWt),
#                                     TotalWt = paste(TotalWt))
#     processed_data <- na.omit(processed_data)
#
#     results <- expand.grid(Time = unique(processed_data$Time),
#                            Mean_W = NA,
#                            var_W = NA,
#                            Mean_F= NA,
#                            var_F= NA)
#
#     for(i in unique(processed_data$Time)){
#       tmp <- dplyr::filter(processed_data, Time == i)
#       TotalWt <- tmp$TotalWt
#       results[which(results$Time == i),"Mean_W"]  <- mean(TotalWt)
#       results[which(results$Time == i),"var_W"]  <- var(TotalWt)
#
#       if(is.null(GonadFrWt)){
#         prediction_wt <- TotalWt
#       }else{
#         prediction_wt <- tmp$GonadFrWt
#       }
#
#       # remove data that may not be in grams
#       prediction_wt <- prediction_wt[prediction_wt >1]
#
#       # The fecundity relationship is not time invariant but the mean weight used is.
#       Fec_results <- Estimate_Batch_Fecundity(fecundity.data, parameters,
#                                               prediction.int = prediction_wt,
#                                               verbose = FALSE)
#       results[which(results$Time == i),"Mean_F"] <- mean(Fec_results$Fecundity)
#       results[which(results$Time == i),"var_F"] <- var(Fec_results$Fecundity)
#     }
#
#   } else if(!is.null(Time) & !is.null(Region)) {
#
#     processed_data <- dplyr::select(weight.data,
#                                     Time = paste(Time),
#                                     Region = paste(Region),
#                                     GonadFrWt =  paste(GonadFrWt),
#                                     TotalWt = paste(TotalWt))
#     processed_data <- na.omit(processed_data)
#
#     results <-  tidyr::expand(processed_data, tidyr::nesting(Time = Time, Region = Region),
#                               Mean_W = NA,
#                               var_W = NA,
#                               Mean_F= NA,
#                               var_F= NA)
#
#     for(i in unique(processed_data$Time)){
#       for(j in unique(processed_data$Region)){
#         tmp <- dplyr::filter(processed_data, Time == i, Region == j)
#         TotalWt <- tmp$TotalWt
#         results[which(results$Time == i & results$Region == j),"Mean_W"]  <- mean(TotalWt)
#         results[which(results$Time == i & results$Region == j),"var_W"]  <- var(TotalWt)
#
#         if(is.null(GonadFrWt)){
#           prediction_wt <- TotalWt
#         }else{
#           prediction_wt <- tmp$GonadFrWt
#         }
#
#         # remove data that may not be in grams
#         prediction_wt <- prediction_wt[prediction_wt >1]
#
#         # The fecundity relationship is not time invariant but the mean weight used is.
#         Fec_results <- Estimate_Batch_Fecundity(fecundity.data, parameters,
#                                                 prediction.int = prediction_wt,
#                                                 verbose = FALSE)
#         results[which(results$Time == i & results$Region == j),"Mean_F"] <- mean(Fec_results$Fecundity)
#         results[which(results$Time == i & results$Region == j),"var_F"] <- var(Fec_results$Fecundity)
#       }
#     }
#
#   } else {
#     processed_data <- dplyr::select(weight.data,
#                                     GonadFrWt =  paste(GonadFrWt),
#                                     TotalWt = paste(TotalWt))
#
#     processed_data <- na.omit(processed_data)
#
#     TotalWt <- processed_data$TotalWt
#     Mean_W  <- mean(TotalWt)
#     var_W  <- var(TotalWt)
#
#     if(is.null(GonadFrWt)){
#       prediction_wt <- TotalWt
#     }else{
#       prediction_wt <- tmp$GonadFrWt
#     }
#
#     # remove data that may not be in grams
#     prediction_wt <- prediction_wt[prediction_wt >1]
#
#     # The fecundity relationship is not time invariant but the mean weight used is.
#     Fec_results <- Estimate_Batch_Fecundity(fecundity.data, parameters,
#                                             prediction.int = prediction_wt,
#                                             verbose = FALSE)
#     Mean_F <- mean(Fec_results$Fecundity)
#     var_F <- var(Fec_results$Fecundity)
#
#     results <- data.frame(Mean_W, var_W, Mean_F,var_F )
#   }
#
#   return(results)
#
# }

