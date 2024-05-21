
###################################################################################################
# Function name: ps_calculation
# Purpose: Calculating propensity scores, IPT weights, and generating visualization
# Input: 1. Treatment indicator (string): Yes has been coded as 1 and No as 0 (factorized, labeled, 
#           and set ideal reference level)
#        2. Function formula for the propensity score estimation (string): treatment ~ confounders 
#           (default: using logistic regression model)
#           - Confounders should all be cleaned as: numerical or categorical (factorized, labeled, 
#             and set ideal reference level)
#           - Need to at least include all confounders in the propensity score model
#        3. Data
#        4. If estimation of propensity score 
#        5. If estimation of weights (outside propensity score is needed) 
#        6. If plot distribution of the estimated propensity score for two groups
#        7. If plot distribution of the estimated IPT weights for two groups         
# Output: a list, consist of (estimated propensity scores, estimated IPT weights, distribution plot 
#         for propensity scores, distribution plot for weights)
#        1. If estimation_ps = TRUE (treatment, formula, data are required), return the estimated 
#           propensity score based on current treatment status
#        2. If estimation_weight = TRUE (treatment, data, outside propensity score are required), 
#           return the estimated IPT weights based on current treatment status
#        3. If plot_ps = TRUE, print the distribution of propensity score
#        4. If plot_weight = TRUE, print the distribution of IPT weights
# Required package: ggplot2
###################################################################################################

###################################################################################################
library(ggplot2)

ps_calculation <- function(treatment, data, formula = NA, ps_outside = NA, weight_outside = NA,
                           estimation_ps = FALSE, estimation_weight = FALSE, 
                           plot_ps = FALSE, plot_weight = FALSE) {
  # when propensity score estimation 
  if (estimation_ps == TRUE) {
    ps_model <- glm(formula, data = data, family = binomial)
    ps_score <- predict(ps_model, type = "response")
  } else {
    ps_score <- NA
  }
  
  # when weights calculation
  if (estimation_weight == TRUE & estimation_ps == TRUE) {
    weight <- ifelse(data[,treatment] == levels(data[,treatment])[2], 1/(ps_score), 1/(1-ps_score))
  } else if (estimation_weight == TRUE & estimation_ps == FALSE) {
    weight <- ifelse(data[,treatment] == levels(data[,treatment])[2], 1/(ps_outside), 1/(1-ps_outside))
  } else if (estimation_weight == FALSE) {
    weight <- NA
  }
  
  # when plot the distribution of propensity score
  if (plot_ps == TRUE & estimation_ps == TRUE) {
    ps_dist <- ggplot(data.frame(Ps_score = ps_score, Group = data[,treatment]), aes(x = Ps_score)) + 
      geom_density(aes(color = Group, fill = Group), alpha = 0.4) +
      scale_color_manual(values = c("#00AFBB","#E7B800")) +
      scale_fill_manual(values = c("#00AFBB","#E7B800")) + xlim(0,1) + 
      xlab("Estimated Propensity Score") + ylab("Density") 
  } else if (plot_ps == TRUE & estimation_ps == FALSE) {
    ps_dist <- ggplot(data.frame(Ps_score = ps_outside, Group = data[,treatment]), aes(x = Ps_score)) + 
      geom_density(aes(color = Group, fill = Group), alpha = 0.4) +
      scale_color_manual(values = c("#00AFBB","#E7B800")) +
      scale_fill_manual(values = c("#00AFBB","#E7B800")) + xlim(0,1) + 
      xlab("Estimated Propensity Score") + ylab("Density") 
  } else if (plot_ps == FALSE) {
    ps_dist <- NA
  }
  
  # when plot the distribution of weights
  if (plot_weight == TRUE & estimation_weight == TRUE) {
    weight_dist <- ggplot(data.frame(IPT_weight = weight, Group = data[,treatment]), aes(x = IPT_weight)) + 
      geom_density(aes(color = Group, fill = Group), alpha = 0.4) +
      scale_color_manual(values = c("#00AFBB","#E7B800")) +
      scale_fill_manual(values = c("#00AFBB","#E7B800")) + 
      xlab("Estimated IPT weights") + ylab("Density") 
  } else if (plot_weight == TRUE & estimation_weight == FALSE) {
    weight_dist <- ggplot(data.frame(IPT_weight = weight_outside, Group = data[,treatment]), aes(x = IPT_weight)) + 
      geom_density(aes(color = Group, fill = Group), alpha = 0.4) +
      scale_color_manual(values = c("#00AFBB","#E7B800")) +
      scale_fill_manual(values = c("#00AFBB","#E7B800")) + 
      xlab("Estimated IPT weight") + ylab("Density") 
  } else if (plot_weight == FALSE) {
    weight_dist <- NA
  }
  return(list(ps_score,weight,ps_dist,weight_dist))
}
###################################################################################################


###################################################################################################
# Function name: smd_calculation
# Purpose: Calculating the SMD for the unweighted and weighted population
# Input: 1. Treatment indicator (string): Yes has been coded as 1 and No as 0 (factorized, labeled, 
#           and set ideal reference level)
#        2. Weight: IPT weights, overlap weights, or other weights to create a weighted population
#        3. Names of variables that should be examined balance (vector)
#        4. Data
#        5. If weighted population or the original population 
#        6. If use weighted variance to standardize when weighted population
#        7. type: continuous, binary, categorical
# Output: a dataframe contains SMD
#        1. If weighted_pop = TRUE (weight are required), return a data frame with a list of 
#           variables and its SMD value for the weighted population
#        2. If weighted_pop = FALSE (weight can be NA), return a data frame with a list of 
#           variables and its SMD value for the original population
###################################################################################################

###################################################################################################
smd_calculation <- function(treatment, weight = NA, var_list, data, weighted_pop = FALSE, weighted_variance = FALSE, type = "continuous") {
  
  # define the functions for calculating weighted variance
  weighted.var.continuous <- function(x, w) {
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w) / sum(w)
    (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2)
  }
  
  weighted.var.binary <- function(x, w) {
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w) / sum(w)
    (sum.w^2 / (sum.w^2 - sum.w2)) * mean.w * (1-mean.w)
  }
  
  ind_trt <- as.numeric(data[,treatment]) - 1 == 1
  ind_ctrl <- as.numeric(data[,treatment]) - 1 == 0
  
  weighted_n_trt <- sum(weight[ind_trt])
  weighted_n_ctrl <- sum(weight[ind_ctrl])
  
  smd <- NA
  
  if (type == "continuous") {
    if (weighted_pop == FALSE) {
      for (i in 1:length(var_list)) {
        smd[i] <- abs(mean(data[ind_trt, var_list[i]]) - mean(data[ind_ctrl, var_list[i]]))/
          sqrt((var(data[ind_trt, var_list[i]]) + var(data[ind_ctrl, var_list[i]]))/2)
      }
    }
    if (weighted_pop == TRUE & weighted_variance == FALSE) {
      for (i in 1:length(var_list)) {
        smd[i] <- abs(sum(data[ind_trt, var_list[i]] * weight[ind_trt])/weighted_n_trt  - 
                        sum(data[ind_ctrl, var_list[i]] * weight[ind_ctrl])/weighted_n_ctrl)/
          sqrt((var(data[ind_trt, var_list[i]]) + var(data[ind_ctrl, var_list[i]]))/2)
      }
    }
    if (weighted_pop == TRUE & weighted_variance == TRUE) {
      for (i in 1:length(var_list)) {
        smd[i] <- abs(sum(data[ind_trt, var_list[i]] * weight[ind_trt])/weighted_n_trt  - 
                        sum(data[ind_ctrl, var_list[i]] * weight[ind_ctrl])/weighted_n_ctrl)/
          sqrt((weighted.var.continuous(data[ind_trt, var_list[i]], weight[ind_trt])+
                  weighted.var.continuous(data[ind_ctrl, var_list[i]],weight[ind_ctrl]))/2)
      }
    }
    return(data.frame(Variable = var_list, SMD = smd))
  }
  
  ########
  
  categories <- NA
  k <- 0
  
  if (type == "binary") {
    if (weighted_pop == FALSE) {
      for (i in 1:length(var_list)) {
        j <- levels(data[,var_list[i]])[2]
        k <- k+1
        temp <- rep(0, nrow(data)) # indicator for not in this class of the given variable
        temp[data[,var_list[i]] == j] <- 1
        smd[k] <- abs(mean(temp[ind_trt]) - mean(temp[ind_ctrl]))/
          sqrt((mean(temp[ind_trt])*(1-mean(temp[ind_trt])) + 
                  mean(temp[ind_ctrl])*(1-mean(temp[ind_ctrl])))/2)
        categories[k] <- paste(var_list[i],"-",j)
      }
    }
    if (weighted_pop == TRUE & weighted_variance == FALSE) {
      for (i in 1:length(var_list)) {
        j <- levels(data[,var_list[i]])[2]
        k <- k+1
        temp <- rep(0, nrow(data)) # indicator for not in this class of the given variable
        temp[data[,var_list[i]] == j] <- 1
        smd[k] <- abs(sum(temp[ind_trt] * weight[ind_trt])/weighted_n_trt - 
                        sum(temp[ind_ctrl] * weight[ind_ctrl])/weighted_n_ctrl)/
          sqrt((mean(temp[ind_trt])*(1-mean(temp[ind_trt])) + 
                  mean(temp[ind_ctrl])*(1-mean(temp[ind_ctrl])))/2)
        categories[k] <- paste(var_list[i],"-",j)
      }
    }
    if (weighted_pop == TRUE & weighted_variance == TRUE) {
      for (i in 1:length(var_list)) {
        j <- levels(data[,var_list[i]])[2]
        k <- k+1
        temp <- rep(0, nrow(data)) # indicator for not in this class of the given variable
        temp[data[,var_list[i]] == j] <- 1
        smd[k] <- abs(sum(temp[ind_trt] * weight[ind_trt])/weighted_n_trt  - 
                        sum(temp[ind_ctrl] * weight[ind_ctrl])/weighted_n_ctrl)/
          sqrt((weighted.var.binary(temp[ind_trt], weight[ind_trt])+
                  weighted.var.binary(temp[ind_ctrl],weight[ind_ctrl]))/2)
        categories[k] <- paste(var_list[i],"-",j)
      }    
    }
    return(data.frame(Variable = categories, SMD = smd))
  }
  
  ########
  
  if (type == "categorical") {
    if (weighted_pop == FALSE) {
      for (i in 1:length(var_list)) {
        for (j in levels(data[,var_list[i]])) {
          k <- k+1
          temp <- rep(0, nrow(data)) # indicator for not in this class of the given variable
          temp[data[,var_list[i]] == j] <- 1
          smd[k] <- abs(mean(temp[ind_trt]) - mean(temp[ind_ctrl]))/
            sqrt((mean(temp[ind_trt])*(1-mean(temp[ind_trt])) + 
                    mean(temp[ind_ctrl])*(1-mean(temp[ind_ctrl])))/2)
          categories[k] <- paste(var_list[i],"-",j)
        }
      }
    }
    if (weighted_pop == TRUE & weighted_variance == FALSE) {
      for (i in 1:length(var_list)) {
        for (j in levels(data[,var_list[i]])) {
          k <- k+1
          temp <- rep(0, nrow(data)) # indicator for not in this class of the given variable
          temp[data[,var_list[i]] == j] <- 1
          smd[k] <- abs(sum(temp[ind_trt] * weight[ind_trt])/weighted_n_trt  - 
                          sum(temp[ind_ctrl] * weight[ind_ctrl])/weighted_n_ctrl)/
            sqrt((mean(temp[ind_trt])*(1-mean(temp[ind_trt])) + 
                    mean(temp[ind_ctrl])*(1-mean(temp[ind_ctrl])))/2)
          categories[k] <- paste(var_list[i],"-",j)
        }
      }
    }
    if (weighted_pop == TRUE & weighted_variance == TRUE) {
      for (i in 1:length(var_list)) {
        for (j in levels(data[,var_list[i]])) {
          k <- k+1
          temp <- rep(0, nrow(data)) # indicator for not in this class of the given variable
          temp[data[,var_list[i]] == j] <- 1
          smd[k] <- abs(sum(temp[ind_trt] * weight[ind_trt])/weighted_n_trt  - 
                          sum(temp[ind_ctrl] * weight[ind_ctrl])/weighted_n_ctrl)/
            sqrt((weighted.var.binary(temp[ind_trt], weight[ind_trt])+
                    weighted.var.binary(temp[ind_ctrl],weight[ind_ctrl]))/2)
          categories[k] <- paste(var_list[i],"-",j)
        }
      }    
    }
    return(data.frame(Variable = categories, SMD = smd))
  }
}
###################################################################################################