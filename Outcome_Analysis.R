
###################################################################################################
# Input data: data_final
# Key variables:
#             data_final$survtime: survival time in days, which is the time difference between 
#                                  time of first surgery and min(time of death, time of censoring)
#             data_final$status: censoring indicator, 1 means time of death <= time of censoring
#             data_final$ad_chemo: if there is adjuvant chemotherapy 
#             data_final$ad_chemo: if there is adjuvant radiotherapy
#             data_final$main_surgery: treatment indicator
# Covariates:
#             data_final$agedx, data_final$sex, data_final$race, data_final$marriage, 
#             data_final$morphology, data_final$stage, data_final$primary_site, data_final$size
#             data_final$MEDhistory
#
# Required package: survival, survminer, ggplot2
# Required function from the other file: ps_calculation, smd_calculation
###################################################################################################


###########################
### read processed data ###
###########################
data_final <- read.csv(".../sample_final.csv", row.names = 1)


###############################
### subgroup classification ###
###############################
# data_final <- data_final[data_final$ad_chemo == 0,]
data_final <- data_final[data_final$ad_radiation == 0,]


#########################
### survival analysis ###
#########################
library(survival)
library(survminer)

fit <- survfit(Surv(survtime/365,status) ~ main_surgery, data_final)
ggsurvplot(fit, data = data_final, legend.title = "Groups", 
           legend.labs = c("Lobectomy","Limited Resection"), 
           palette = c("#00AFBB","#E7B800"), conf.int = TRUE, 
           ggtheme = theme_bw(), censor = FALSE) + xlab("Years since surgery")

# log rank test
survdiff(Surv(survtime,status) ~ main_surgery, data_final)


#########################################################
### propensity score analysis (use stabilized weight) ###
#########################################################
out <- ps_calculation(treatment = "main_surgery", 
                      formula = "main_surgery ~ agedx + sex + race + marriage + morphology + 
                                stage + primary_site + size + MEDhistory", 
                      data = data_final, 
                      estimation_ps = TRUE, estimation_weight = FALSE, 
                      plot_ps = TRUE, plot_weight = FALSE)

propensity_score <- out[[1]]

# plot the distribution of the estimated propensity score
out[[3]] 

# weight calculation
weight <- ifelse(data_final$main_surgery == levels(data_final$main_surgery)[2],
                 sum(data_final$main_surgery == levels(data_final$main_surgery)[2])/nrow(data_final)/(propensity_score), 
                 sum(data_final$main_surgery == levels(data_final$main_surgery)[1])/nrow(data_final)/(1-propensity_score))

# truncate weight at upper 99.5% quantile
weight[weight > quantile(weight,0.995)] <- quantile(weight,0.995)

out2 <- ps_calculation(treatment = "main_surgery", formula = NA, data = data_final, ps_outside = NA, 
                       weight_outside = weight, estimation_ps = FALSE, estimation_weight = FALSE, 
                       plot_ps = FALSE, plot_weight = TRUE)

# plot the distribution of the estimated weights
out2[[4]] # 600*400

# calculate smd values
smd1 <- smd_calculation(treatment = "main_surgery", weight = weight, var_list = c("sex","race","marriage","stage"), 
                data = data_final, weighted_pop = FALSE, weighted_variance = FALSE, type = "binary")

smd2 <- smd_calculation(treatment = "main_surgery", weight = weight, var_list = c("morphology","primary_site"), 
                data = data_final, weighted_pop = FALSE, weighted_variance = FALSE, type = "categorical")

smd3 <- smd_calculation(treatment = "main_surgery", weight = weight, var_list = c("agedx","size","MEDhistory"), 
                data = data_final, weighted_pop = FALSE, weighted_variance = FALSE, type = "continuous")

smd4 <- smd_calculation(treatment = "main_surgery", weight = weight, var_list = c("sex","race","marriage","stage"), 
                data = data_final, weighted_pop = TRUE, weighted_variance = TRUE, type = "binary")

smd5 <- smd_calculation(treatment = "main_surgery", weight = weight, var_list = c("morphology","primary_site"), 
                data = data_final, weighted_pop = TRUE, weighted_variance = TRUE, type = "categorical")

smd6 <- smd_calculation(treatment = "main_surgery", weight = weight, var_list = c("agedx","size","MEDhistory"), 
                data = data_final, weighted_pop = TRUE, weighted_variance = TRUE, type = "continuous")

data_smd <- data.frame(rbind(smd1,smd2,smd3,smd4,smd5,smd6), 
                       Weighting = c(rep("Before",nrow(rbind(smd1,smd2,smd3))), rep("After",nrow(rbind(smd1,smd2,smd3)))))

data_smd$Weighting <- factor(data_smd$Weighting, levels = c("Before","After"))

data_smd$Variable <- c("Sex - male","Race - other","Marital status - other",
                       "Stage - stage II","Morphology - adenocarcinoma","Morphology - large cell carcinoma",
                       "Morphology - squamous cell carcinoma","Morphology - other specified carcinoma",
                       "Morphology - unspecified NSCLC","Primary site - upper lobe","Primary site - middle lobe",
                       "Primary site - lower lobe","Primary site - others","Age at diagnosis","Size","Modified comorbidity score",
                       "Sex - male","Race - other","Marital status - other",
                       "Stage - stage II","Morphology - adenocarcinoma","Morphology - large cell carcinoma",
                       "Morphology - squamous cell carcinoma","Morphology - other specified carcinoma",
                       "Morphology - unspecified NSCLC","Primary site - upper lobe","Primary site - middle lobe",
                       "Primary site - lower lobe","Primary site - others","Age at diagnosis","Size","Modified comorbidity score")

# plot the smd values before and after weighting
ggplot(data_smd, aes(Variable, SMD)) + geom_point(aes(colour = Weighting)) + coord_flip() + 
  geom_hline(yintercept = 0.1, color = "navy", size = 1, linetype = "dashed") +
  theme(axis.title.y = element_blank())


#########################
### Survival analysis ###
#########################
fit_cox_unweighted1 <- coxph(Surv(survtime,status) ~ main_surgery, data = data_final)
fit_cox_unweighted2 <- coxph(Surv(survtime,status) ~ main_surgery + agedx + sex + race + marriage + morphology + 
                              stage + size + primary_site + MEDhistory, data = data_final)
fit_cox_weighted1 <- coxph(Surv(survtime,status) ~ main_surgery, data = data_final, weight = weight)
fit_cox_weighted2 <- coxph(Surv(survtime,status) ~ main_surgery + agedx + sex + race + marriage + morphology + 
                   stage + size + primary_site + MEDhistory, data = data_final, weight = weight)
round(summary(fit_cox_unweighted1)$conf.int,2)
round(summary(fit_cox_unweighted2)$conf.int,2)
round(summary(fit_cox_weighted1)$conf.int,2)
round(summary(fit_cox_weighted2)$conf.int,2)




