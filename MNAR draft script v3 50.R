#--------------------------------------------------------------#
#### MNAR v3 (indicators in outcome simulation) 50% missing ####
#--------------------------------------------------------------#

# loading packages
library(bindata)
library(MASS)
library(tidyverse)
library(lme4)
library(lmerTest)
library(naniar)
library(mice)
library(broom.mixed)
library(pROC)
library(DescTools)
#library(table1)
library(missForest)
library(writexl)

#--------------------------------------------------------------#

# defining expit function, will be used in defining probabilities
expit <- function(x) {
  1 / (1 + exp(-x))
}

#--------------------------------------------------------------#

# setting up matrix that will store results of simulation

n_sims <- 250 # number of simulations

m_results <- matrix(nrow = n_sims, ncol = 15) # defining matrix

colnames(m_results) <- c(
  "PFC Pain Med Ind", "PFC Dep Med Ind", "NRMSE GS Ind", "NRMSE SLB Ind",
  "PFC Pain Med No Ind", "PFC Dep Med No Ind", "NRMSE GS No Ind", "NRMSE SLB No Ind",
  "AUC Ind", "AUC No Ind", "AUC CCA", "Miss Overall", "Sing Ind", "Sing Noind", "Sing CCA"
) # names of columns in matrix

#--------------------------------------------------------------#

#### beginning of loop #####

for (i in 1:n_sims) {
  #--------------------------------------------------------------#
  
  set.seed(124 + i)
  
  # defining number of patients, visits, and total observations
  n_pat <- 250
  n_vis <- 5
  n_tobs <- n_pat * n_vis
  
  #--------------------------------------------------------------#
  
  # simulating birth sex, diabetes, dementia, hypertension, and ui
  
  # binary correlation matrix
  matrix2 <- matrix(c(
    1, -0.15, 0.15, 0.15, 0.15,
    -0.15, 1, 0.15, 0.15, 0.15,
    0.15, 0.15, 1, 0.1, 0.15,
    0.15, 0.15, 0.1, 1, 0.1,
    0.15, 0.15, 0.15, 0.1, 1
  ), ncol = 5)
  
  
  df_first <- as.data.frame(rmvbin(
    n = n_pat, margprob = c(0.5, 0.22, 0.22, 0.7, 0.15), bincorr = matrix2,
    colnames = c("birth_sex", "diab", "dem", "htn", "ui")
  ))
  
  #--------------------------------------------------------------#
  
  # simulating age
  
  # creating the variables n_chronic, a count of how many chronic diseases a patient has
  df_first <- df_first %>% mutate(n_chronic = diab + dem + htn + ui)
  
  # obtaining counts - how many patients have x number of chronic diseases
  n_zero <- length(df_first$n_chronic[df_first$n_chronic == 0])
  n_one <- length(df_first$n_chronic[df_first$n_chronic == 1])
  n_two <- length(df_first$n_chronic[df_first$n_chronic == 2])
  n_three <- length(df_first$n_chronic[df_first$n_chronic == 3])
  n_four <- length(df_first$n_chronic[df_first$n_chronic == 4])
  
  # creating age variable; distribution based on n_chronic
  df_first$age <- ifelse(df_first$n_chronic == 0, rnorm(n = n_zero, mean = 76, sd = 4),
                         ifelse(df_first$n_chronic == 1, rnorm(n = n_one, mean = 78, sd = 4),
                                ifelse(df_first$n_chronic == 2, rnorm(n = n_two, mean = 80, sd = 4),
                                       ifelse(df_first$n_chronic == 3, rnorm(n = n_three, mean = 82, sd = 4),
                                              rnorm(n = n_four, mean = 84, sd = 4)
                                       )
                                )
                         )
  )
  
  
  #--------------------------------------------------------------#
  
  # creates patient id and switching to longitudinal format
  df_first$patient_id <- 1:n_pat
  df_long <- df_first[rep(seq_len(nrow(df_first)), each = n_vis), ]
  df_long$visit <- NA
  df_long$visit <- rep(1:n_vis, n_pat)
  
  #--------------------------------------------------------------#
  
  # simulating BMI
  
  bmi_ran_int_mean <- 0
  bmi_ran_int_sd <- 2
  bmi_ran_error_mean <- 0
  bmi_ran_error_sd <- 1
  
  u_0i_bmi <- rnorm(n = n_pat, mean = bmi_ran_int_mean, sd = bmi_ran_int_sd)
  df_long$u_0i_bmi <- rep(u_0i_bmi, each = n_vis) # random intercept
  
  df_long$e_ij_bmi <- rnorm(n = n_tobs, mean = bmi_ran_error_mean, sd = bmi_ran_error_sd) 
  # random error
  
  df_long$BMI <- NA
  df_long$BMI <- 27 + df_long$u_0i_bmi + 0.01 * df_long$age + 0.9 * df_long$diab +
    0.5 * df_long$htn - 0.3 * df_long$birth_sex + df_long$e_ij_bmi # defining variable
  
  #--------------------------------------------------------------#
  
  # simulating GS
  
  gs_ran_int_mean <- 0
  gs_ran_int_sd <- 0.18
  gs_ran_error_mean <- 0
  gs_ran_error_sd <- 0.09
  
  u_0i_gs <- rnorm(n = n_pat, mean = gs_ran_int_mean, sd = gs_ran_int_sd)
  df_long$u_0i_gs <- rep(u_0i_gs, each = n_vis) # random intercept
  
  df_long$e_ij_gs <- rnorm(n = n_tobs, mean = gs_ran_error_mean, sd = gs_ran_error_sd) 
  # random error
  
  df_long$GS <- NA
  df_long$GS <- 1.35 + df_long$u_0i_gs - 0.001 * df_long$age - 0.01 * df_long$BMI -
    0.2 * df_long$diab - 0.2 * df_long$birth_sex + df_long$e_ij_gs # defining variable

  # setting bounds - only lower bound
  df_long$GS <- ifelse(df_long$GS < 0, 0, df_long$GS)
  
  
  #--------------------------------------------------------------#
  
  # simulating SLB
  
  slb_ran_int_mean <- 0
  slb_ran_int_sd <- 4
  slb_ran_error_mean <- 0
  slb_ran_error_sd <- 2
  
  u_0i_slb <- rnorm(n = n_pat, mean = slb_ran_int_mean, sd = slb_ran_int_sd)
  df_long$u_0i_slb <- rep(u_0i_slb, each = n_vis) # random intercept
  
  
  df_long$e_ij_slb <- rnorm(n = n_tobs, mean = slb_ran_error_mean, sd = slb_ran_error_sd) 
  # random error
  
  df_long$SLB <- NA
  df_long$SLB <- 17 + df_long$u_0i_slb - 0.01 * df_long$age - 0.01 * df_long$BMI -
    0.5 * df_long$diab - 0.5 * df_long$dem - 0.2 * df_long$birth_sex + df_long$e_ij_slb # defining variable
 
  # setting bounds - only lower bound
  df_long$SLB <- ifelse(df_long$SLB < 0, 0, df_long$SLB)
 
  
  #--------------------------------------------------------------#
  
  # removing unnecessary variables
  df_long <- df_long %>% select(
    -n_chronic, -u_0i_slb, -u_0i_gs, -u_0i_bmi,
    -e_ij_slb, -e_ij_gs, -e_ij_bmi
  )
  
  # rounding BMI, GS, and BMI to 2 decimal places
  df_long$age <- round(df_long$age, 0)
  df_long$BMI <- round(df_long$BMI, 2)
  df_long$GS <- round(df_long$GS, 2)
  df_long$SLB <- round(df_long$SLB, 2)
  
  #--------------------------------------------------------------#
  
  # creating pain medication and depression medication binary variables,
  # probabilities depend on age, sex, and diabetes
  
  # defining random intercept parameters
  pain_meds_ran_int_mean <- 0
  pain_meds_ran_int_sd <- .5
  
  dep_meds_ran_int_mean <- 0
  dep_meds_ran_int_sd <- .5
  
  # creating random intercept
  # pain meds
  u_0i_pain_meds <- rnorm(
    n = n_pat, mean = pain_meds_ran_int_mean,
    sd = pain_meds_ran_int_sd
  )
  df_long$u_0i_pain_meds <- rep(u_0i_pain_meds, each = n_vis)
  
  # depression meds
  u_0i_dep_meds <- rnorm(
    n = n_pat, mean = dep_meds_ran_int_mean,
    sd = dep_meds_ran_int_sd
  )
  df_long$u_0i_dep_meds <- rep(u_0i_dep_meds, each = n_vis)
  
  
  df_long$prob_pain_meds <- expit(-1.5 + df_long$u_0i_pain_meds +
                                    0.7 * df_long$birth_sex + 0.01 * df_long$age +
                                    0.9 * df_long$diab) 
  # defining pain medication probability
  
  
  df_long$prob_dep_meds <- expit(-1.5 + df_long$u_0i_dep_meds +
                                   0.7 * df_long$birth_sex + 0.9 * df_long$dem +
                                   0.01 * df_long$age) 
  # defining depression medication probability
  
  
  df_long$pain_meds <- rbinom(n = n_tobs, size = 1, prob = df_long$prob_pain_meds) 
  # defining pain medication variable
  df_long$dep_meds <- rbinom(n = n_tobs, size = 1, prob = df_long$prob_dep_meds) 
  # defining depression medication variable
  
  #--------------------------------------------------------------#
  
  # creating junk variables, parameters chosen at random
  df_long$junk1 <- rnorm(n = n_tobs, mean = 8, sd = 4)
  df_long$junk2 <- rnorm(n = n_tobs, mean = 30, sd = 15)
  df_long$junk3 <- rnorm(n = n_tobs, mean = 100, sd = 10)
  df_long$junk4 <- rnorm(n = n_tobs, mean = 0, sd = 1)
  df_long$junk5 <- rnorm(n = n_tobs, mean = 80, sd = 50)
  
  #--------------------------------------------------------------#
  
  # inducing missingness in level 1 variables
  
  # moving forward with data set called df_MNAR
  df_MNAR <- df_long
  
  df_MNAR$prob_missing_GS = ifelse(df_MNAR$SLB > quantile(df_MNAR$SLB, 0.25), 0, 0.7)
  df_MNAR$prob_missing_SLB = ifelse(df_MNAR$GS > quantile(df_MNAR$GS, 0.25), 0, 0.7)
  df_MNAR$prob_missing_pm = ifelse(df_MNAR$pain_meds == 0, 0, 0.4)
  df_MNAR$prob_missing_dm = ifelse(df_MNAR$dep_meds == 0, 0, 0.4)
  
  # creating four separate indicator variables
  df_MNAR$m_gs <- rbinom(n = n_tobs, size = 1, prob = df_MNAR$prob_missing_GS)
  df_MNAR$m_slb <- rbinom(n = n_tobs, size = 1, prob = df_MNAR$prob_missing_SLB)
  df_MNAR$m_pm <- rbinom(n = n_tobs, size = 1, prob = df_MNAR$prob_missing_pm)
  df_MNAR$m_dm <- rbinom(n = n_tobs, size = 1, prob = df_MNAR$prob_missing_dm)

  #--------------------------------------------------------------#
  
  # simulating outcome
  
  # parameters for random intercept
  falling_ran_int_mean <- 0
  falling_ran_int_sd <- 0.6
  
  
  # creating random intercept
  u_0i_falling <- rnorm(
    n = n_pat, mean = falling_ran_int_mean,
    sd = falling_ran_int_sd
  )
  df_MNAR$u_0i_falling <- rep(u_0i_falling, each = n_vis)
  
  # creating linear predictor
  RHS <- -6.3 + df_MNAR$u_0i_falling +
    0.05 * df_MNAR$age + 0.1 * df_MNAR$birth_sex +
    0.4 * df_MNAR$diab + 0.4 * df_MNAR$dem +
    0.2 * df_MNAR$htn + 0.4 * df_MNAR$ui + 0.1 * df_MNAR$BMI -
    0.4 * df_MNAR$GS - 0.05 * df_MNAR$SLB + 0.2 * df_MNAR$pain_meds +
    0.2 * df_MNAR$dep_meds + 0.3*df_MNAR$m_gs +
    0.3*df_MNAR$m_slb + 0.3*df_MNAR$m_pm + 0.3*df_MNAR$m_dm
  
  
  # defining probability
  df_MNAR$prob_falling <- expit(RHS)
  
  # creating outcome Y
  df_MNAR$Y <- rbinom(n = n_tobs, size = 1, prob = df_MNAR$prob_falling)
  
  
  #--------------------------------------------------------------#
  
  # changing binary variables to be factors
  df_MNAR$birth_sex <- as.factor(df_MNAR$birth_sex)
  df_MNAR$diab <- as.factor(df_MNAR$diab)
  df_MNAR$dem <- as.factor(df_MNAR$dem)
  df_MNAR$htn <- as.factor(df_MNAR$htn)
  df_MNAR$ui <- as.factor(df_MNAR$ui)
  #df_MNAR$Y <- as.factor(df_MNAR$Y)
  df_MNAR$m_gs <- as.factor(df_MNAR$m_gs)
  df_MNAR$m_slb <- as.factor(df_MNAR$m_slb)
  df_MNAR$m_pm <- as.factor(df_MNAR$m_pm)
  df_MNAR$m_dm <- as.factor(df_MNAR$m_dm)
  
  # will be used later for calculating statistics on imputation
  df_MNAR_true <- df_MNAR
  df_MNAR_mis <- df_MNAR
  # move forward with df_MNAR_mis so df_MNAR_true has the actual values
  
  # based on whether the indicators are 1, setting values of gait speed, single leg balance, pain meds, and dep meds to missing
  df_MNAR_mis$GS <- ifelse(df_MNAR_mis$m_gs == 1, NA, df_MNAR_mis$GS)
  df_MNAR_mis$SLB <- ifelse(df_MNAR_mis$m_slb == 1, NA, df_MNAR_mis$SLB)
  df_MNAR_mis$pain_meds <- ifelse(df_MNAR_mis$m_pm == 1, NA, df_MNAR_mis$pain_meds) %>%
    as.factor()
  df_MNAR_mis$dep_meds <- ifelse(df_MNAR_mis$m_dm == 1, NA, df_MNAR_mis$dep_meds) %>%
    as.factor()
  
  # used to calculate overall percent missingness
  miss <- df_MNAR_mis %>% na.omit()
  miss_overall <- 1 - dim(miss)[1] / 1250
  
  
  #--------------------------------------------------------------#
  
  # imputing with indicator variables
  
  # creating separate dataframe for imputation with indicator,
  # removing unnecessary variables
  df_MNAR_ind <- df_MNAR_mis %>% select(
    -u_0i_falling, -prob_falling, -prob_missing_GS,-prob_missing_SLB,-prob_missing_pm,
    -prob_missing_dm,
    -u_0i_pain_meds,  -u_0i_dep_meds, -prob_pain_meds, -prob_dep_meds
  )
  
  
  # doing this to obtain a predictor matrix
  ini_ind <- mice(df_MNAR_ind, m = 1, maxit = 0)
  
  # now we can look at the predictor matrix
  # for our continuous missing variables, want everything to be 2 except the
  # cluster variable (-2) and the variable itself (0)
  # for binary missing variables, want everything as 1 except the variable itself 
  # and cluster variable (-2); also set indicator foe variable to 0
  pred_ind <- ini_ind$predictorMatrix
  
  pred_ind["GS", ] <- c(
    2, 2, 2, 2, 2, 2, -2, 2, 2, 0, 2, 2, 2, 2, 2, 2,
    2, 2, 0, 2, 2, 2, 2
  ) # 2l.norm
  pred_ind["SLB", ] <- c(
    2, 2, 2, 2, 2, 2, -2, 2, 2, 2, 0, 2, 2, 2, 2,
    2, 2, 2, 2, 0, 2, 2, 2
  ) # 2l.norm
  pred_ind["pain_meds", ] <- c(
    1, 1, 1, 1, 1, 1, -2, 1, 1, 1, 1, 0, 1, 1,
    1, 1, 1, 1, 1, 1, 0, 1, 1
  ) # 2l.bin
  pred_ind["dep_meds", ] <- c(
    1, 1, 1, 1, 1, 1, -2, 1, 1, 1, 1, 1, 0, 1,
    1, 1, 1, 1, 1, 1, 1, 0, 1
  ) # 2l.bin
  
  meth_ind <- c(
    "", "", "", "", "", "", "", "", "", "2l.norm", "2l.norm",
    "2l.bin", "2l.bin", "", "", "", "", "", "", "", "", "", ""
  )
  
  # performing the imputation
  imp_ind <- mice(df_MNAR_ind, m = 1, pred = pred_ind, meth = meth_ind, print = FALSE)
  
  # using the complete function to create the dataframe with imputed values
  df_MNAR_ind <- complete(imp_ind)
  
  # cap lower values
  df_MNAR_ind$GS <- ifelse(df_MNAR_ind$GS < 0, 0, df_MNAR_ind$GS)
  df_MNAR_ind$SLB <- ifelse(df_MNAR_ind$SLB < 0, 0, df_MNAR_ind$SLB)
  
  #--------------------------------------------------------------#
  
  # calculating NRMSE and PFC
  
  # setting pain_meds and dep_meds to be factors in df_MNAR,
  # necessary for calculations
  df_MNAR$pain_meds <- as.factor(df_MNAR$pain_meds)
  df_MNAR$dep_meds <- as.factor(df_MNAR$dep_meds)
  
  # PFC for pain peds
  xmis <- df_MNAR_mis %>% select(pain_meds)
  ximp <- df_MNAR_ind %>% select(pain_meds)
  xtrue <- df_MNAR %>% select(pain_meds)
  pfc_pm_ind <- missForest::mixError(xmis = xmis, ximp = ximp, xtrue = xtrue)
  
  # PFC for depression meds
  xmis <- df_MNAR_mis %>% select(dep_meds)
  ximp <- df_MNAR_ind %>% select(dep_meds)
  xtrue <- df_MNAR %>% select(dep_meds)
  pfc_dm_ind <- missForest::mixError(xmis = xmis, ximp = ximp, xtrue = xtrue)
  
  # NRMSE for GS
  xmis <- df_MNAR_mis %>% select(GS)
  ximp <- df_MNAR_ind %>% select(GS)
  xtrue <- df_MNAR %>% select(GS)
  nrmse_gs_ind <- nrmse(xmis = xmis, ximp = ximp, xtrue = xtrue)
  
  # NRMSE for SLB
  xmis <- df_MNAR_mis %>% select(SLB)
  ximp <- df_MNAR_ind %>% select(SLB)
  xtrue <- df_MNAR %>% select(SLB)
  nrmse_slb_ind <- nrmse(xmis = xmis, ximp = ximp, xtrue = xtrue)
  
  #--------------------------------------------------------------#
  
  # running the model with indicator variables
  lme.test.Y_mar_ind <- glmer(
    Y ~ scale(age) + scale(BMI) + scale(GS) + scale(SLB) + birth_sex + diab + dem +
      htn + ui + pain_meds + dep_meds +
      scale(junk1) + scale(junk2) + scale(junk3) + scale(junk4) + scale(junk5) +
      m_gs + m_slb + m_pm + m_dm + (1 | patient_id),
    data = df_MNAR_ind,
    control = glmerControl(optimizer = "bobyqa"),
    family = "binomial"
  )
  
  
  #--------------------------------------------------------------#
  
  # AUC calculation
  prediction_ind <- predict(lme.test.Y_mar_ind, type = "response")
  roc_ind <- roc(df_MNAR_ind$Y, prediction_ind)
  auc_ind <- roc_ind$auc
  
  #--------------------------------------------------------------#
  
  # imputation without indicator
  
  # creating separate dataframe for imputation without indicator,
  # removing unnecessary variables
  df_MNAR_noind <- df_MNAR_mis %>% select(
    -u_0i_falling, -prob_falling, -prob_missing_GS,-prob_missing_SLB,-prob_missing_pm,
    -prob_missing_dm,
    -u_0i_pain_meds,  -u_0i_dep_meds, -prob_pain_meds, -prob_dep_meds,  -m_gs, -m_slb,
    -m_pm, -m_dm
  )
  
  # imputation
  
  # doing this to obtain a predictor matrix
  ini_noind <- mice(df_MNAR_noind, m = 1, maxit = 0)
  
  # same setup as with indicator
  pred_noind <- ini_noind$predictorMatrix
  
  pred_noind["GS", ] <- c(
    2, 2, 2, 2, 2, 2, -2, 2, 2, 0, 2, 2,
    2, 2, 2, 2, 2, 2, 2
  ) # 2l.norm
  pred_noind["SLB", ] <- c(
    2, 2, 2, 2, 2, 2, -2, 2, 2, 2, 0,
    2, 2, 2, 2, 2, 2, 2, 2
  ) # 2l.norm
  pred_noind["pain_meds", ] <- c(
    1, 1, 1, 1, 1, 1, -2, 1, 1,
    1, 1, 0, 1, 1, 1, 1, 1, 1, 1
  ) # 2l.bin
  pred_noind["dep_meds", ] <- c(
    1, 1, 1, 1, 1, 1, -2, 1, 1,
    1, 1, 1, 0, 1, 1, 1, 1, 1, 1
  ) # 2l.bin
  
  meth_noind <- c(
    "", "", "", "", "", "", "", "", "", "2l.norm", "2l.norm",
    "2l.bin", "2l.bin", "", "", "", "", "", ""
  )
  
  # performing the imputation
  imp_noind <- mice(df_MNAR_noind,
                    m = 1, pred = pred_noind, meth = meth_noind,
                    print = FALSE
  )
  
  # using the complete function to create the dataframe with imputed values
  df_MNAR_noind <- complete(imp_noind)
  
  # cap lower values
  df_MNAR_noind$GS <- ifelse(df_MNAR_noind$GS < 0, 0, df_MNAR_noind$GS)
  df_MNAR_noind$SLB <- ifelse(df_MNAR_noind$SLB < 0, 0, df_MNAR_noind$SLB)
  
  #--------------------------------------------------------------#
  
  # calculating NRMSE and PFC
  
  # PFC for pain meds
  xmis <- df_MNAR_mis %>% select(pain_meds)
  ximp <- df_MNAR_noind %>% select(pain_meds)
  xtrue <- df_MNAR %>% select(pain_meds)
  pfc_pm_noind <- missForest::mixError(xmis = xmis, ximp = ximp, xtrue = xtrue)
  
  # PFC for dep meds
  xmis <- df_MNAR_mis %>% select(dep_meds)
  ximp <- df_MNAR_noind %>% select(dep_meds)
  xtrue <- df_MNAR %>% select(dep_meds)
  pfc_dm_noind <- missForest::mixError(xmis = xmis, ximp = ximp, xtrue = xtrue)
  
  # NRMSE for GS
  xmis <- df_MNAR_mis %>% select(GS)
  ximp <- df_MNAR_noind %>% select(GS)
  xtrue <- df_MNAR %>% select(GS)
  nrmse_gs_noind <- nrmse(xmis = xmis, ximp = ximp, xtrue = xtrue)
  
  # NRMSE for SLB
  xmis <- df_MNAR_mis %>% select(SLB)
  ximp <- df_MNAR_noind %>% select(SLB)
  xtrue <- df_MNAR %>% select(SLB)
  nrmse_slb_noind <- nrmse(xmis = xmis, ximp = ximp, xtrue = xtrue)
  
  #--------------------------------------------------------------#
  
  # running model
  lme.test.Y_mar_noind <- glmer(
    Y ~ scale(age) + scale(BMI) + scale(GS) + scale(SLB) + birth_sex + diab + dem +
      htn + ui + pain_meds + dep_meds +
      scale(junk1) + scale(junk2) + scale(junk3) + scale(junk4) + scale(junk5) +
      (1 | patient_id),
    data = df_MNAR_noind,
    control = glmerControl(optimizer = "bobyqa"),
    family = "binomial"
  )
  
  #--------------------------------------------------------------#
  
  # calculating AUC
  prediction_noind <- predict(lme.test.Y_mar_noind, type = "response")
  roc_noind <- roc(df_MNAR_noind$Y, prediction_noind)
  auc_noind <- roc_noind$auc
  
  #--------------------------------------------------------------#
  
  # Complete Case Analysis
  
  # defining dataset with all observations with missing values removed
  df_cca <- na.omit(df_MNAR_mis)
  
  # running model
  lme.test.Y_cca <- glmer(
    Y ~ scale(age) + scale(BMI) + scale(GS) + scale(SLB) + birth_sex + diab +
      dem + htn + ui + pain_meds + dep_meds +
      scale(junk1) + scale(junk2) + scale(junk3) + scale(junk4) + scale(junk5) +
      (1 | patient_id),
    data = df_cca,
    control = glmerControl(optimizer = "bobyqa"),
    family = "binomial"
  )
  
  # calculating AUC
  prediction_cca <- predict(lme.test.Y_cca, type = "response")
  roc_cca <- roc(df_cca$Y, prediction_cca)
  auc_cca <- roc_cca$auc
  
  # singularity
  sing_ind = ifelse(isSingular(lme.test.Y_mar_ind), 1, 0)
  sing_noind = ifelse(isSingular(lme.test.Y_mar_noind), 1, 0)
  sing_cca = ifelse(isSingular(lme.test.Y_cca), 1, 0)
  
  #--------------------------------------------------------------#
  
  # saving results in matrix
  
  m_results[i, 1] <- pfc_pm_ind
  m_results[i, 2] <- pfc_dm_ind
  m_results[i, 3] <- nrmse_gs_ind
  m_results[i, 4] <- nrmse_slb_ind
  m_results[i, 5] <- pfc_pm_noind
  m_results[i, 6] <- pfc_dm_noind
  m_results[i, 7] <- nrmse_gs_noind
  m_results[i, 8] <- nrmse_slb_noind
  m_results[i, 9] <- auc_ind
  m_results[i, 10] <- auc_noind
  m_results[i, 11] <- auc_cca
  m_results[i, 12] <- miss_overall
  m_results[i, 13] <- sing_ind
  m_results[i, 14] <- sing_noind
  m_results[i, 15] <- sing_cca
  
}

#### end of loop ####
m_results <- as.data.frame(m_results)

write.csv(m_results, "MNAR v3 50.csv")

