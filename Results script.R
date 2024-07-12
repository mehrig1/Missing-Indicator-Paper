#--------------------------------------------------------------#
#### Results script
#--------------------------------------------------------------#

library(tidyverse)
library(tidyr)
library(dplyr)
library(gridExtra)
library(ggpubr)

setwd( "C:/Users/mehrig/OneDrive - Wake Forest Baptist Health/Missing data paper files/Results")

############# MNAR v2 ######################

# reading in data from simulation
data20 <- read.csv("MNAR v2 20.csv")
data50 <- read.csv("MNAR v2 50.csv")

mean(data20$Miss.Overall)
sd(data20$Miss.Overall)
mean(data50$Miss.Overall)
sd(data50$Miss.Overall)

table(data20$Sing.CCA)
table(data20$Sing.Ind)
table(data20$Sing.Noind)

# AUC plot

# 20 percent data prep

data_auc20 <- data20 %>% select(10:12) # selecting columns with AUC information

# reshaping data
data_long_auc20 <- gather(data_auc20, method, auc, AUC.Ind:AUC.CCA, factor_key = TRUE)

# creates new data set with mean AUC for each method along with lower and upper bound of
# confidence intervals
data_sum_auc20 <- data_long_auc20 %>%
  group_by(method) %>%
  summarize(
    auc_mean = mean(auc), LB = quantile(auc, probs = 0.025), UB = quantile(auc, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# renaming factor levels of method for plot
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.Ind"] <-
  "Indicator"
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.No.Ind"] <-
  "No Indicator"
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.CCA"] <-
  "Complete Case"

# 50 percent data prep


data_auc50 <- data50 %>% select(10:12) # selecting columns with AUC information

# reshaping data
data_long_auc50 <- gather(data_auc50, method, auc, AUC.Ind:AUC.CCA, factor_key = TRUE)

# creates new data set with mean AUC for each method along with lower and upper bound of
# confidence intervals
data_sum_auc50 <- data_long_auc50 %>%
  group_by(method) %>%
  summarize(
    auc_mean = mean(auc), LB = quantile(auc, probs = 0.025), UB = quantile(auc, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# renaming factor levels of method for plot
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.Ind"] <-
  "Indicator"
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.No.Ind"] <-
  "No Indicator"
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.CCA"] <-
  "Complete Case"


# merge 20 and 50

data_sum_auc = rbind.data.frame(data_sum_auc50, data_sum_auc20)

pd <- position_dodge(.7)

# making plot
plot_auc_MNAR_v2 <- ggplot(data_sum_auc, aes(
  x = method, y = auc_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(auc_mean, 3)), vjust = 1, hjust = -0.2, show.legend = FALSE,
            position = pd, size = 2.75) +
  ylab("AUC") +
  xlab("Method") +
  ggtitle("MNAR: Indicators not included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing")

plot_auc_MNAR_v2

pd <- position_dodge(1)

#--------------------------------------------------------------#

# PFC plot

# 20

data_pfc20 <- data20 %>% select(2:3, 6:7) # selecting columns with PFC

# reshaping data
data_long_pfc20 <- gather(data_pfc20, method, pfc,
                        PFC.Pain.Med.Ind:PFC.Dep.Med.No.Ind,
                        factor_key = TRUE
) # still not in right form for plot

# reshaping again so data is in ideal shape for plot
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Med_type <- c(
  rep("Pain", 50), rep("Depression", 50),
  rep("Pain", 50), rep("Depression", 50)
)
data_long_pfc20 <- cbind.data.frame(Method, Med_type, data_long_pfc20) %>%
  select(-method)

# creates new data set with mean pfc, lower bound, and upper bound for each combination of
# medication and method
data_sum_pfc20 <- data_long_pfc20 %>%
  group_by(Method, Med_type) %>%
  summarize(
    pfc_mean = mean(pfc), LB = quantile(pfc, probs = 0.025), UB = quantile(pfc, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# 50

data_pfc50 <- data50 %>% select(2:3, 6:7) # selecting columns with PFC

# reshaping data
data_long_pfc50 <- gather(data_pfc50, method, pfc,
                          PFC.Pain.Med.Ind:PFC.Dep.Med.No.Ind,
                          factor_key = TRUE
) 
# still not in right form for plot

# reshaping again so data is in ideal shape for plot
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Med_type <- c(
  rep("Pain", 50), rep("Depression", 50),
  rep("Pain", 50), rep("Depression", 50)
)
data_long_pfc50 <- cbind.data.frame(Method, Med_type, data_long_pfc50) %>%
  select(-method)

# creates new data set with mean pfc, lower bound, and upper bound for each combination of
# medication and method
data_sum_pfc50 <- data_long_pfc50 %>%
  group_by(Method, Med_type) %>%
  summarize(
    pfc_mean = mean(pfc), LB = quantile(pfc, probs = 0.025), UB = quantile(pfc, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# merge

data_sum_pfc = rbind.data.frame(data_sum_pfc50, data_sum_pfc20)

# creates plot
plot_pfc_MNAR_v2 <- ggplot(data_sum_pfc, aes(
  x = Method, y = pfc_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(pfc_mean, 3)), vjust = 1, hjust = -0.2, show.legend = FALSE,
            position = pd, size = 2.75) +
  ylab("PFC") +
  xlab("Method") +
  ggtitle("MNAR: Indicators not included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing") +
  facet_wrap(~Med_type)

plot_pfc_MNAR_v2

#--------------------------------------------------------------#

# NRMSE plot

# 20

data_nrmse20 <- data20 %>% select(4:5, 8:9) # selecting columns with NRMSE

# similar to with PFC - reshapes data
data_long_nrmse20 <- gather(data_nrmse20, method, nrmse,
                          NRMSE.GS.Ind:NRMSE.SLB.No.Ind,
                          factor_key = TRUE
)
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Measure <- c(
  rep("GS", 50), rep("SLB", 50),
  rep("GS", 50), rep("SLB", 50)
)
data_long_nrmse20 <- cbind.data.frame(Method, Measure, data_long_nrmse20) %>%
  select(-method)

# creates new data set with mean nrmse, lower bound, and upper bound for each combination of
# variable and method
data_sum_nrmse20 <- data_long_nrmse20 %>%
  group_by(Method, Measure) %>%
  summarize(
    nrmse_mean = mean(nrmse), LB = quantile(nrmse, probs = 0.025),
    UB = quantile(nrmse, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# 50

data_nrmse50 <- data50 %>% select(4:5, 8:9) # selecting columns with NRMSE

# similar to with PFC - reshapes data
data_long_nrmse50 <- gather(data_nrmse50, method, nrmse,
                            NRMSE.GS.Ind:NRMSE.SLB.No.Ind,
                            factor_key = TRUE
)
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Measure <- c(
  rep("GS", 50), rep("SLB", 50),
  rep("GS", 50), rep("SLB", 50)
)
data_long_nrmse50 <- cbind.data.frame(Method, Measure, data_long_nrmse50) %>%
  select(-method)

# creates new data set with mean nrmse, lower bound, and upper bound for each combination of
# variable and method
data_sum_nrmse50 <- data_long_nrmse50 %>%
  group_by(Method, Measure) %>%
  summarize(
    nrmse_mean = mean(nrmse), LB = quantile(nrmse, probs = 0.025),
    UB = quantile(nrmse, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# merge data

data_sum_nrmse = rbind.data.frame(data_sum_nrmse50, data_sum_nrmse20)

# creates plot
plot_nrmse_MNAR_v2 <- ggplot(data_sum_nrmse, aes(
  x = Method, y = nrmse_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(nrmse_mean, 3)), vjust = 1, 
            hjust = -0.2, show.legend = FALSE, position = pd, size = 2.75) +
  ylab("NRMSE") +
  xlab("Method") +
  ggtitle("MNAR: Indicators not included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing") +
  facet_wrap(~Measure)

plot_nrmse_MNAR_v2

data_sum_nrmse <- data_sum_nrmse[, c(1,2,6,3,4,5)]
data_sum_nrmse$nrmse_mean = round(data_sum_nrmse$nrmse_mean, 3)
data_sum_nrmse$LB = round(data_sum_nrmse$LB, 3)
data_sum_nrmse$UB = round(data_sum_nrmse$UB, 3)
colnames(data_sum_nrmse) = c("Method", "Measure", "Missing Percentage", 
                             "NRMSE Mean", "2.5th Quantile", "97.5th Quantile")
MNARv2_tab = flextable(data_sum_nrmse) %>% 
  autofit() %>% bold(part = "header") %>% set_caption(caption = "MNAR v2")

################## MNAR v3 #########################

# reading in data from simulation
data20 <- read.csv("MNAR v3 20.csv")
data50 <- read.csv("MNAR v3 50.csv")

mean(data20$Miss.Overall)
sd(data20$Miss.Overall)
mean(data50$Miss.Overall)
sd(data50$Miss.Overall)

table(data20$Sing.CCA)
table(data20$Sing.Ind)
table(data20$Sing.Noind)

# AUC plot

# 20 percent data prep

data_auc20 <- data20 %>% select(10:12) # selecting columns with AUC information

# reshaping data
data_long_auc20 <- gather(data_auc20, method, auc, AUC.Ind:AUC.CCA, factor_key = TRUE)

# creates new data set with mean AUC for each method along with lower and upper bound of
# confidence intervals
data_sum_auc20 <- data_long_auc20 %>%
  group_by(method) %>%
  summarize(
    auc_mean = mean(auc), LB = quantile(auc, probs = 0.025), UB = quantile(auc, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# renaming factor levels of method for plot
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.Ind"] <-
  "Indicator"
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.No.Ind"] <-
  "No Indicator"
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.CCA"] <-
  "Complete Case"

# 50 percent data prep


data_auc50 <- data50 %>% select(10:12) # selecting columns with AUC information

# reshaping data
data_long_auc50 <- gather(data_auc50, method, auc, AUC.Ind:AUC.CCA, factor_key = TRUE)

# creates new data set with mean AUC for each method along with lower and upper bound of
# confidence intervals
data_sum_auc50 <- data_long_auc50 %>%
  group_by(method) %>%
  summarize(
    auc_mean = mean(auc), LB = quantile(auc, probs = 0.025), UB = quantile(auc, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# renaming factor levels of method for plot
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.Ind"] <-
  "Indicator"
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.No.Ind"] <-
  "No Indicator"
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.CCA"] <-
  "Complete Case"


# merge 20 and 50

data_sum_auc = rbind.data.frame(data_sum_auc50, data_sum_auc20)

pd <- position_dodge(.7)

# making plot
plot_auc_MNAR_v3 <- ggplot(data_sum_auc, aes(
  x = method, y = auc_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(auc_mean, 3)), vjust = 1, hjust = -0.2, show.legend = FALSE,
            position = pd, size = 2.75) +
  ylab("AUC") +
  xlab("Method") +
  ggtitle("MNAR: Indicators included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing")

plot_auc_MNAR_v3

pd <- position_dodge(1)
#--------------------------------------------------------------#

# PFC plot

# 20

data_pfc20 <- data20 %>% select(2:3, 6:7) # selecting columns with PFC

# reshaping data
data_long_pfc20 <- gather(data_pfc20, method, pfc,
                          PFC.Pain.Med.Ind:PFC.Dep.Med.No.Ind,
                          factor_key = TRUE
) # still not in right form for plot

# reshaping again so data is in ideal shape for plot
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Med_type <- c(
  rep("Pain", 50), rep("Depression", 50),
  rep("Pain", 50), rep("Depression", 50)
)
data_long_pfc20 <- cbind.data.frame(Method, Med_type, data_long_pfc20) %>%
  select(-method)

# creates new data set with mean pfc, lower bound, and upper bound for each combination of
# medication and method
data_sum_pfc20 <- data_long_pfc20 %>%
  group_by(Method, Med_type) %>%
  summarize(
    pfc_mean = mean(pfc), LB = quantile(pfc, probs = 0.025), UB = quantile(pfc, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# 50

data_pfc50 <- data50 %>% select(2:3, 6:7) # selecting columns with PFC

# reshaping data
data_long_pfc50 <- gather(data_pfc50, method, pfc,
                          PFC.Pain.Med.Ind:PFC.Dep.Med.No.Ind,
                          factor_key = TRUE
) 
# still not in right form for plot

# reshaping again so data is in ideal shape for plot
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Med_type <- c(
  rep("Pain", 50), rep("Depression", 50),
  rep("Pain", 50), rep("Depression", 50)
)
data_long_pfc50 <- cbind.data.frame(Method, Med_type, data_long_pfc50) %>%
  select(-method)

# creates new data set with mean pfc, lower bound, and upper bound for each combination of
# medication and method
data_sum_pfc50 <- data_long_pfc50 %>%
  group_by(Method, Med_type) %>%
  summarize(
    pfc_mean = mean(pfc), LB = quantile(pfc, probs = 0.025), UB = quantile(pfc, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# merge

data_sum_pfc = rbind.data.frame(data_sum_pfc50, data_sum_pfc20)

# creates plot
plot_pfc_MNAR_v3 <- ggplot(data_sum_pfc, aes(
  x = Method, y = pfc_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(pfc_mean, 3)), vjust = 1, hjust = -0.2, show.legend = FALSE,
            position = pd, size = 2.75) +
  ylab("PFC") +
  xlab("Method") +
  ggtitle("MNAR: Indicators included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing") +
  facet_wrap(~Med_type)

plot_pfc_MNAR_v3

#--------------------------------------------------------------#

# NRMSE plot

# 20

data_nrmse20 <- data20 %>% select(4:5, 8:9) # selecting columns with NRMSE

# similar to with PFC - reshapes data
data_long_nrmse20 <- gather(data_nrmse20, method, nrmse,
                            NRMSE.GS.Ind:NRMSE.SLB.No.Ind,
                            factor_key = TRUE
)
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Measure <- c(
  rep("GS", 50), rep("SLB", 50),
  rep("GS", 50), rep("SLB", 50)
)
data_long_nrmse20 <- cbind.data.frame(Method, Measure, data_long_nrmse20) %>%
  select(-method)

# creates new data set with mean nrmse, lower bound, and upper bound for each combination of
# variable and method
data_sum_nrmse20 <- data_long_nrmse20 %>%
  group_by(Method, Measure) %>%
  summarize(
    nrmse_mean = mean(nrmse), LB = quantile(nrmse, probs = 0.025),
    UB = quantile(nrmse, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# 50

data_nrmse50 <- data50 %>% select(4:5, 8:9) # selecting columns with NRMSE

# similar to with PFC - reshapes data
data_long_nrmse50 <- gather(data_nrmse50, method, nrmse,
                            NRMSE.GS.Ind:NRMSE.SLB.No.Ind,
                            factor_key = TRUE
)
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Measure <- c(
  rep("GS", 50), rep("SLB", 50),
  rep("GS", 50), rep("SLB", 50)
)
data_long_nrmse50 <- cbind.data.frame(Method, Measure, data_long_nrmse50) %>%
  select(-method)

# creates new data set with mean nrmse, lower bound, and upper bound for each combination of
# variable and method
data_sum_nrmse50 <- data_long_nrmse50 %>%
  group_by(Method, Measure) %>%
  summarize(
    nrmse_mean = mean(nrmse), LB = quantile(nrmse, probs = 0.025),
    UB = quantile(nrmse, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# merge data

data_sum_nrmse = rbind.data.frame(data_sum_nrmse50, data_sum_nrmse20)

# creates plot
plot_nrmse_MNAR_v3 <- ggplot(data_sum_nrmse, aes(
  x = Method, y = nrmse_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(nrmse_mean, 3)), vjust = 1, 
            hjust = -0.2, show.legend = FALSE, position = pd, size = 2.75) +
  ylab("NRMSE") +
  xlab("Method") +
  ggtitle("MNAR: Indicators included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing") +
  facet_wrap(~Measure)

plot_nrmse_MNAR_v3

################ MAR v2 ########################

# reading in data from simulation
data20 <- read.csv("MAR v2 20.csv")
data50 <- read.csv("MAR v2 50.csv")

mean(data20$Miss.Overall)
sd(data20$Miss.Overall)
mean(data50$Miss.Overall)
sd(data50$Miss.Overall)

table(data20$Sing.CCA)
table(data20$Sing.Ind)
table(data20$Sing.Noind)

# AUC plot

# 20 percent data prep

data_auc20 <- data20 %>% select(10:12) # selecting columns with AUC information

# reshaping data
data_long_auc20 <- gather(data_auc20, method, auc, AUC.Ind:AUC.CCA, factor_key = TRUE)

# creates new data set with mean AUC for each method along with lower and upper bound of
# confidence intervals
data_sum_auc20 <- data_long_auc20 %>%
  group_by(method) %>%
  summarize(
    auc_mean = mean(auc), LB = quantile(auc, probs = 0.025), UB = quantile(auc, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# renaming factor levels of method for plot
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.Ind"] <-
  "Indicator"
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.No.Ind"] <-
  "No Indicator"
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.CCA"] <-
  "Complete Case"

# 50 percent data prep


data_auc50 <- data50 %>% select(10:12) # selecting columns with AUC information

# reshaping data
data_long_auc50 <- gather(data_auc50, method, auc, AUC.Ind:AUC.CCA, factor_key = TRUE)

# creates new data set with mean AUC for each method along with lower and upper bound of
# confidence intervals
data_sum_auc50 <- data_long_auc50 %>%
  group_by(method) %>%
  summarize(
    auc_mean = mean(auc), LB = quantile(auc, probs = 0.025), UB = quantile(auc, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# renaming factor levels of method for plot
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.Ind"] <-
  "Indicator"
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.No.Ind"] <-
  "No Indicator"
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.CCA"] <-
  "Complete Case"


# merge 20 and 50

data_sum_auc = rbind.data.frame(data_sum_auc50, data_sum_auc20)

pd <- position_dodge(.7)

# making plot
plot_auc_MAR_v2 <- ggplot(data_sum_auc, aes(
  x = method, y = auc_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(auc_mean, 3)), vjust = 1, hjust = -0.2, show.legend = FALSE,
            position = pd, size = 2.75) +
  ylab("AUC") +
  xlab("Method") +
  ggtitle("MAR: Indicators not included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing")

plot_auc_MAR_v2

pd <- position_dodge(1)
#--------------------------------------------------------------#

# PFC plot

# 20

data_pfc20 <- data20 %>% select(2:3, 6:7) # selecting columns with PFC

# reshaping data
data_long_pfc20 <- gather(data_pfc20, method, pfc,
                          PFC.Pain.Med.Ind:PFC.Dep.Med.No.Ind,
                          factor_key = TRUE
) # still not in right form for plot

# reshaping again so data is in ideal shape for plot
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Med_type <- c(
  rep("Pain", 50), rep("Depression", 50),
  rep("Pain", 50), rep("Depression", 50)
)
data_long_pfc20 <- cbind.data.frame(Method, Med_type, data_long_pfc20) %>%
  select(-method)

# creates new data set with mean pfc, lower bound, and upper bound for each combination of
# medication and method
data_sum_pfc20 <- data_long_pfc20 %>%
  group_by(Method, Med_type) %>%
  summarize(
    pfc_mean = mean(pfc), LB = quantile(pfc, probs = 0.025), UB = quantile(pfc, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# 50

data_pfc50 <- data50 %>% select(2:3, 6:7) # selecting columns with PFC

# reshaping data
data_long_pfc50 <- gather(data_pfc50, method, pfc,
                          PFC.Pain.Med.Ind:PFC.Dep.Med.No.Ind,
                          factor_key = TRUE
) 
# still not in right form for plot

# reshaping again so data is in ideal shape for plot
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Med_type <- c(
  rep("Pain", 50), rep("Depression", 50),
  rep("Pain", 50), rep("Depression", 50)
)
data_long_pfc50 <- cbind.data.frame(Method, Med_type, data_long_pfc50) %>%
  select(-method)

# creates new data set with mean pfc, lower bound, and upper bound for each combination of
# medication and method
data_sum_pfc50 <- data_long_pfc50 %>%
  group_by(Method, Med_type) %>%
  summarize(
    pfc_mean = mean(pfc), LB = quantile(pfc, probs = 0.025), UB = quantile(pfc, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# merge

data_sum_pfc = rbind.data.frame(data_sum_pfc50, data_sum_pfc20)

# creates plot
plot_pfc_MAR_v2 <- ggplot(data_sum_pfc, aes(
  x = Method, y = pfc_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(pfc_mean, 3)), vjust = 1, hjust = -0.2, show.legend = FALSE,
            position = pd, size = 2.75) +
  ylab("PFC") +
  xlab("Method") +
  ggtitle("MAR: Indicators not included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing") +
  facet_wrap(~Med_type)

plot_pfc_MAR_v2

#--------------------------------------------------------------#

# NRMSE plot

# 20

data_nrmse20 <- data20 %>% select(4:5, 8:9) # selecting columns with NRMSE

# similar to with PFC - reshapes data
data_long_nrmse20 <- gather(data_nrmse20, method, nrmse,
                            NRMSE.GS.Ind:NRMSE.SLB.No.Ind,
                            factor_key = TRUE
)
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Measure <- c(
  rep("GS", 50), rep("SLB", 50),
  rep("GS", 50), rep("SLB", 50)
)
data_long_nrmse20 <- cbind.data.frame(Method, Measure, data_long_nrmse20) %>%
  select(-method)

# creates new data set with mean nrmse, lower bound, and upper bound for each combination of
# variable and method
data_sum_nrmse20 <- data_long_nrmse20 %>%
  group_by(Method, Measure) %>%
  summarize(
    nrmse_mean = mean(nrmse), LB = quantile(nrmse, probs = 0.025),
    UB = quantile(nrmse, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# 50

data_nrmse50 <- data50 %>% select(4:5, 8:9) # selecting columns with NRMSE

# similar to with PFC - reshapes data
data_long_nrmse50 <- gather(data_nrmse50, method, nrmse,
                            NRMSE.GS.Ind:NRMSE.SLB.No.Ind,
                            factor_key = TRUE
)
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Measure <- c(
  rep("GS", 50), rep("SLB", 50),
  rep("GS", 50), rep("SLB", 50)
)
data_long_nrmse50 <- cbind.data.frame(Method, Measure, data_long_nrmse50) %>%
  select(-method)

# creates new data set with mean nrmse, lower bound, and upper bound for each combination of
# variable and method
data_sum_nrmse50 <- data_long_nrmse50 %>%
  group_by(Method, Measure) %>%
  summarize(
    nrmse_mean = mean(nrmse), LB = quantile(nrmse, probs = 0.025),
    UB = quantile(nrmse, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# merge data

data_sum_nrmse = rbind.data.frame(data_sum_nrmse50, data_sum_nrmse20)

# creates plot
plot_nrmse_MAR_v2 <- ggplot(data_sum_nrmse, aes(
  x = Method, y = nrmse_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(nrmse_mean, 3)), vjust = 1, 
            hjust = -0.2, show.legend = FALSE, position = pd, size = 2.75) +
  ylab("NRMSE") +
  xlab("Method") +
  ggtitle("MAR: Indicators not included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing") +
  facet_wrap(~Measure)

plot_nrmse_MAR_v2


############### MAR v3 #######################

# reading in data from simulation
data20 <- read.csv("MAR v3 20.csv")
data50 <- read.csv("MAR v3 50.csv")

mean(data20$Miss.Overall)
sd(data20$Miss.Overall)
mean(data50$Miss.Overall)
sd(data50$Miss.Overall)

table(data20$Sing.CCA)
table(data20$Sing.Ind)
table(data20$Sing.Noind)

# AUC plot

# 20 percent data prep

data_auc20 <- data20 %>% select(10:12) # selecting columns with AUC information

# reshaping data
data_long_auc20 <- gather(data_auc20, method, auc, AUC.Ind:AUC.CCA, factor_key = TRUE)

# creates new data set with mean AUC for each method along with lower and upper bound of
# confidence intervals
data_sum_auc20 <- data_long_auc20 %>%
  group_by(method) %>%
  summarize(
    auc_mean = mean(auc), LB = quantile(auc, probs = 0.025), UB = quantile(auc, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# renaming factor levels of method for plot
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.Ind"] <-
  "Indicator"
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.No.Ind"] <-
  "No Indicator"
levels(data_sum_auc20$method)[levels(data_sum_auc20$method) == "AUC.CCA"] <-
  "Complete Case"

# 50 percent data prep


data_auc50 <- data50 %>% select(10:12) # selecting columns with AUC information

# reshaping data
data_long_auc50 <- gather(data_auc50, method, auc, AUC.Ind:AUC.CCA, factor_key = TRUE)

# creates new data set with mean AUC for each method along with lower and upper bound of
# confidence intervals
data_sum_auc50 <- data_long_auc50 %>%
  group_by(method) %>%
  summarize(
    auc_mean = mean(auc), LB = quantile(auc, probs = 0.025), UB = quantile(auc, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# renaming factor levels of method for plot
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.Ind"] <-
  "Indicator"
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.No.Ind"] <-
  "No Indicator"
levels(data_sum_auc50$method)[levels(data_sum_auc50$method) == "AUC.CCA"] <-
  "Complete Case"


# merge 20 and 50

data_sum_auc = rbind.data.frame(data_sum_auc50, data_sum_auc20)

pd <- position_dodge(0.7)

# making plot
plot_auc_MAR_v3 <- ggplot(data_sum_auc, aes(
  x = method, y = auc_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(auc_mean, 3)), vjust = 1, hjust = -0.2, show.legend = FALSE,
            position = pd, size = 2.75) +
  ylab("AUC") +
  xlab("Method") +
  ggtitle("MAR: Indicators included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing")

plot_auc_MAR_v3

pd <- position_dodge(1)

#--------------------------------------------------------------#

# PFC plot

# 20

data_pfc20 <- data20 %>% select(2:3, 6:7) # selecting columns with PFC

# reshaping data
data_long_pfc20 <- gather(data_pfc20, method, pfc,
                          PFC.Pain.Med.Ind:PFC.Dep.Med.No.Ind,
                          factor_key = TRUE
) # still not in right form for plot

# reshaping again so data is in ideal shape for plot
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Med_type <- c(
  rep("Pain", 50), rep("Depression", 50),
  rep("Pain", 50), rep("Depression", 50)
)
data_long_pfc20 <- cbind.data.frame(Method, Med_type, data_long_pfc20) %>%
  select(-method)

# creates new data set with mean pfc, lower bound, and upper bound for each combination of
# medication and method
data_sum_pfc20 <- data_long_pfc20 %>%
  group_by(Method, Med_type) %>%
  summarize(
    pfc_mean = mean(pfc), LB = quantile(pfc, probs = 0.025), UB = quantile(pfc, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# 50

data_pfc50 <- data50 %>% select(2:3, 6:7) # selecting columns with PFC

# reshaping data
data_long_pfc50 <- gather(data_pfc50, method, pfc,
                          PFC.Pain.Med.Ind:PFC.Dep.Med.No.Ind,
                          factor_key = TRUE
) 
# still not in right form for plot

# reshaping again so data is in ideal shape for plot
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Med_type <- c(
  rep("Pain", 50), rep("Depression", 50),
  rep("Pain", 50), rep("Depression", 50)
)
data_long_pfc50 <- cbind.data.frame(Method, Med_type, data_long_pfc50) %>%
  select(-method)

# creates new data set with mean pfc, lower bound, and upper bound for each combination of
# medication and method
data_sum_pfc50 <- data_long_pfc50 %>%
  group_by(Method, Med_type) %>%
  summarize(
    pfc_mean = mean(pfc), LB = quantile(pfc, probs = 0.025), UB = quantile(pfc, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# merge

data_sum_pfc = rbind.data.frame(data_sum_pfc50, data_sum_pfc20)

# creates plot
plot_pfc_MAR_v3 <- ggplot(data_sum_pfc, aes(
  x = Method, y = pfc_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(pfc_mean, 3)), vjust = 1, hjust = -0.2, show.legend = FALSE,
            position = pd, size = 2.75) +
  ylab("PFC") +
  xlab("Method") +
  ggtitle("MAR: Indicators included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing") +
  facet_wrap(~Med_type)

plot_pfc_MAR_v3

#--------------------------------------------------------------#

# NRMSE plot

# 20

data_nrmse20 <- data20 %>% select(4:5, 8:9) # selecting columns with NRMSE

# similar to with PFC - reshapes data
data_long_nrmse20 <- gather(data_nrmse20, method, nrmse,
                            NRMSE.GS.Ind:NRMSE.SLB.No.Ind,
                            factor_key = TRUE
)
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Measure <- c(
  rep("GS", 50), rep("SLB", 50),
  rep("GS", 50), rep("SLB", 50)
)
data_long_nrmse20 <- cbind.data.frame(Method, Measure, data_long_nrmse20) %>%
  select(-method)

# creates new data set with mean nrmse, lower bound, and upper bound for each combination of
# variable and method
data_sum_nrmse20 <- data_long_nrmse20 %>%
  group_by(Method, Measure) %>%
  summarize(
    nrmse_mean = mean(nrmse), LB = quantile(nrmse, probs = 0.025),
    UB = quantile(nrmse, probs = 0.975)
  ) %>% mutate(perc_miss = 20)

# 50

data_nrmse50 <- data50 %>% select(4:5, 8:9) # selecting columns with NRMSE

# similar to with PFC - reshapes data
data_long_nrmse50 <- gather(data_nrmse50, method, nrmse,
                            NRMSE.GS.Ind:NRMSE.SLB.No.Ind,
                            factor_key = TRUE
)
Method <- c(rep("Indicator", 100), rep("No Indicator", 100))
Measure <- c(
  rep("GS", 50), rep("SLB", 50),
  rep("GS", 50), rep("SLB", 50)
)
data_long_nrmse50 <- cbind.data.frame(Method, Measure, data_long_nrmse50) %>%
  select(-method)

# creates new data set with mean nrmse, lower bound, and upper bound for each combination of
# variable and method
data_sum_nrmse50 <- data_long_nrmse50 %>%
  group_by(Method, Measure) %>%
  summarize(
    nrmse_mean = mean(nrmse), LB = quantile(nrmse, probs = 0.025),
    UB = quantile(nrmse, probs = 0.975)
  ) %>% mutate(perc_miss = 50)

# merge data

data_sum_nrmse = rbind.data.frame(data_sum_nrmse50, data_sum_nrmse20)

# creates plot
plot_nrmse_MAR_v3 <- ggplot(data_sum_nrmse, aes(
  x = Method, y = nrmse_mean, color = as.factor(perc_miss),
)) +
  geom_point(size = 3, position = pd) +
  geom_linerange(aes(ymin = LB, ymax = UB), position = pd) +
  geom_text(aes(label = round(nrmse_mean, 3)), vjust = 1, 
            hjust = -0.2, show.legend = FALSE, position = pd, size = 2.75) +
  ylab("NRMSE") +
  xlab("Method") +
  ggtitle("MAR: Indicators included in outcome simulation") +
  theme_classic() +
  scale_color_discrete(name = "Percent Missing") +
  facet_wrap(~Measure)

plot_nrmse_MAR_v3


## AUC plot

tiff("auc_plot_test.tiff", units = "in", width = 10, height = 7, res = 300)
ggarrange(plot_auc_MAR_v2, plot_auc_MAR_v3, plot_auc_MNAR_v2, plot_auc_MNAR_v3, 
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

## PFC plot

tiff("auc_pfc_test.tiff", units = "in", width = 10, height = 7, res = 300)
ggarrange(plot_pfc_MAR_v2, plot_pfc_MAR_v3, plot_pfc_MNAR_v2, plot_pfc_MNAR_v3, 
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

## NRMSE plot

tiff("auc_nrmse_test.tiff", units = "in", width = 10, height = 7, res = 300)
ggarrange(plot_nrmse_MAR_v2, plot_nrmse_MAR_v3, plot_nrmse_MNAR_v2, plot_nrmse_MNAR_v3, 
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()