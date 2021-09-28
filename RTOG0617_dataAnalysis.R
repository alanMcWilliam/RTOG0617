#### analysis of RTOG0617 data


library(survival)
library(KMsurv)
library(ggplot2)
library(ggpubr)
library(survminer)
library(dplyr)
library(sjPlot)
library(summarytools)

#############################################
### load in data and select variables for analysis

RTOGpatientgData <- read.csv("C:\\Users\\alan_\\Desktop\\RTOG0617\\dataSheets\\NCT00533949-D1-Dataset.csv")
RTOGheartDose <- read.csv("C:\\Users\\alan_\\Desktop\\RTOG0617\\heartDose.csv")

RTOGpatientgData <- merge(RTOGpatientgData, RTOGheartDose, by = 'patid')

summary(RTOGpatientgData$heartRegion)


view_df(RTOGpatientgData)
View(dfSummary(RTOGpatientgData))

# select variables for multi-variable analysis
# some patients have GTV volume, others ITV - search and keep a combined colum for these
# high dose arms 2 and 4 -- keep high and low dsoe arms for now, ignoring chemo in dose arms currently
# zubrod - all patients 0 or 1 for trial inclusion
##########
# other heart stats available in datasheets - v5_heart and v30_heart
RTOGmv <- RTOGpatientgData %>%
  mutate(tumourVolume = case_when(volume_gtv > 0 ~ volume_gtv,
                                  is.na(volume_gtv) ~ volume_itv)) %>%
  mutate(doseLevel = case_when(arm == 1 ~ 1,# 'low',
                               arm == 2 ~ 2, #'high',
                               arm == 3 ~ 1, #'low',
                               arm == 4 ~ 2 )) %>%# 'high')) %>%
  mutate(chemo = case_when(arm == 1 ~ 0, #'no',
                           arm == 2 ~ 0, #'no',
                           arm == 3 ~ 1, #'yes',
                           arm == 4 ~ 1 )) %>% #'yes')) %>%
  select(patid, survival_months, survival_status, heartRegion, tumourVolume, age, gender, doseLevel, chemo, zubrod, rt_technique, smoke_hx, race, v5_heart, v30_heart, nonsquam_squam, pet_staging, dmean_lung, volume_ptv, grade3_toxicity)
### ajcc_stage_grp (pt IIIA or IIIB)
RTOGmv <- RTOGmv[RTOGmv$race != 4, ]


### some quick summary stats
View(RTOGmv)

summary(RTOGmv)
summary(RTOGmv$tumourVolume)

# dosimetric summary stats
summary(RTOGmv$heartRegion)
summary(RTOGmv$v5_heart)
summary(RTOGmv$v30_heart)
summary(RTOGmv$dmean_lung)

# summary of dataframe
view(dfSummary(RTOGmv))

######################################################################
#### 1. testing the correct functional form of tumour volume for cox model
#### Paper used volume_ptv
#########################################################################

uniCoxPTV <- coxph(Surv(time = survival_months, event = survival_status)~volume_ptv, data = RTOGmv)
uniCoxPTV
uniCoxPTVlog <- coxph(Surv(time = survival_months, event = survival_status)~log(volume_ptv), data = RTOGmv)
uniCoxPTVlog

AICPTV <- AIC(uniCoxPTV)
AICPTV
AICPTVlog <- AIC(uniCoxPTVlog)
AICPTVlog


uniCoxTumour <- coxph(Surv(time = survival_months, event = survival_status)~tumourVolume, data = RTOGmv)
uniCoxTumour
uniCoxTumourlog <- coxph(Surv(time = survival_months, event = survival_status)~log(tumourVolume), data = RTOGmv)
uniCoxTumourlog

AICtumour <- AIC(uniCoxTumour)
AICtumour
AICtumourlog <- AIC(uniCoxTumourlog)
AICtumourlog

# tumour volume and PTV volume plot
plot(RTOGmv$tumourVolume, RTOGmv$volume_ptv)

### plot histogram of PTV and log(PTV) volume
ggplot(RTOGmv, aes(x=volume_ptv)) + 
  geom_histogram(binwidth=25, color="black", fill="white") +
  theme_classic()

ggplot(RTOGmv, aes(x=log(volume_ptv))) + 
  geom_histogram(binwidth=0.1, color="black", fill="white") +
  theme_classic()


### plot histogram of tumour and log(tumour) volume
ggplot(RTOGmv, aes(x=tumourVolume)) + 
  geom_histogram(binwidth=25, color="black", fill="white") +
  theme_classic()

ggplot(RTOGmv, aes(x=log(tumourVolume))) + 
  geom_histogram(binwidth=0.1, color="black", fill="white") +
  theme_classic()

#### test function form for inclusion in the cox model
# ggcoxfunctional does not like missing data - check for complete cases
RTOGmv2 <- RTOGmv %>% 
  mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
  na.omit()
View(RTOGmv2)

# PTV volume
uniCoxFunc_ptv <- coxph(Surv(time = survival_months, event = survival_status)~volume_ptv + log(volume_ptv), data = RTOGmv2)
ggcoxfunctional(uniCoxFunc_ptv, data = RTOGmv2, point.col = "blue", point.alpha = 0.5)
# tumour volume
uniCoxFunc_tumour <- coxph(Surv(time = survival_months, event = survival_status)~tumourVolume + log(tumourVolume), data = RTOGmv2)
ggcoxfunctional(uniCoxFunc_tumour, data = RTOGmv2, point.col = "blue", point.alpha = 0.5)

## test reduction in the deviance of the model
## defined as minus two times the log-likelihood 
## (chi-squared with one degree of freedom  chi2(1)
anova(uniCoxTumour, uniCoxTumourlog, test = "Chisq")
anova(uniCoxPTV, uniCoxPTVlog, test = "Chisq")

#### test difference in multi-variable models
### 1. recreate moel in paper
### remove v5 for clarity here?
multiCox_Orig <- coxph(Surv(time = survival_months, event = survival_status)~factor(doseLevel) + grade3_toxicity + volume_ptv + v5_heart + factor(zubrod) + factor(pet_staging) + gender + factor(nonsquam_squam) + factor(smoke_hx),  data = RTOGmv)
summary(multiCox_Orig)
### 2. test models with log PTV and tumour or log tomur volumes
multiCox_Orig_PTVlog <- coxph(Surv(time = survival_months, event = survival_status)~factor(doseLevel) + grade3_toxicity + log(volume_ptv) + v5_heart + factor(zubrod) + factor(pet_staging) + gender + factor(nonsquam_squam) + factor(smoke_hx),  data = RTOGmv)
summary(multiCox_Orig_PTVlog)

multiCox_Orig_tumour <- coxph(Surv(time = survival_months, event = survival_status)~factor(doseLevel) + grade3_toxicity + tumourVolume + v5_heart + factor(zubrod) + factor(pet_staging) + gender + factor(nonsquam_squam) + factor(smoke_hx),  data = RTOGmv)
summary(multiCox_Orig_tumour)
multiCox_Orig_tumourlog <- coxph(Surv(time = survival_months, event = survival_status)~factor(doseLevel) + grade3_toxicity + log(tumourVolume) + v5_heart + factor(zubrod) + factor(pet_staging) + gender + factor(nonsquam_squam) + factor(smoke_hx),  data = RTOGmv)
summary(multiCox_Orig_tumourlog)

#AIC of each model
AIC_mv_orig <- AIC(multiCox_Orig)
AIC_mv_PTVlog <- AIC(multiCox_Orig_PTVlog)
AIC_mv_tumour <- AIC(multiCox_Orig_tumour)
AIC_mv_tumourlog <- AIC(multiCox_Orig_tumourlog)

AIC_mv_orig
AIC_mv_PTVlog
AIC_mv_tumour
AIC_mv_tumourlog

anova(multiCox_Orig, multiCox_Orig_PTVlog, test = "Chisq")
anova(multiCox_Orig_tumour, multiCox_Orig_tumourlog, test = "Chisq")

#########################################
####### 2. Heart doses
########################################
## need to include lung mean dose

###AGE?
## factor(chemo) +  factor(rt_technique)


## heart stats 
summary(RTOGmv$heartRegion)
summary(RTOGmv$v5_heart)
summary(RTOGmv$v30_heart)

tapply(RTOGmv$tumourVolume, RTOGmv$doseLevel, summary)

tapply(RTOGmv$heartRegion, RTOGmv$doseLevel, summary)
tapply(RTOGmv$v5_heart, RTOGmv$doseLevel, summary)
tapply(RTOGmv$v30_heart, RTOGmv$doseLevel, summary)
tapply(RTOGmv$dmean_lung, RTOGmv$doseLevel, summary)

t.test(RTOGmv$v30_heart~RTOGmv$doseLevel)
t.test(RTOGmv$v5_heart~RTOGmv$doseLevel)
t.test(RTOGmv$heartRegion~RTOGmv$doseLevel)

wilcox.test(RTOGmv$v30_heart~RTOGmv$doseLevel)
wilcox.test(RTOGmv$v5_heart~RTOGmv$doseLevel)
wilcox.test(RTOGmv$heartRegion~RTOGmv$doseLevel)
wilcox.test(RTOGmv$dmean_lung~RTOGmv$doseLevel)


# uni-variable heart analysis
uniCox <- coxph(Surv(time = survival_months, event = survival_status)~heartRegion, data = RTOGmv)
uniCox <- coxph(Surv(time = survival_months, event = survival_status)~v5_heart, data = RTOGmv)
uniCox <- coxph(Surv(time = survival_months, event = survival_status)~v30_heart, data = RTOGmv)
uniCox


multiCox_Orig_tumourlog_heartRegion <- coxph(Surv(time = survival_months, event = survival_status)~dmean_lung + age + factor(doseLevel) + grade3_toxicity + log(tumourVolume) + heartRegion + factor(zubrod) + factor(pet_staging) + gender + factor(nonsquam_squam) + factor(smoke_hx),  data = RTOGmv)
summary(multiCox_Orig_tumourlog_heartRegion)
multiCox_Orig_tumourlog_heartV5 <- coxph(Surv(time = survival_months, event = survival_status)~dmean_lung + age + factor(doseLevel) + grade3_toxicity + log(tumourVolume) + v5_heart + factor(zubrod) + factor(pet_staging) + gender + factor(nonsquam_squam) + factor(smoke_hx),  data = RTOGmv)
summary(multiCox_Orig_tumourlog_heartV5)
multiCox_Orig_tumourlog_heartV30 <- coxph(Surv(time = survival_months, event = survival_status)~dmean_lung + age + factor(doseLevel) + grade3_toxicity + log(tumourVolume) + v30_heart + factor(zubrod) + factor(pet_staging) + gender + factor(nonsquam_squam) + factor(smoke_hx),  data = RTOGmv)
summary(multiCox_Orig_tumourlog_heartV30)

multiCox_Orig_tumourlog_baseline <- coxph(Surv(time = survival_months, event = survival_status)~dmean_lung + age + factor(doseLevel) + grade3_toxicity + log(tumourVolume) + factor(zubrod) + factor(pet_staging) + gender + factor(nonsquam_squam) + factor(smoke_hx),  data = RTOGmv)
summary(multiCox_Orig_tumourlog_baseline)

AIC(multiCox_baseline)
AIC(multiCox_Orig_tumourlog_heartV5)
AIC(multiCox_Orig_tumourlog_heartV30)
AIC(multiCox_Orig_tumourlog_heartRegion)


## all heart volume stats
multiCox_Orig_tumourlog_heartAll <- coxph(Surv(time = survival_months, event = survival_status)~dmean_lung + age + factor(doseLevel) + grade3_toxicity + log(tumourVolume) + heartRegion + v5_heart + v30_heart + factor(zubrod) + factor(pet_staging) + gender + factor(nonsquam_squam) + factor(smoke_hx),  data = RTOGmv)
summary(multiCox_Orig_tumourlog_heartAll)




################################################
### variable selection
### try LASSO / elastic net
################################################
### need to bootstrap...


pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)


#need to select and shape dataframe
RTOGmv_EN <- RTOGmv
RTOGmv_EN$FU <- RTOGmv$survival_months
RTOGmv_EN$stat <- RTOGmv$survival_status
RTOGmv_EN <- RTOGmv_EN %>%
  select(heartRegion, tumourVolume, age, gender,  zubrod, doseLevel,chemo, rt_technique, smoke_hx, race, v5_heart, v30_heart, nonsquam_squam, pet_staging, dmean_lung, stat)
#doseLevel, chemo,
#grade3_toxicity,


##check for complete cases
RTOGmv_EN$comp <- complete.cases(RTOGmv_EN)
RTOGmv_EN <- RTOGmv_EN[RTOGmv_EN$comp == TRUE,]
RTOGmv_EN <- RTOGmv_EN[-ncol(RTOGmv_EN)]
View(RTOGmv_EN)



bootstrap_r <- function(ds, B) {
  ##ds <- all2_EN
  # Preallocate storage for statistics
  boot_stat_LASSO <- matrix(NA, nrow = B, ncol = 16) #number for ncol - 30 sub-structures collecting coefficients
  boot_stat_rigid <- matrix(NA, nrow = B, ncol = 16) #number for ncol - 30 sub-structures collecting coefficients
  boot_stat_elasticNet <- matrix(NA, nrow = B, ncol = 16) #number for ncol - 30 sub-structures collecting coefficients
  
  # Number of observations
  n <- nrow(ds)
  
  # Perform bootstrap
  for(i in seq_len(B)) {
    print(i)
    # Sample initial data
    gen_data <- ds[ sample(nrow(ds), replace=TRUE), ]
    # Calculate sample data mean and SD
    coefOut2 <- performElasticNet(gen_data)
    #print(coefOut2)
    
    
    boot_stat_LASSO[i,] <- coefOut2[1,]
    boot_stat_rigid[i,] <- coefOut2[2,]
    boot_stat_elasticNet[i,] <- coefOut2[3,]
  }
  
  boot_stat <- rbind(boot_stat_LASSO, boot_stat_rigid, boot_stat_elasticNet)
  
  # Return bootstrap result
  return(boot_stat)
}


performElasticNet <- function(dataAll){
  #set coefficents for glm 
  mdlY <- as.factor(as.matrix(dataAll["stat"]))
  mdlX <- as.matrix(dataAll[-ncol(dataAll)])
  
  coefOut <- matrix(NA, nrow = 3, ncol = 16)
  
  #full LASSO
  cv1 <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = 1)
  md1 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 1)
  #print(coef(md1))
  tmp_LASSO <- data.frame(coef.name = dimnames(coef(md1))[[1]], coef.value = matrix(coef(md1)))
  coefOut[1,] <- tmp_LASSO[,2]
  
  #rigid
  md2 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 0) 
  #print(coef(md2))
  tmp_rigid <- data.frame(coef.name = dimnames(coef(md2))[[1]], coef.value = matrix(coef(md2)))
  coefOut[2,] <- tmp_rigid[,2]
  
  #elastic net LASSO
  a <- seq(0.1, 0.9, 0.05) #change final number for fine tuning to be faster was 0.01
  search <- foreach(i = a, .combine = rbind) %dopar% {
    cv <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
  }
  cv3 <- search[search$cvm == min(search$cvm), ]
  md3 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
  #coef(md3)
  tmp_elasticNet <- data.frame(coef.name = dimnames(coef(md3))[[1]], coef.value = matrix(coef(md3)))
  coefOut[3,] <- tmp_elasticNet[,2]
  
  #coefOut <- c(tmp_LASSO[,2], tmp_rigid[,2], tmp_elasticNet[,2])
  
  #print(coefOut)
  return (coefOut)
}

#set a seed for bootstrapping
set.seed(2017)
b = 200
resultsAll <- bootstrap_r(RTOGmv_EN, b)

#split into results for each model, add colum names and convert to data frames
boot_LASSO <- resultsAll[1:b,]
boot_rigid <- resultsAll[(b+1):(2*b),]
boot_elastic <- resultsAll[(2*b+1):nrow(resultsAll),]

colnames(boot_rigid) <- c("Intercept", "heartRegion", "tumourVolume", "age", "gender","zubrod", "doseLevel", "chemo", "rt_technique", "smoke_hx", "race", "v5_heart", "v30_heart", "nonsquam_squam", "pet_staging", "dmean_lung")
colnames(boot_LASSO) <- c("Intercept", "heartRegion", "tumourVolume", "age", "gender","zubrod", "doseLevel", "chemo", "rt_technique", "smoke_hx", "race", "v5_heart", "v30_heart", "nonsquam_squam", "pet_staging", "dmean_lung")
colnames(boot_elastic) <- c("Intercept", "heartRegion", "tumourVolume", "age", "gender","zubrod", "doseLevel", "chemo", "rt_technique", "smoke_hx", "race", "v5_heart", "v30_heart", "nonsquam_squam", "pet_staging", "dmean_lung")

boot_rigid <- as.data.frame(boot_rigid)
boot_LASSO <- as.data.frame(boot_LASSO)
boot_elastic <- as.data.frame(boot_elastic)

View(boot_rigid)
View(boot_LASSO)
View(boot_elastic)


generateStats <- function(df){
  for(i in colnames(df)){
    print(i)
    #print(mean(df[[i]]))
    #print(summary(df[[i]]))
    #print(quantile(df[[i]], probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)))
    print(sum(df[[i]] > 0))
    
    #plot histogram
    ##LASSO y 0, 150, seq -0.05, 0.05
    #  ggplot(data=df, aes(df[[i]])) + 
    #    geom_histogram(breaks=seq(-0.05,0.05, by = 0.001),
    #                 col = "skyblue", fill = "lightblue") +
    #    labs(title = i, x = "coefficent" ) +
    #    ylim(0,50) +
    #    theme(panel.background = element_blank())
    
    #ggsave(paste("C:\\Users\\alan_\\Desktop\\RTOG0617\\bootstrapFigures\\", i, ".jpg", sep=""))
  }
}

generateStats(boot_elastic)

#write/read data in from bootstarpping
write.csv(boot_rigid, "C:\\Users\\alan_\\Desktop\\RTOG0617\\bootstrapped\\bootRigid500.csv")
write.csv(boot_LASSO, "C:\\Users\\alan_\\Desktop\\RTOG0617\\bootstrapped\\bootLASSO500.csv")
write.csv(boot_elastic, "C:\\Users\\alan_\\Desktop\\RTOG0617\\bootstrapped\\bootElastic500.csv")

boot_rigid <- read.csv("C:\\Users\\alan_\\Desktop\\RTOG0617\\bootstrapped\\bootRigid500.csv")
boot_LASSO <- read.csv("C:\\Users\\alan_\\Desktop\\RTOG0617\\bootstrapped\\bootLASSO500.csv")
boot_elastic <- read.csv("C:\\Users\\alan_\\Desktop\\RTOG0617\\bootstrapped\\bootElastic500.csv")




##### without bootstrapping
#set coefficents for glm 
mdlY <- as.factor(as.matrix(RTOGmv_EN["stat"]))
mdlX <- as.matrix(RTOGmv_EN[-ncol(RTOGmv_EN)])


#full LASSO
cv1 <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = 1)
md1 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 1)
coef(md1)
#write.table(md1, "C:\\Users\\alan_\\Desktop\\templateHeart\\results\\md1.txt", sep="\t") 
#tmp <- roc(newY, as.numeric(predict(md1, newX, type = "response")))
#plot(tmp)

#rigid
md2 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv1$lambda.1se, alpha = 0) 
coef(md2)
tmp2 <- roc(newY, as.numeric(predict(md2, newX, type = "response")))
plot(tmp2)
tmp2

#elastic net LASSO
a <- seq(0.1, 0.9, 0.01)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(mdlX, mdlY, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(mdlX, mdlY, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
coef(md3)



################################################
## make some plots

plot(RTOGmv$heartRegion, RTOGmv$v5_heart)
plot(RTOGmv$heartRegion, RTOGmv$v30_heart)
plot(RTOGmv$v5_heart, RTOGmv$v30_heart)


## tumour with heart stats
plot(RTOGmv$v5_heart, log(RTOGmv$tumourVolume))
plot(RTOGmv$v30_heart, log(RTOGmv$tumourVolume))
plot(RTOGmv$heartRegion, log(RTOGmv$tumourVolume))

cor.test(RTOGmv$v5_heart, log(RTOGmv$tumourVolume))
cor.test(RTOGmv$v30_heart, log(RTOGmv$tumourVolume))
cor.test(RTOGmv$heartRegion, log(RTOGmv$tumourVolume))



## lung dose with heart stats
plot(RTOGmv$dmean_lung, RTOGmv$v5_heart)
plot(RTOGmv$dmean_lung, RTOGmv$v30_heart)
plot(RTOGmv$dmean_lung, RTOGmv$heartRegion)

cor.test(RTOGmv$dmean_lung, RTOGmv$v5_heart)
cor.test(RTOGmv$dmean_lung, RTOGmv$v30_heart)
cor.test(RTOGmv$dmean_lung, RTOGmv$heartRegion)


## 1 = 3D conformal, 2 = IMRT
summary(factor(RTOGmv$rt_technique))
tapply(RTOGmv$heartRegion, RTOGmv$rt_technique, summary)

ggplot(RTOGmv, mapping = aes(x = factor(rt_technique), y = heartRegion)) + geom_boxplot() + theme_classic() +
  labs(title = "Rt technique", x = "Technique", y = "Heart Region Dose (Gy)", fill = "") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
  )#legend.text = element_text(size = 20),
#  ) + stat_compare_means(method = "wilcox.test", aes(group = rt_technique, label = paste0("p = ",..p.format..)), label.x = 2.7, label.y = 90, size = 6)
##fill = rt_technique


## boxplot for heart region, high and low dose arms
summary(factor(RTOGmv$doseLevel))
tapply(RTOGmv$heartRegion, RTOGmv$doseLevel, summary)

ggplot(RTOGmv, mapping = aes(x = factor(doseLevel), y = heartRegion)) + geom_boxplot() + theme_classic() +
  labs(title = "Treatment dose arm", x = "Dose arm", y = "Heart Region Dose (Gy)", fill = "") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
  ) + stat_compare_means(method = "wilcox.test", aes(group = doseLevel, label = paste0("p = ",..p.format..)), label.x = 2.7, label.y = 90, size = 6)

## boxplot for heart v5, high and low dose arms
tapply(RTOGmv$v5_heart, RTOGmv$doseLevel, summary)

ggplot(RTOGmv, mapping = aes(x = factor(doseLevel), y = v5_heart)) + geom_boxplot() + theme_classic() +
  labs(title = "Treatment dose arm", x = "Dose arm", y = "Heart v5 (Gy)", fill = "") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
  ) + stat_compare_means(method = "wilcox.test", aes(group = doseLevel, label = paste0("p = ",..p.format..)), label.x = 2.7, label.y = 90, size = 6)


## boxplot for heart v5, high and low dose arms
tapply(RTOGmv$v30_heart, RTOGmv$doseLevel, summary)

ggplot(RTOGmv, mapping = aes(x = factor(doseLevel), y = v30_heart)) + geom_boxplot() + theme_classic() +
  labs(title = "Treatment dose arm", x = "Dose arm", y = "Heart v30 (Gy)", fill = "") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
  ) + stat_compare_means(method = "wilcox.test", aes(group = doseLevel, label = paste0("p = ",..p.format..)), label.x = 2.7, label.y = 90, size = 6)



## boxplot patients dead and alive - heart region dose
tapply(RTOGmv$heartRegion, RTOGmv$survival_status, summary)

ggplot(RTOGmv, mapping = aes(x = factor(survival_status), y = heartRegion)) + geom_boxplot() + theme_classic() +
  labs(title = "Treatment dose arm", x = "Dose arm", y = "Heart Region Dose (Gy)", fill = "") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
  ) + stat_compare_means(method = "wilcox.test", aes(group = survival_status, label = paste0("p = ",..p.format..)), label.x = 2.7, label.y = 90, size = 6)



########################################
plot_models(multiCox_region, multiCox_v5, grid = TRUE)


###########################################
## look at reported toxicities
###### Not that many events and nothing coming out significant here 

reportedTox <- read.csv("C:\\Users\\alan_\\Desktop\\RTOG0617\\dataSheets\\NCT00533949-D2-Dataset.csv")

cardiacTox <- c('Hypotension', 'Hypertension', 'Cardiac disorder', 'Cardiopulmonary arrest', 'Myocardial ischemia', 'Pericardial effusion', 'Pericarditis', 'Cardiac general', 'Cardiac arrhythmia', 'Cardiac disorder', 'Cardiac pain')

reportedCardiac <- filter(reportedTox, toxicity %in% cardiacTox)
View(reportedCardiac)


## multiple entries per patient, first keep highest grade cardiac tox for each patient, but this will keep all entries of same grade
reportedCardiac_sorted <- reportedCardiac %>% 
  group_by(patid) %>%
  #top_n()
  top_n(1, grade)
View(reportedCardiac_sorted)
## next remove duplicate pat with grade so only one entry left per patient
reportedCardiac_sorted <- reportedCardiac_sorted[!duplicated(reportedCardiac_sorted[c("patid","grade")]),]
View(reportedCardiac_sorted)

summary(factor(reportedCardiac_sorted$grade))
reportedCardiac_merged <- merge(reportedCardiac_sorted, RTOGmv, by = 'patid')
reportedCardiac_merged$toxHigh <- reportedCardiac_merged$grade > 2

ggplot(reportedCardiac_merged, mapping = aes(x = factor(grade), y = heartRegion)) + geom_boxplot() + theme_classic() +
  stat_compare_means(method = "anova", aes(group = grade, label = paste0("p = ",..p.format..)), label.x = 2.7, label.y = 90, size = 6)
ggplot(reportedCardiac_merged, mapping = aes(x = factor(grade), y = heartRegion)) + geom_boxplot() + theme_classic() +
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "1")  

summary(reportedCardiac_merged$toxHigh)
ggplot(reportedCardiac_merged, mapping = aes(x = factor(toxHigh), y = heartRegion)) + geom_boxplot() + theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = toxHigh, label = paste0("p = ",..p.format..)), label.x = 2, label.y = 90, size = 6)

card_fit <- glm(toxHigh ~ heartRegion + v5_heart + v30_heart + age + doseLevel + gender + log(tumourVolume) + dmean_lung, data = reportedCardiac_merged, family = binomial)
summary(card_fit)

ggplot(reportedCardiac_merged, mapping = aes(x = factor(grade), y = v5_heart)) + geom_boxplot() + theme_classic()
ggplot(reportedCardiac_merged, mapping = aes(x = factor(grade), y = v30_heart)) + geom_boxplot() + theme_classic()

anova()