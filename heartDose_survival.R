
### playing with heart region data
### question - do patients on the high dose arm do better if we control the heart region dose


summary(RTOGmv$heartRegion)
tapply(RTOGmv$tumourVolume, RTOGmv$doseLevel, summary)

RTOGmv_limited <- RTOGmv %>%
  filter(heartRegion < 600) 

summary(RTOGmv_limited$heartRegion)

doseArm_limited <- survfit(Surv(time = survival_months, event = survival_status)~doseLevel, data = RTOGmv_limited)
ggsurvplot(doseArm_limited, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)

# complete case analysis
RTOGmv_limited$comp <- complete.cases(RTOGmv_limited)
RTOGmv2 <- RTOGmv_limited[RTOGmv_limited$comp == TRUE,]
RTOGmv2 <- RTOGmv2[-ncol(RTOGmv2)]

library(Matching)
# another method using match best approach so far
### 
psmodel<-glm(doseLevel ~  heartRegion + log(tumourVolume) + dmean_lung, family = binomial(), data = RTOGmv2)
summary(psmodel)
pscore<-psmodel$fitted.values        
summary(pscore)

psmatch <- Match(Tr=RTOGmv2$doseLevel, M=1, X= log(pscore), replace = FALSE, caliper = 0.1)
summary(psmatch)
matched_data <- RTOGmv2[unlist(psmatch[c("index.treated" , "index.control")]),]

tapply(matched_data$heartRegion, matched_data$doseLevel, summary)
tapply(matched_data$dmean_lung, matched_data$doseLevel, summary)
tapply(matched_data$tumourVolume, matched_data$doseLevel, summary)
tapply(matched_data$age, matched_data$doseLevel, summary)

doseArm_matched <- survfit(Surv(time = survival_months, event = survival_status)~doseLevel, data = matched_data)
ggsurvplot(doseArm_matched, risk.table = TRUE, conf.int = TRUE, surv.median.line = "hv", pval = TRUE, ncensor.plot = FALSE)


#write.csv(matched_data, file = "MatchedM=1.csv") #csv file with the matched data :)
