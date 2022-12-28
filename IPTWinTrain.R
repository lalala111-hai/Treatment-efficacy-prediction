## 训练集IPTW

setwd("E:\\UAI Program\\3-ChuanBeiHospital\\4-刘颖-肝癌治疗方案选择\\0-组学方法做方案选择和生存分析\\3-TreatmentDffect_SurvialAnalysis\\1-Tableone")

library(tableone);library(magrittr)
library(glmnet)
library(Matching);library(survey)

## 对训练集进行IPTW
TrainDat <- read.csv("3-afterSMD-beforeIPTW-train-data.csv")
str(TrainDat);table(TrainDat[,4])
TrainDat[which(TrainDat[,4]==1),4] <- 1
TrainDat[which(TrainDat[,4]==2),4] <- 1
TrainDat[which(TrainDat[,4]==3),4] <- 2
TrainDat[which(TrainDat[,4]==4),4] <- 2
TrainDat[which(TrainDat[,4]==5),4] <- 3
TrainDat[which(TrainDat[,4]==6),4] <- 3


ClinicalDataNew <- TrainDat[,c(1,4,7,8,9,15,16,26)]
str(ClinicalDataNew);table(ClinicalDataNew[,1])
#ClinicalDataNew[,8] <- scale(ClinicalDataNew[,8])
## 倾向性分数 
Formular <- paste0(names(ClinicalDataNew)[1],"~",paste0(names(ClinicalDataNew)[-1],collapse="+")) %>% as.formula()
PSmodel <- glm(Formular,family  = binomial(link = "logit"),data=ClinicalDataNew)

## Predicted probability of being assigned to Tace and no tace
#ClinicalDataNew$PSTace <- predict(PSmodel,newdata=ClinicalDataNew, type = "response") %>% as.numeric()
#ClinicalDataNew$PSNoTace <- 1 - ClinicalDataNew$PSTace

## IPTW
pRhc <- predict(PSmodel,newdata=ClinicalDataNew, type = "response") %>% as.numeric()


TrainDat$mw1 <- ifelse(TrainDat$Label==1,1/(pRhc),1/(1-pRhc))
## Weighted data
IPTWTrain <- svydesign(ids = ~ 1, data = TrainDat, weights = ~ mw1,nest=T)#[[7]]
## Construct a table (This is a bit slow.)
tabWeighted <- svyCreateTableOne(vars = names(TrainDat)[-c(1,ncol((TrainDat)))], strata = "Label", factorVars=names(TrainDat)[c(2,4:8,12:21)],
	data = IPTWTrain, test = T)
print(tabWeighted, smd = TRUE)

#write.csv(print(tabWeighted, smd = TRUE),"5-afterPS&IPTW-smdTrain-生存分析.csv")


## 未处理数据分布查看
unWeighted <- CreateTableOne(vars = names(TrainDat)[-c(1,ncol((TrainDat)))], strata = "Label", factorVars=names(TrainDat)[c(2,4:8,12:21)],
	data = TrainDat, test = T)
print(unWeighted, smd = TRUE)


## IPTW前后的tableone对比
IPTWCompari <- cbind(print(unWeighted, smd = TRUE),print(tabWeighted, smd = TRUE))
# 插入分组
IPTWCompari <- rbind(Group=rep(c("Operation","Tace","P","Test","SMD"),2),IPTWCompari)
write.csv(IPTWCompari ,"2-smd_after_IPTW-Train.csv")

write.csv(TrainDat,"4-afterSMD-afterIPTW-train-data.csv",row.names=F)


