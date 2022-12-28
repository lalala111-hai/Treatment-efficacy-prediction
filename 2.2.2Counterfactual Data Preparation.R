## 反事实预测结果的反事实数据准备
setwd("E:\\UAI Program\\3-ChuanBeiHospital\\4-刘颖-肝癌治疗方案选择\\0-组学方法做方案选择和生存分析\\3-TreatmentDffect_SurvialAnalysis\\4-治疗效应分析\\4.3 反事实预测")
## 手术组：随机森林模型，预测特征：Survival time、Portal.Vein.Tumor.Thrombus、Cirrhosis、Fusion.lesions、
## Intratumoral.necrosis、Tumor.number、Tumor.diameter、age、ALB、Radscore_A3及Radscore_P

## TACE组：逻辑回归模型，预测特征：Survival time、Portal.Vein.Tumor.Thrombus、ALT、age、ALB及Radscore_A5
library(magrittr)
# 临床指标读取
TrainClinical <- read.csv("3-Clincal_afterSMD-afterIPTW-train-data.csv")[,-31]
TestClinical <- read.csv("3-Clincal_afterSMD-test-data.csv")

table(names(TrainClinical)==names(TestClinical))
Clinical <- rbind(TrainClinical,TestClinical)
# 临床指标保存
ReverseOperation <- Clinical[which(Clinical[,1]==1),c("Status","time","Portal.Vein.Tumor.Thrombus","ALT","age","ALB")]
ReverseTAce <- Clinical[which(Clinical[,1]==0),c("Status","time","Portal.Vein.Tumor.Thrombus","age","ALB","Cirrhosis","Fusion.lesions",
	"Intratumoral.necrosis","Tumor.number","Tumor.diameter")]

### Radscore计算并保存
# 原始特征读取
RadioTrainFiles <- list.files()[4:6]
RadioTestFiles <- list.files()[1:3]
RiomicsFeatures <- list()
for(i in 1:length(RadioTrainFiles)){
	TrainTmp <- read.csv(RadioTrainFiles[i])
	TestTmp <- read.csv(RadioTestFiles[i])
	RiomicsFeatures[[i]] <- rbind(TrainTmp ,TestTmp )
}
names(RiomicsFeatures ) <- c("A3","A5","P")

# 特征选择结果读取
RadioFeatureROITaceA5 <- read.csv("4.2-FeatureSelecton_Radio_Lasso_Tace-A5.csv")
RadioFeatureROILRA3 <- read.csv("4-FeatureSelecton_Radio_RF_LR-A3.csv")
RadioFeatureROILRP <- read.csv("4-FeatureSelecton_Radio_RF_LR-P.csv")

## 反事实Radscore计算
# 手术治疗反事实，手术治疗的患者，采用Tace模型的特征
ReverseOperationA5 <- RiomicsFeatures[[2]][which(RiomicsFeatures[[2]][,1]==1),RadioFeatureROITaceA5[-1,1]]
ReverseOperationA5 <- apply(ReverseOperationA5,2,scale)
ReverseOperationA5Label <- RiomicsFeatures[[2]][which(RiomicsFeatures[[2]][,1]==1),3]
ReverseOperationA5Radscore <- ReverseOperationA5 %*% RadioFeatureROITaceA5[-1,3] + RadioFeatureROITaceA5[1,3]

ReverseOperationA5 <- cbind(ReverseOperation,Radscore_A5=ReverseOperationA5Radscore)
names(ReverseOperationA5)[1] <- "Label"
str(ReverseOperationA5)

# Tace治疗反事实，Tace治疗的患者，采用手术模型的特征
ReverseLRA3 <- RiomicsFeatures[[1]][which(RiomicsFeatures[[1]][,1]==0),RadioFeatureROILRA3[-1,1]]
ReverseLRA3 <- apply(ReverseLRA3,2,scale)
ReverseLRP <- RiomicsFeatures[[3]][which(RiomicsFeatures[[3]][,1]==0),RadioFeatureROILRP[-1,1]]
ReverseLRP <- apply(ReverseLRP,2,scale)

ReverseLRA3Radscore <- ReverseLRA3 %*% RadioFeatureROILRA3[-1,3] + RadioFeatureROILRA3[1,3]
ReverseLRPRadscore <- ReverseLRP %*% RadioFeatureROILRP[-1,3] + RadioFeatureROILRP[1,3]

ReverseLRA3P <- cbind(ReverseTAce,Radscore_A3=ReverseLRA3Radscore,Radscore_P=ReverseLRPRadscore)
names(ReverseLRA3P)[1] <- "Label"
str(ReverseLRA3P)


## Tace数据的顺序是变过了的
TaceOrder <- read.csv("6.1-TaceData_glmnet_Uploader_Model_Clinical+A5.csv")
ReverseLRA3POrder <- c()
for(i in 1:nrow(TaceOrder)){
	Tar <- TaceOrder[i,c(1,3,4,6,5)]
	for(j in 1:nrow(ReverseLRA3P)){
		Tmp <- ReverseLRA3P[j,1:5]
		if(sum(Tar)==sum(Tmp) & sum(Tar[c(2,4,5)])==sum(Tmp[c(2,4,5)])){
			ReverseLRA3POrder <- c(ReverseLRA3POrder,j)
		}
	}
}
length(ReverseLRA3POrder)

ReverseLRA3P <- ReverseLRA3P[ReverseLRA3POrder,]

## 输出数据，并于平台上验证
write.csv(ReverseOperationA5,"7.手术治疗患者的反事实数据.csv",row.names=F)
write.csv(ReverseLRA3P,"7.TAce治疗患者的反事实数据.csv",row.names=F)

## 数据标准化
ReverseOperationA5[,2:7] <- apply(ReverseOperationA5[,2:7],2,scale)
ReverseLRA3P[,2:12] <- apply(ReverseLRA3P[,2:12],2,scale)
## 输出数据，并于平台上验证
write.csv(ReverseOperationA5,"7.手术治疗患者的反事实数据Zscore.csv",row.names=F)
write.csv(ReverseLRA3P,"7.TAce治疗患者的反事实数据Zscore.csv",row.names=F)

