## 临床特征筛选，分类模型和生存分析模型
setwd("E:\\UAI Program\\3-ChuanBeiHospital\\4-刘颖-肝癌治疗方案选择\\0-组学方法做方案选择和生存分析\\3-TreatmentDffect_SurvialAnalysis\\4-治疗效应分析\\4.1 特征选择")

library(glmnet);library(magrittr)
library(survival);library(survminer);library(ezcox)
library(caret)

### 组学特征数据准备
TrainFiles <- list.files()[7:12]
TestFiles <- list.files()[1:6]

TrainLR <- list();TrainTace <- list();TestLR <- list();TestTace <- list()
for(i in 1:length(TrainFiles)){
	TrainData <- read.csv(TrainFiles[i])
	TestData <- read.csv(TestFiles[i])
	TrainLR[[i]] <- TrainData[which(TrainData[,1]==1),-c(1,2)]
	TrainTace[[i]] <- TrainData[which(TrainData[,1]==0),-c(1,2)]
	TestLR[[i]] <- TestData[which(TestData[,1]==1),-c(1,2)]
	TestTace[[i]] <- TestData[which(TestData[,1]==0),-c(1,2)]
}

str(TrainLR[[1]][,1:20])


##### 特征筛选        ----------------------randomForest----------------------
###  随机森林特征选择  ### 随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
DataTrain <- TrainLR[[1]];Threshold <- 70
RFclinical <- function(DataTrain,Threshold){
	RFData <- DataTrain
	# 去除自相关特征
	CorMatrix <- cor(RFData[,-1])
	MultiColin <- findCorrelation(CorMatrix,cutoff = 0.9)
	TarCol <- setdiff(1:(ncol(RFData)-1),MultiColin)
	if(length(TarCol)>0){
		RFData <- RFData[,c(1,TarCol+1)]
	}
	else{
		RFData <- RFData
	}

	RFData[,1] <- as.factor(RFData[,1])
	Setseed <- seq(100,10000,by=100)
	Variables <- c();Coeff <- c()
	for(i in 1:length(Setseed)){
		set.seed(Setseed[i])
		RF <- randomForest(Status~., data =RFData, importance = TRUE)
		Importance <- RF$importance
		variables <- row.names(Importance)[which(round(Importance[,3],3)>0)]
		Variables <- c(Variables,variables)
		Coeff <- c(Coeff,Importance[which(round(Importance[,3],3)>0),2])
	}
	VariablesLast <- data.frame(Variables,Coeff)
	VariablesUnique <- unique(VariablesLast[,1]);CoefMean <- c();OccurNum <- c()
	for(i in 1:length(VariablesUnique)){
		OccurNum <- c(OccurNum,which(VariablesLast[,1]==VariablesUnique[i]) %>% length())
		CoefMean <- c(CoefMean,VariablesLast[which(VariablesLast[,1]==VariablesUnique[i]),2] %>% mean())
	}
	VariablesLast <- data.frame(VariablesUnique,OccurNum,CoefMean)
	VariablesLast <- VariablesLast[which(VariablesLast[,2]>=Threshold),]
	VariablesLast 
}

VariablesLastLR <- list();VariablesLastTace <- list()
for(i in 1:length(TrainLR)){
	VariablesLastLR[[i]] <- RFclinical(TrainLR[[i]],80)
	VariablesLastTace[[i]] <- RFclinical(TrainTace[[i]],80)
}


# 保存特征选择结果，并计算Radscore
Names <- c("A","A3","A5","P","P3","P5")
for(i in 1:length(VariablesLastLR)){
	# LR
	write.csv(VariablesLastLR[[i]],paste0("4-FeatureSelecton_Radio_RF_LR-",Names[i],".csv"),row.names=F)
	TmpDataTrainLR <- TrainLR[[i]][,c("Status",VariablesLastLR[[i]][,1])]
	TmpDataTestLR <- TestLR[[i]][,c("Status",VariablesLastLR[[i]][,1])]
	Coeff <- VariablesLastLR[[i]][,3]
	RadscoreTrainLR <- apply(TmpDataTrainLR[,-1],2,scale) %*% Coeff
	RadscoreTestLR <- apply(TmpDataTestLR[,-1],2,scale) %*% Coeff
	RadioLRTrain <- cbind(TmpDataTrainLR,Radscore=RadscoreTrainLR)
	RadioLRTest <- cbind(TmpDataTestLR,Radscore=RadscoreTestLR)
	RadioLR <- rbind(RadioLRTrain,RadioLRTest)
	RadioLR <- cbind(Group=c(rep("Train",nrow(RadioLRTrain)),rep("Test",nrow(RadioLRTest))),RadioLR)
	write.csv(RadioLR,row.names=F,paste0("5-RadiomicsData_LR-",Names[i],".csv"))
	# Tace
	write.csv(VariablesLastTace[[i]],paste0("4-FeatureSelecton_Radio_RF_Tace-",Names[i],".csv"),row.names=F)
	TmpDataTrainTace <- TrainTace[[i]][,c("Status",VariablesLastTace[[i]][,1])]
	TmpDataTestTace <- TestTace[[i]][,c("Status",VariablesLastTace[[i]][,1])]
	Coeff <- VariablesLastTace[[i]][,3]
	RadscoreTrainTace <- apply(TmpDataTrainTace[,-1],2,scale) %*% Coeff
	RadscoreTestTace <- apply(TmpDataTestTace[,-1],2,scale) %*% Coeff
	RadioTaceTrain <- cbind(TmpDataTrainTace,Radscore=RadscoreTrainTace)
	RadioTaceTest <- cbind(TmpDataTestTace,Radscore=RadscoreTestTace)
	RadioTace <- rbind(RadioTaceTrain,RadioTaceTest)
	RadioTace <- cbind(Group=c(rep("Train",nrow(RadioTaceTrain)),rep("Test",nrow(RadioTaceTest))),RadioTace)
	write.csv(RadioTace,row.names=F,paste0("5-RadiomicsData_Tace-",Names[i],".csv"))
}
	


