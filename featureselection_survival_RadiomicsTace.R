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

str(TrainTace[[1]][,1:20])

## Tace的患者，重新分组为训练集与验证集
TaceTrain <- list();TaceTest <- list()
for(i in 1:6){
	TrainData <- TrainTace[[i]]
	TestData <- TestTace[[i]]
	AllData <- rbind(TrainData,TestData)
	## 按生存状态分层分组
	DeadRow <- which(AllData[,1]==1)
	SurviRow <- which(AllData[,1]==0)
	set.seed(123)
	DeadSample <- sample(DeadRow,length(DeadRow)*0.7)
	set.seed(123)
	SurviSample <- sample(SurviRow,length(SurviRow)*0.7)
	TaceTrain[[i]] <- AllData[c(DeadSample,SurviSample),]
	TaceTest[[i]] <- AllData[c(setdiff(DeadRow,DeadSample),setdiff(SurviRow,SurviSample)),]
}




##### 特征筛选        ----------------------randomForest----------------------
###  随机森林特征选择  ### 随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
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

VariablesLastTace <- list()
for(i in 1:length(TrainLR)){
	VariablesLastTace[[i]] <- RFclinical(TaceTrain[[i]],80)
}



# 保存特征选择结果，并计算Radscore
Names <- c("A","A3","A5","P","P3","P5")
for(i in 1:length(VariablesLastTace)){
	# Tace
	write.csv(VariablesLastTace[[i]],paste0("4.1-FeatureSelecton_Radio_RF_Tace-",Names[i],".csv"),row.names=F)
	TmpDataTrainTace <- TrainTace[[i]][,c("Status",VariablesLastTace[[i]][,1])]
	TmpDataTestTace <- TestTace[[i]][,c("Status",VariablesLastTace[[i]][,1])]
	Coeff <- VariablesLastTace[[i]][,3]
	RadscoreTrainTace <- apply(TmpDataTrainTace[,-1],2,scale) %*% Coeff
	RadscoreTestTace <- apply(TmpDataTestTace[,-1],2,scale) %*% Coeff
	RadioTaceTrain <- cbind(TmpDataTrainTace,Radscore=RadscoreTrainTace)
	RadioTaceTest <- cbind(TmpDataTestTace,Radscore=RadscoreTestTace)
	RadioTace <- rbind(RadioTaceTrain,RadioTaceTest)
	RadioTace <- cbind(Group=c(rep("Train",nrow(RadioTaceTrain)),rep("Test",nrow(RadioTaceTest))),RadioTace)
	RadioTace <- RadioTace[,c(2,1,3:ncol(RadioTace ))]
	names(RadioTace)[1] <- "Label"
	write.csv(RadioTace,row.names=F,paste0("5.1-RadiomicsData_Tace-",Names[i],".csv"))
}
	



##### Lasso        ----------------------glmnet----------------------
###  Lasso特征选择  ###
DataTrain <- TaceTrain[[1]];Threshold <- 80
Lassoclinical <- function(DataTrain,Threshold){
	RFData <- DataTrain
	# 去除自相关特征
	CorMatrix <- cor(RFData[,-1])
	MultiColin <- findCorrelation(CorMatrix,cutoff = 0.9)
	TarCol <- setdiff(1:(ncol(RFData)-1),MultiColin)
	if(length(TarCol)>0){
		LassoData <- RFData[,c(1,TarCol+1)]
	}
	else{
		LassoData <- RFData
	}

	#RFData[,1] <- as.factor(RFData[,1])
	Setseed <- seq(100,10000,by=100)
	Variables <- c();Coeff <- c()
	for(i in 1:length(Setseed)){
		set.seed(Setseed[i])
		Plotcv <- cv.glmnet(as.matrix(LassoData[,-1]),as.matrix(LassoData[,1]),
			family='binomial',type.measure="auc",nfolds=10)
		Lambda <- Plotcv$lambda.min
		Variable <- coef(Plotcv$glmnet.fit,s=Lambda ,exact = F)@Dimnames[[1]]
		Coeffsave <- coef(Plotcv$glmnet.fit,s=Lambda,exact = F)@x
		coefficients <- coef(Plotcv$glmnet.fit,s=Lambda  ,exact = F)
		COl <- which(coefficients!=0)     #系数不为0的特征索引
		#VariableLast <- data.frame(Factor=Variable[COl],Coeff)
		Variables <- c(Variables,Variable[COl])
		Coeff <- c(Coeff,Coeffsave)
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

VariablesLastTaceLasso <- list()
for(i in 1:length(TaceTrain)){
	VariablesLastTaceLasso[[i]] <- Lassoclinical(TaceTrain[[i]],80)
}



# 保存特征选择结果，并计算Radscore
Names <- c("A","A3","A5","P","P3","P5")
for(i in 1:length(VariablesLastTaceLasso)){
	# Tace
	write.csv(VariablesLastTaceLasso[[i]],paste0("4.2-FeatureSelecton_Radio_Lasso_Tace-",Names[i],".csv"),row.names=F)
	TmpDataTrainTace <- TrainTace[[i]][,c("Status",VariablesLastTaceLasso[[i]][-1,1])]
	TmpDataTestTace <- TestTace[[i]][,c("Status",VariablesLastTaceLasso[[i]][-1,1])]
	Coeff <- VariablesLastTaceLasso[[i]][-1,3]
	Intercept <- VariablesLastTaceLasso[[i]][1,3]
	if(ncol(TmpDataTestTace)>2){
		RadscoreTrainTace <- apply(TmpDataTrainTace[,-1],2,scale) %*% Coeff
		RadscoreTestTace <- apply(TmpDataTestTace[,-1],2,scale) %*% Coeff
	}
	if(ncol(TmpDataTestTace)<=2){
		RadscoreTrainTace <- scale(TmpDataTrainTace[,-1]) %*% Coeff
		RadscoreTestTace <- scale(TmpDataTestTace[,-1]) %*% Coeff
	}
	RadioTaceTrain <- cbind(TmpDataTrainTace,Radscore=RadscoreTrainTace)
	RadioTaceTest <- cbind(TmpDataTestTace,Radscore=RadscoreTestTace)
	RadioTace <- rbind(RadioTaceTrain,RadioTaceTest)
	RadioTace <- cbind(Group=c(rep("Train",nrow(RadioTaceTrain)),rep("Test",nrow(RadioTaceTest))),RadioTace)
	RadioTace <- RadioTace[,c(2,1,3:ncol(RadioTace ))]
	names(RadioTace)[1] <- "Label"
	write.csv(RadioTace,row.names=F,paste0("5.2-RadiomicsDataLasso_Tace-",Names[i],".csv"))
}
	
