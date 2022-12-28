setwd("F:/#Users/Desktop/论文资料/论著资料/文章初稿/TACE 和 手术切除（影像组学）/材料/3-TreatmentDffect_SurvialAnalysis1/4-治疗效应分析/4.2 随机森林模型预测")

library(magrittr);library(randomForest)



# 临床特征读取
ClinciaLRTest <- read.csv("5-Clincial_LR_Test.csv")
ClinciaLRTrain <- read.csv("5-Clincial_LR_Train.csv")
ClinciaTaceTest <- read.csv("5-Clincial_Tace_Test.csv")
ClinciaTaceTrain <- read.csv("5-Clincial_Tace_Train.csv")

ClinicalLR <- rbind(ClinciaLRTrain,ClinciaLRTest)
GroupLR <- c(rep("Train",nrow(ClinciaLRTrain)),rep("Test",nrow(ClinciaLRTest)))
ClinicalLR <- cbind(Group=GroupLR,ClinicalLR)

ClinicalTace <- rbind(ClinciaTaceTrain,ClinciaTaceTest)
GroupTace <- c(rep("Train",nrow(ClinciaTaceTrain)),rep("Test",nrow(ClinciaTaceTest)))
ClinicalTace <- cbind(Group=GroupTace,ClinicalTace)
table(ClinicalLR[,1])

str(ClinicalLR)

# 影像组学特征读取
FilesLR <- list.files()[5:10]
FilesTace <- list.files()[11:16]


LRRad <- list()
TaceRad <- list()
for(i in 1:6){
	LRRad[[i]] <- read.csv(FilesLR[i])[,-ncol(read.csv(FilesLR[i]))]
	TaceRad[[i]] <- read.csv(FilesTace[i])[,-ncol(read.csv(FilesTace[i]))]
}

names(LRRad) <- paste0("Radscore_",c("A","A3","A5","P","P3","P5"))
names(TaceRad) <- names(LRRad)

# 特征融合
LRRad[[7]] <- ClinicalLR
TaceRad[[7]] <- ClinicalTace
names(LRRad)[7] <- "Clinical"
names(TaceRad)[7] <- "Clinical"

str(TaceRad[[7]])
# daochu
for(i in 1:length(LRRad)){
	Data <- LRRad[[i]][,c(2,1,3:ncol(LRRad[[i]]))]
	names(Data)[1] <- "Label"
	#write.csv(Data,paste0("LR-features-",names(LRRad)[i],".csv"),row.names=F)
	DataTace <- TaceRad[[i]][,c(2,1,3:ncol(TaceRad[[i]]))]
	names(DataTace)[1] <- "Label"
	#write.csv(DataTace,paste0("Tace-features-",names(TaceRad)[i],".csv"),row.names=F)
}

######### LR 模型建立 ##########
Data <- LRRad;ModleColLR <- c(7)
RFClassfi <- function(Data,ModleColLR){
	if(length(ModleColLR)>1){
		RFData <- Data[[ModleColLR[1]]]
		for(i in 2:length(ModleColLR)){
			RFData <- cbind(RFData,Data[[ModleColLR[i]]][,-c(1:2)])
		}
	}
	else{
		RFData <- Data[[ModleColLR[1]]]
	}
	RFData[,2] <- as.factor(RFData[,2])
	TrainData <- RFData[which(RFData[,1]=="Train"),-1]
	TestData <- RFData[which(RFData[,1]=="Test"),-1]
	set.seed(123)
	rf_ntree <- randomForest(Status~ ., data=TrainData, ntree=100,proximity=TRUE)
	PredictTrain <- predict(rf_ntree, newdata=TrainData,type="prob")[,2]
	TrainTable <- table(TrainData[,1],ifelse(PredictTrain>=0.65,1,0))
	PredictTest <- predict(rf_ntree, newdata=TestData,type="prob")[,2]
	TestTable <- table(TestData[,1],ifelse(PredictTest>=0.65,1,0))
	ACCTrain <- sum(diag(TrainTable))/sum(TrainTable)
	ACCTest <- sum(diag(TestTable))/sum(TestTable)
	data.frame(ACCTrain,ACCTest)
}



ModleColLR <- list(7,1,2,3,4,5,6,
	c(1,4),c(1,5),c(1,6),c(2,4),c(2,5),c(2,6),c(3,4),c(3,5),c(3,6),
	c(1,7),c(2,7),c(3,7),c(4,7),c(5,7),c(6,7),
	c(1,7,4),c(1,7,5),c(1,7,6),
	c(2,7,4),c(2,7,5),c(2,7,6),
	c(3,7,4),c(3,7,5),c(3,7,6))
names(ModleColLR) <- c("Clincial","A","A3","A5","P","P3","P5",
	"A+P","A+P3","A+P5","A3+P","A3+P3","A3+P5","A5+P","A5+P3","A5+P5",
	"Clinical+A","Clinical+A3","Clinical+A5","Clinical+P","Clinical+P3","Clinical+P5",
	"Clinical+A+P","Clinical+A+P3","Clinical+A+P5",
	"Clinical+A3+P","Clinical+A3+P3","Clinical+A3+P5",
	"Clinical+A5+P","Clinical+A5+P3","Clinical+A5+P5")

ACCSingleLR <- data.frame()
for(i in 1:length(ModleColLR)){
	ACCSingleLR <- rbind(ACCSingleLR,RFClassfi(LRRad,ModleColLR[[i]]))
}
row.names(ACCSingleLR) <- names(ModleColLR);ACCSingleLR  

## Tace
ACCSingleTace <- data.frame()
for(i in 1:length(ModleColLR)){
	ACCSingleTace <- rbind(ACCSingleTace,RFClassfi(TaceRad,ModleColLR[[i]]))
}
row.names(ACCSingleTace) <- names(ModleColLR);ACCSingleTace
