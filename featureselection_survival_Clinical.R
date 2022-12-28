## 临床特征筛选，分类模型和生存分析模型
setwd("E:\\UAI Program\\3-ChuanBeiHospital\\4-刘颖-肝癌治疗方案选择\\0-组学方法做方案选择和生存分析\\3-TreatmentDffect_SurvialAnalysis\\4-治疗效应分析\\4.1 特征选择")

library(randomForest);library(magrittr)
library(randomForestSRC) # 随机生存森林模型
library(caret)
## 生存分析特征筛选
TrainData <- read.csv("3-Clincal_afterSMD-afterIPTW-train-data.csv")
TestData <- read.csv("3-Clincal_afterSMD-test-data.csv")
str(TrainData )
names(TrainData) == names(TestData)

# 对训练集的特征做单因素cox回归分析
FeaSelec_Clinical <- TrainData[,c(3,2,1,4:8,12:21,9:11,22:ncol(TrainData))]
FeaSelec_ClinicalTest <- TestData[,c(3,2,1,4:8,12:21,9:11,22:ncol(TestData))]
str(FeaSelec_Clinical)
table(FeaSelec_Clinical[,3])
## NA列
NAInfor <- which(is.na(FeaSelec_Clinical),arr.ind=T)

# 缺失值处理
FeaSelec_Clinical[c(44,56,76,155,260,271),28] <- round(mean(FeaSelec_Clinical[,28],na.rm=T),1)
FeaSelec_Clinical[222,29] <- round(mean(FeaSelec_Clinical[-222,29]),1)
FeaSelec_Clinical[222,30] <- round(mean(FeaSelec_Clinical[-222,30]),2)

#FeaSelec_Clinical[,c(19:30)] <- apply(FeaSelec_Clinical[,c(19:30)],2,scale) %>% as.data.frame()
#FeaSelec_ClinicalTest[,c(19:30)] <- apply(FeaSelec_ClinicalTest[,c(19:30)],2,scale) %>% as.data.frame()
str(FeaSelec_Clinical)


### 两个模型的数据
ResctionTrain <- FeaSelec_Clinical[which(FeaSelec_Clinical[,3]==1),-c(3,4,31)]
TaceTrain <- FeaSelec_Clinical[which(FeaSelec_Clinical[,3]==0),-c(3,4,31)]
ResctionTest <- FeaSelec_ClinicalTest[which(FeaSelec_ClinicalTest[,3]==1),-c(3,4)]
TaceTest <- FeaSelec_ClinicalTest[which(FeaSelec_ClinicalTest[,3]==0),-c(3,4)]

# 第一列因子化
ResctionTrain[,2] <- as.factor(ResctionTrain[,2])
TaceTrain[,2] <- as.factor(TaceTrain[,2])


## Tace的患者，重新分组为训练集与验证集
#for(i in 1:6){
	TrainData <- TaceTrain
	TestData <- TaceTest
	AllData <- rbind(TrainData,TestData)
	## 按生存状态分层分组
	DeadRow <- which(AllData[,2]==1)
	SurviRow <- which(AllData[,2]==0)
	set.seed(123)
	DeadSample <- sample(DeadRow,length(DeadRow)*0.7)
	set.seed(123)
	SurviSample <- sample(SurviRow,length(SurviRow)*0.7)
	TaceTrain <- AllData[c(DeadSample,SurviSample),]
	TaceTest <- AllData[c(setdiff(DeadRow,DeadSample),setdiff(SurviRow,SurviSample)),]
#}

which(is.na(TaceTrain),arr.ind=T)
TaceTrain <- TaceTrain[,-c(20,26,27,28)]
TaceTest <- TaceTest[,-c(20,26,27,28)]
table(TaceTest[,2])

###  随机森林特征选择  ### 随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
DataTrain <- TaceTrain;Threshold <- 80
RFclinical <- function(DataTrain,Threshold){
	RFData <- DataTrain
	RFData <- RFData[,c(2,1,3:ncol(RFData))]
	# 去除自相关特征
	CorMatrix <- cor(RFData[,-1])
	MultiColin <- findCorrelation(CorMatrix,cutoff = 0.7)
	TarCol <- setdiff(1:(ncol(RFData)-1),MultiColin)
	if(length(TarCol)>0){
		RFData <- RFData[,c(1,TarCol+1)]
	}
	else{
		RFData <- RFData
	}
	# 随机森林
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

#ResectionVariables <- RFclinical(ResctionTrain,90)
TaceVariables <- RFclinical(TaceTrain,80)

# 输出结果
Clincial_LR_Train <- ResctionTrain[,c(2,which(names(ResctionTrain) %in% ResectionVariables[,1]))]
Clincial_LR_Test <- ResctionTest[,c(2,which(names(ResctionTrain) %in% ResectionVariables[,1]))]

Clincial_Tace_Train <- TaceTrain[,c(2,which(names(TaceTrain) %in% TaceVariables[,1]))]
Clincial_Tace_Test <- TaceTest[,c(2,which(names(TaceTrain) %in% TaceVariables[,1]))]

write.csv(Clincial_LR_Train,"5-Clincial_LR_Train.csv",row.names=F)
write.csv(Clincial_LR_Test,"5-Clincial_LR_Test.csv",row.names=F)
write.csv(Clincial_Tace_Train,"5.1-Clincial_Tace_Train.csv",row.names=F)
write.csv(Clincial_Tace_Test,"5.1-Clincial_Tace_Test.csv",row.names=F)



### TrainTace的Lasso ####
TaceTrain <- TaceTrain[,c(2,1,3:ncol(TaceTrain))]
TaceTest <- TaceTest[,c(2,1,3:ncol(TaceTest))]

LRData <- TaceTrain#[,c(1:5)]
LRData[,c(2,17:ncol(LRData))] <- apply(LRData[,c(2,17:ncol(LRData))],2,scale)
## 单因素逻辑回归
Bvalue <- c();Orvalue <- c();CIRange <- c();Pvalue <- c()
for(i in 2:ncol(LRData)){
	Formula <- as.formula(paste0("Status~",names(LRData)[i]))
	fit <- glm(Formula,data=LRData,family=binomial(link = "logit"))
	B <- coef(fit) %>% as.data.frame();B <- B[2,1]
	Bvalue <- c(Bvalue ,B)
	Orvalue <- c(Orvalue ,exp(B))
	CI <- exp(confint(fit)) %>% as.data.frame()
	CIRange <- c(CIRange ,paste0(round(CI[-1,1],3),"-",round(CI[-1,2],3)))
	Pvalue <- c(Pvalue ,summary(fit)$coefficients[,4][2] %>% as.numeric())
}

Univariate <- data.frame(Bvalue ,Orvalue ,CIRange,Pvalue)
Univariate <- cbind(Clinical_Factor=names(LRData)[-1],Univariate )
Univariate[,c(2,3,5)] <- round(Univariate[,c(2,3,5)],5)
OR_95CI <- paste0(Univariate[,3],"(",Univariate[,4],")")
Univariate <- cbind(Univariate[,c(1,2,5)],OR_95CI)
write.csv(Univariate ,"Clinical_Selection_Univariate_Allpatients.csv",row.names=F)#+Vali+Test

##### 直接多因素逻辑回归
#MulData <- LRData[,c(1,which(Univariate[,3]<=0.05)+1)];str(MulData)
MulData <- LRData
fit <- glm(Status~.,data=MulData,family=binomial())
Summary <- summary(fit)

Clincial_Tace_Train_LR <- TaceTrain[,c("Status","time","Portal.Vein.Tumor.Thrombus","ALB","age","ALT")]
Clincial_Tace_Test_LR <- TaceTest[,c("Status","time","Portal.Vein.Tumor.Thrombus","ALB","age","ALT")]

write.csv(Clincial_Tace_Train_LR,"5.2-Clincial_Tace_Train_LR.csv",row.names=F)
write.csv(Clincial_Tace_Test_LR,"5.2-Clincial_Tace_Test_LR.csv",row.names=F)







