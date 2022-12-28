setwd("E:\\UAI Program\\3-ChuanBeiHospital\\4-刘颖-肝癌治疗方案选择\\0-组学方法做方案选择和生存分析\\3-TreatmentDffect_SurvialAnalysis\\3-模型建立")

library(survival);library(Hmisc);library(rms);library(survminer)
library(ggplot2);library(randomForestSRC)
library(mlr3); library(mlr3proba); library(mlr3pipelines)

# 读入数据
AllData <- read.csv(list.files()[2])
#AllData[,1] <- ifelse(AllData[,2]==1,"Resection","Tace") %>% as.factor()
#AllData[,6] <- ifelse(AllData[,6]==1,"Present","Absent") %>% as.factor()
#AllData[,5] <- as.factor(AllData[,5])

CoxTrain <- AllData[which(AllData[,1]=="Train"),-1]
CoxVali <- AllData[which(AllData[,1]=="Test"),-1]

# 权重
Weights <- read.csv(list.files()[1])[,"mw1"]


#### weibull分布及其KM曲线制作，ROC曲线制作
### 中位生存时间计算
Varibales <- list(c(1:3,4:7),c(2:3,8),c(2:3,9),c(2:3,10),c(2:3,11),c(2:3,12),c(2:3,13),
	c(2:3,8,11),c(2:3,8,13),c(2:3,8,12),c(2:3,10,11),c(2:3,10,13),c(2:3,10,12),
	c(1:3,4:7,8),c(1:3,4:7,9),c(1:3,4:7,10),c(1:3,4:7,11),c(1:3,4:7,12),c(1:3,4:7,13),
	c(1:3,4:7,8,11),c(1:3,4:7,8,12),c(1:3,4:7,8,13),c(1:3,4:7,10,11),c(1:3,4:7,10,12),c(1:3,4:7,10,13))
names(Varibales) <- c("Clinical","A","A3","A5","P","P3","P5",
	"A+P","A+P5","A+P3","A5+P","A5+P5","A5+P3",
	"Clinical+A","Clinical+A3","Clinical+A5","Clinical+P","Clinical+P3","Clinical+P5",
	"Clinical+A+P","Clinical+A+P3","Clinical+A+P5","Clinical+A5+P","Clinical+A5+P3","Clinical+A5+P5")

HighMedian <- c();LowMedian <- c();DeadACCTrain <- c();DeadACCTest <- c()
for(i in 1:length(Varibales)){#
	TmpDataTrain <- CoxTrain[,Varibales[[i]]]
	TmpDataTest <- CoxVali[,Varibales[[i]]]
	# 模型建立
	surv.rfsrc <- survreg(Surv(time, Status) ~ ., data=TmpDataTrain, dist="weibull",weights=Weights)
	#ModelDetail <- summary(surv.rfsrc)$table
	#write.csv(ModelDetail,"8-Clincal+P5模型结果.csv")
	# 模型预测结果
	Predict <- predict(surv.rfsrc, newdata = TmpDataTrain)
	PreMedian <- surv.rfsrc$linear.predictors %>% median()
	Risk <- ifelse(surv.rfsrc$linear.predictors >= PreMedian,"Low risk","High risk")
	#table(Risk,CoxTrain[,3])
	## 高风险和低风险生存中位数计算
	HighLow <- cbind(TmpDataTrain[,c("Status","time")],Risk=Risk)
	HighMedian <- c(HighMedian,HighLow[which(HighLow[,3]=="High risk"),2] %>% median())
	LowMedian <- c(LowMedian,HighLow[which(HighLow[,3]=="Low risk"),2] %>% median())
	## KM曲线
	fittmp <- survfit(Surv(time,Status) ~ Risk ,data=HighLow)
	pdf(paste0("7-KM curve-",names(Varibales)[i]," Train.pdf"),family="Times")
		print(ggsurvplot(fittmp,pval=T,conf.int = TRUE,surv.median.line = "hv",size=1.5,risk.table = TRUE))
	dev.off()
	## 测试集结果
	PredictTest <- predict(surv.rfsrc, newdata = TmpDataTest,type="lp")
	#print(PredictTest)
	Risk <- ifelse(PredictTest >= PreMedian,"Low risk","High risk")
	KMData <- cbind(TmpDataTest[,c("Status","time")],Risk=Risk)
	fittmpTest <- survfit(Surv(time,Status) ~ Risk ,data=KMData)
	pdf(paste0("7-KM curve-",names(Varibales)[i]," validation.pdf"),family="Times")
		print(ggsurvplot(fittmpTest,pval=T,conf.int = TRUE,surv.median.line = "hv",size=1.5,risk.table = TRUE))
	dev.off()

	# 死亡预测精度计算
	HighLow[,3] <- ifelse(HighLow[,3]=="High risk",1,0)
	KMData[,3] <- ifelse(KMData[,3]=="High risk",1,0)
	DeadACCTrain <- c(DeadACCTrain,sum(diag(table(HighLow[,1],HighLow[,3])))/sum(table(HighLow[,1],HighLow[,3])))
	DeadACCTest <- c(DeadACCTest,sum(diag(table(KMData[,1],KMData[,3])))/sum(table(KMData[,1],KMData[,3])))
}
Reslut <- data.frame(HighMedian,LowMedian,DeadACCTrain=round(DeadACCTrain,3),DeadACCTest=round(DeadACCTest,3))
row.names(Reslut) <- names(Varibales)
Reslut
write.csv(Reslut,"7-AFT_Accuracy.csv")


### 时间ROC曲线

### ROC曲线
# 1:5年的预测

Variables <- Varibales[[22]]

### 1，2，3年的ROC曲线
#创建和绘制AUC（展示AUC随时间的变化情况）
library(risksetROC)
CombinationModelTrain <- CoxTrain[,Variables]
CombinationModelVali <- CoxVali[,Variables]
str(CombinationModelTrain )
#for(i in c(1,4,5)){
#	CombinationModelTrain[,i] <- factor(CombinationModelTrain[,i],levels=unique(CombinationModelTrain[,i]))
#	CombinationModelVali[,i] <- factor(CombinationModelVali[,i],levels=unique(CombinationModelVali[,i]))
#}

#CombinationModelVali[,6:7] <- apply(CombinationModelVali[,6:7],2,scale)
#CombinationModelTrain[,2] <- CombinationModelTrain[,2]*30
#CombinationModelVali[,2] <- CombinationModelVali[,2]*30

ddist <- datadist(CombinationModelTrain)
options(datadist='ddist')

Formula <- paste0("Surv(time,Status) ~ ",paste0(names(CombinationModelTrain)[-c(1:2)],collapse="+")) %>% as.formula()
fit <- survreg(Surv(time, Status) ~ ., data=ROCTrain, dist="weibull",weights=Weights)
eta <- predict(fit,newdata=CombinationModelTrain,type="lp")#$predicted
etaVali <- predict(fit,newdata=CombinationModelVali,type="lp")#$predicted



library(timeROC);library(survivalROC)
ROC<-timeROC(T=CombinationModelTrain$time,
             delta=CombinationModelTrain$Status,
             marker=eta,
             #other_markers=as.matrix(bp),
             cause=1,
             weighting="marginal",
             times=quantile(CombinationModelTrain$time,probs=seq(0.1,0.9,0.1)),
             ROC = TRUE,
             iid = TRUE)
ROC
confint(ROC)
#通过约登指数获取用于生存分析分组的最佳阈值：（以cut作为fit的阈值将数据分为两组，然后进行生存分析。）
AUCOne <- survivalROC(Stime=CombinationModelTrain$time, status=CombinationModelTrain$Status,     
                    marker=eta,predict.time =365, method="KM")
AUCTwo <- survivalROC(Stime=CombinationModelTrain$time, status=CombinationModelTrain$Status,     
                    marker=eta,predict.time =365*2, method="KM")
AUCThree <- survivalROC(Stime=CombinationModelTrain$time, status=CombinationModelTrain$Status,     
                    marker=eta,predict.time =365*3, method="KM")
AUCFour <- survivalROC(Stime=CombinationModelTrain$time, status=CombinationModelTrain$Status,     
                    marker=eta,predict.time =365*4, method="KM")
AUCFive <- survivalROC(Stime=CombinationModelTrain$time, status=CombinationModelTrain$Status,     
                    marker=eta,predict.time =365*5, method="KM")


pdf("9-3 year ROC train Clinical+A+P5.pdf",family="Times",height=6,width=6)
plot(AUCOne$TP,AUCOne$FP,type="l",col="red",lwd=2,xlab="False Positive Rate",ylab="True Positive Rate")
lines(sort(AUCTwo$TP,decreasing=T),sort(AUCTwo$FP,decreasing=T),type="l",col="blue",lwd=2)
lines(sort(AUCThree$TP,decreasing=T),sort(AUCThree$FP,decreasing=T),type="l",col="green",lwd=2)
lines(sort(AUCFour$TP,decreasing=T),sort(AUCFour$FP,decreasing=T),type="l",col="purple",lwd=2)
lines(sort(AUCFive$TP,decreasing=T),sort(AUCFive$FP,decreasing=T),type="l",col="orange",lwd=2)
abline(a=0,b=1,col="black")
legend("bottomright", #图例的位置
     	 		legend = paste0(paste0(c("One","Two","Three","Four","Five"),"-year AUC="),
				c(round(1-AUCOne$AUC,3),round(1-AUCTwo$AUC,3),round(1-AUCThree$AUC,3),
					round(1-AUCFour$AUC,3),round(1-AUCFive$AUC,3))), #图例文字
     	  		col =c("red","blue","green","purple","orange"), #图例线的颜色，与文字对应
     	  		lwd = 2,#图例中线的粗细
     	  		cex = 1,#图例字体大小
			bty = "n")
dev.off()

AUCOneVali <- survivalROC(Stime=CombinationModelVali$time, status=CombinationModelVali$Status,     
                    marker=etaVali,predict.time =365, method="KM")
AUCTwoVali <- survivalROC(Stime=CombinationModelVali$time, status=CombinationModelVali$Status,     
                    marker=etaVali,predict.time =365*2, method="KM")
AUCThreeVali <- survivalROC(Stime=CombinationModelVali$time, status=CombinationModelVali$Status,     
                    marker=etaVali,predict.time =365*3, method="KM")
AUCFoueVali <- survivalROC(Stime=CombinationModelVali$time, status=CombinationModelVali$Status,     
                    marker=etaVali,predict.time =365*4, method="KM")
AUCFiveVali <- survivalROC(Stime=CombinationModelVali$time, status=CombinationModelVali$Status,     
                    marker=etaVali,predict.time =365*5, method="KM")

pdf("9-3 year ROC Vali Clinical+A+P5.pdf",family="Times",height=6,width=6)
plot(AUCOneVali$TP,AUCOneVali$FP,type="l",col="red",lwd=2,xlab="False Positive Rate",ylab="True Positive Rate")
lines(sort(AUCTwoVali$TP,decreasing=T),sort(AUCTwoVali$FP,decreasing=T),type="l",col="blue",lwd=2)
lines(sort(AUCThreeVali$TP,decreasing=T),sort(AUCThreeVali$FP,decreasing=T),type="l",col="green",lwd=2)
lines(sort(AUCFoueVali$TP,decreasing=T),sort(AUCFoueVali$FP,decreasing=T),type="l",col="purple",lwd=2)
lines(sort(AUCFiveVali$TP,decreasing=T),sort(AUCFiveVali$FP,decreasing=T),type="l",col="orange",lwd=2)
abline(a=0,b=1,col="black")
legend("bottomright", #图例的位置
     	 		legend = paste0(paste0(c("One","Two","Three","Four","Five"),"-year : AUC="),
				c(round(1-AUCOneVali$AUC,3),round(1-AUCTwoVali$AUC,3),round(1-AUCThreeVali$AUC,3),#
				round(1-AUCFoueVali$AUC,3),round(1-AUCFiveVali$AUC,3))), #图例文字,
     	  		col =c("red","blue","green","purple","orange"), #图例线的颜色，与文字对应
     	  		lwd = 2,#图例中线的粗细
     	  		cex = 1,#图例字体大小
			bty = "n")
dev.off()
