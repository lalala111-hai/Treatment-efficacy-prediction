# 设置工作目录
setwd("E:\\UAI Program\\3-ChuanBeiHospital\\4-刘颖-肝癌治疗方案选择\\0-组学方法做方案选择和生存分析\\3-TreatmentDffect_SurvialAnalysis\\1-Tableone")

library(tableone);library(stringr);library(car)
library(glmnet);library(magrittr)

TrainDat <- read.csv("3-afterSMD-beforeIPTW-train-data.csv")
TestDat <- read.csv("4-afterSMD-test-data.csv")

AllDataFortableNew <- rbind(TrainDat,TestDat)
str(AllDataFortableNew)
AllDataFortableNew <- AllDataFortableNew[,c(2,3,1,4:30)]
AllDataFortableNew[which(AllDataFortableNew[,4]==1),4] <- 1
AllDataFortableNew[which(AllDataFortableNew[,4]==2),4] <- 1
AllDataFortableNew[which(AllDataFortableNew[,4]==3),4] <- 2
AllDataFortableNew[which(AllDataFortableNew[,4]==4),4] <- 2
AllDataFortableNew[which(AllDataFortableNew[,4]==5),4] <- 3
AllDataFortableNew[which(AllDataFortableNew[,4]==6),4] <- 3

Group <- c(rep("Train",nrow(TrainDat)),rep("Test",nrow(TestDat)))
AllDataFortableNew <- cbind(Group=Group,AllDataFortableNew)


### Table One



TableoneFun <- function(Data,Numeirc,Category){
	NunericData <- Data[,c(1,Numeirc )]
	CategoryData <- Data[,c(1,Category )]
	# 针对数值型变量，行Kolmogorov-Smirnov正太分布检验，符合正太分布，返回均值和标准差，否则计算中位数和上下四分位数
	RowName <- c();ColRes <- c();TestWay <- c()
	Levels <- unique(NunericData[,1])
	for(i in 2:ncol(NunericData)){
		X <- NunericData[which(NunericData[,1]==Levels[1]),i]
		Y <- NunericData[which(NunericData[,1]==Levels[2]),i]
		KS <- shapiro.test(NunericData[,i])	
		LT <- leveneTest(NunericData[,i],NunericData[,1])
		#KS <- ks.test(NunericData[,i],"pnorm")
		Pvalue <- round(KS$p.value ,2);LTPvalue <- round(LT[[3]][1] ,2)
		if(Pvalue>0.05 & LTPvalue>0.05){
			RowName <- c(RowName,paste0(names(NunericData)[i],"(mean±sd)"))
			for(j in 1:length(Levels)){
				X0 <- NunericData[which(NunericData[,1]==Levels[j]),i]
				Res <- paste0(round(mean(X0),3)," ± ",round(sqrt(var(X0)),3))
				ColRes <- c(ColRes,Res)
			}
			Pdiff <- t.test(X,Y)
			#Pdiff <- t.test(NunericData[,1],NunericData[,i])
			ColRes <- c(ColRes,round(Pdiff$p.value,5))
			TestWay <- c(TestWay,"t.test")
		}
		else{
			RowName <- c(RowName,paste0(names(NunericData)[i],"(Median[Q1~Q3])"))
			for(j in 1:length(Levels)){
				X0 <- as.numeric(summary(NunericData[which(NunericData[,1]==Levels[j]),i]))
				Res <- paste0(round(X0[3],3),"[",round(X0[2],3),"~",round(X0[5],3),"]")
				ColRes <- c(ColRes,Res)
			}
			Pdiff <- wilcox.test(X,Y)
			#Pdiff <- wilcox.test(NunericData[,1],NunericData[,i])
			ColRes <- c(ColRes,round(Pdiff$p.value,7))
			TestWay <- c(TestWay,"wilcox.test")
		}
	}

	NumericDF <- as.data.frame(matrix(ColRes ,ncol=length(Levels)+1 ,byrow=T))
	names(NumericDF ) <-c(Levels,"p") ; row.names(NumericDF ) <-RowName 
	NumericDF <- cbind(level=rep(" ",nrow(NumericDF)),NumericDF,TestWay )
	FactorName <- row.names(NumericDF )
	OutPut <- cbind(FactorName ,NumericDF)

	# 先利用tableone包的内容，将各组占比信息整合起来。然后根据每一个指标，计算其样本量和T值，选择卡方检验/矫正卡方检验和Fisher检验  
	#str(CategoryData)
	TableoneCategory <- CreateTableOne(vars=names(CategoryData)[2:ncol(CategoryData)], 
                 #Vector of variables to summarize
                 data = CategoryData,
                 strata=paste0(names(CategoryData)[1]),
                 factorVars=names(CategoryData)[2:ncol(CategoryData)]) 
	Tableone_csv <- print(TableoneCategory,
              exact =names(CategoryData)[2:ncol(CategoryData)],smd=F,
              showAllLevels = TRUE,quote = FALSE,noSpaces = TRUE,printToggle = FALSE)
	Tableone_csv <- as.data.frame(Tableone_csv)
	CategoryDF <- Tableone_csv[,-5]
	CategoryDF <- cbind(CategoryDF ,TestWay=rep("fisher.test",nrow(CategoryDF )))

	# 如果样本量小于40，直接用Tableone_csv的结果
	if(nrow(Data)<40){
		OutPut <- rbind(CategoryDF ,NumericDF)
		OutPut[1,5] <- " "
		FactorName <- row.names(OutPut)
		NullLocate <- which(stringr::str_detect(row.names(OutPut),"X."))
		FactorName[c(3,NullLocate)] <- " "
		OutPut <- cbind(FactorName ,OutPut)
		OutPut 
	
	}
	else{
		TestWayCat <- c() 
		Pcategory <- c()
		#计算每一个指标的理论频数（T），所有T均大于5，则用卡方检验。T处于1-5之间，用连续矫正的卡方检验
		for(i in 2:ncol(CategoryData)){ #
			DataSub <- table(CategoryData[,1],CategoryData[,i])
			Tvalue <- c();TvalueMore2 <- c()
			for(j in 1:nrow(DataSub)){
				for(k in 1:ncol(DataSub)){
					Value <- sum(DataSub[j,])*sum(DataSub[,k])/sum(DataSub)
					Tvalue <- c(Tvalue ,Value )
					Value2 <- sum(DataSub[j,])*sum(DataSub[,k])/(sum(DataSub)*DataSub[j,k])
					TvalueMore2 <- c(TvalueMore2 ,Value2)
				}
			}
			# 查看是否有频数小于5的
			TBigger5 <- which(Tvalue<5)
			TLower1 <- which(Tvalue<1)
			# 2*2,如果没有，用卡方检验
			if(length(Tvalue)<=4 & length(TBigger5)<1){
				ChisqTest <- chisq.test(DataSub,correct = F)
				Pcategory <- c(Pcategory,round(ChisqTest$p.value,5),rep(" ",ncol(DataSub)-1))
				TestWayCat <- c(TestWayCat,"chisq.test",rep(" ",ncol(DataSub)-1))
			}
			# 2*2,如果有，计算连续矫正卡方值及其P值
			if(length(Tvalue)<=4 & length(TBigger5)>=1 & length(TLower1)<1){
				ChisqTest <- chisq.test(DataSub,correct = T)
				Pcategory <- c(Pcategory,round(ChisqTest$p.value,5),rep(" ",ncol(DataSub)-1))
				TestWayCat <- c(TestWayCat,"correction chisq.test",rep(" ",ncol(DataSub)-1))
			}
			# 2*2,如果有，且有低于1的期望频数，计算Fisher
			if(length(Tvalue)<=4 & length(TLower1)>=1){
				ChisqTest <- fisher.test(DataSub)
				Pcategory <- c(Pcategory,round(ChisqTest$p.value,5),rep(" ",ncol(DataSub)-1))
				TestWayCat <- c(TestWayCat,"Fisher.test",rep(" ",ncol(DataSub)-1))
			}

			# > 2*2,卡方检验
			Ratio <- length(TBigger5)/length(Tvalue)
			if(length(Tvalue)>4 & Ratio < 0.2 & length(TLower1)<1){
				ChisqTest <- chisq.test(DataSub,correct = F)
				Pcategory <- c(Pcategory,round(ChisqTest$p.value,5),rep(" ",ncol(DataSub)-1))
				TestWayCat <- c(TestWayCat,"Fisher.test",rep(" ",ncol(DataSub)-1))
			}
			if(length(Tvalue)>4 & Ratio >= 0.2 & length(TLower1)<1){
				TvalueFor <- as.data.frame(matrix(TvalueMore2,ncol=ncol(DataSub),byrow=T))
				RatioHR <- as.matrix(DataSub*log(TvalueFor)) %>% as.numeric()
				NanNum <- which(RatioHR=="NaN")
				if(length(NanNum)>0){
					RatioHR <- RatioHR[-NanNum]
				}
				else{
					RatioHR <- RatioHR
				}

				Static <- -2*sum(RatioHR)
				Pcategory <- c(Pcategory,round(pchisq(Static,ncol(DataSub)-1),5),rep(" ",ncol(DataSub)-1))
				TestWayCat <- c(TestWayCat,"Likelihood ratio.test",rep(" ",ncol(DataSub)-1))
			}
			if(length(Tvalue)>4 & Ratio < 0.2 & length(TLower1)>=1){
				TvalueFor <- as.data.frame(matrix(TvalueMore2,ncol=ncol(DataSub),byrow=T))
				RatioHR <- as.matrix(DataSub*log(TvalueFor)) %>% as.numeric()
				NanNum <- which(RatioHR=="NaN")
				if(length(NanNum)>0){
					RatioHR <- RatioHR[-NanNum]
				}
				else{
					RatioHR <- RatioHR
				}

				Static <- -2*sum(RatioHR)
				Pcategory <- c(Pcategory,round(pchisq(Static,ncol(DataSub)-1),5),rep(" ",ncol(DataSub)-1))
				TestWayCat <- c(TestWayCat,"Likelihood ratio.test",rep(" ",ncol(DataSub)-1))
			}
			if(length(Tvalue)>4 & Ratio >= 0.2 & length(TLower1)>=1){
				TvalueFor <- as.data.frame(matrix(TvalueMore2,ncol=ncol(DataSub),byrow=T))
				RatioHR <- as.matrix(DataSub*log(TvalueFor)) %>% as.numeric()
				NanNum <- which(RatioHR=="NaN")
				if(length(NanNum)>0){
					RatioHR <- RatioHR[-NanNum]
				}
				else{
					RatioHR <- RatioHR
				}

				Static <- -2*sum(RatioHR)
				Pcategory <- c(Pcategory,round(pchisq(Static,ncol(DataSub)-1),5),rep(" ",ncol(DataSub)-1))
				TestWayCat <- c(TestWayCat,"Likelihood ratio.test",rep(" ",ncol(DataSub)-1))
			}


		}
		CategoryDF[2:nrow(CategoryDF),4] <- Pcategory
		CategoryDF[2:nrow(CategoryDF),5] <- TestWayCat
		OutPut <- rbind(CategoryDF[] ,NumericDF)
		OutPut[1,5] <- " "
		FactorName <- row.names(OutPut)
		NullLocate <- which(stringr::str_detect(row.names(OutPut),"X."))
		FactorName[c(3,NullLocate)] <- " "
		OutPut <- cbind(FactorName ,OutPut)
		OutPut 
	}
}
str(ValiData )
Numeirc <- c(2,9,10,11,22:30);Category <- setdiff(1:30,Numeirc)[-1]
	# 训练集结果
	ValiData <- AllDataFortableNew[c(which(AllDataFortableNew[,1]=="Train")),-1]
	Res <- TableoneFun(ValiData,Numeirc,Category);Res 
	write.csv(Res,row.names=F,paste0("4-Tableone_训练集.csv"))
	# 验证集结果
	ValiData <- AllDataFortableNew[c(which(AllDataFortableNew[,1]=="Test")),-1]
	Res <- TableoneFun(ValiData,Numeirc,Category);Res 
	write.csv(Res,row.names=F,paste0("4-Tableone_验证集.csv"))

	# 训练集比验证集结果
	ValiData <- AllDataFortableNew[,-2]
	Res <- TableoneFun(ValiData,Numeirc,Category);Res 
	write.csv(Res,row.names=F,paste0("4-Tableone_训练集VS验证集.csv"))

	# All患者
	ValiData <- AllDataFortableNew[,-1]
	Res <- TableoneFun(ValiData,Numeirc,Category);Res 
	write.csv(Res,row.names=F,paste0("4-Tableone_All.csv"))