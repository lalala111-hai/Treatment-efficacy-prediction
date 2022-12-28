## TestRows
setwd("E:\\UAI Program\\3-ChuanBeiHospital\\4-刘颖-肝癌治疗方案选择\\0-组学方法做方案选择和生存分析\\1-ForRData\\1-Tableone")
AllData <- read.csv("1-Clinical_features.csv")[,c("BCLC.stage","Age","Survival.time...month","PLT")]
AllData[which(AllData[,1]=="0"),1] <- "A"
AllData <- AllData[which(AllData[,1]!="A"),-1]
TestData <- read.csv("3-afterSMD-test-data.csv")[,2:4]

TestRows <- c()
for(i in 1:nrow(TestData)){
	Tar <- as.numeric(TestData[i,]) 
	for(j in 1:nrow(AllData)){
		Par <- as.numeric(AllData[j,])
		Tab <- which(Tar==Par)
		if(length(Tab)==3){
			TestRows <- c(TestRows,j)
		}
		else{
			next()
		}
	}
}

write.csv(TestRows,"Testrows.csv")