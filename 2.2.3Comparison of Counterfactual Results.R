setwd("E:\\UAI Program\\3-ChuanBeiHospital\\4-刘颖-肝癌治疗方案选择\\0-组学方法做方案选择和生存分析\\3-TreatmentDffect_SurvialAnalysis\\4-治疗效应分析\\4.3 反事实预测")

library(magrittr)

# 数据服务
TaceDaat <- read.csv("8.Tace治疗患者的实际与反事实预测结果ForR.csv",encoding="UTF-8")
LRDaat <- read.csv("8.手术治疗患者的实际与反事实预测结果ForR.csv",encoding="UTF-8")

str(LRDaat)
Treatments=c(rep(0,nrow(TaceDaat)),rep(1,nrow(LRDaat)))
AllData <- rbind(TaceDaat,LRDaat)
AllData <- cbind(AllData,Treatments)

## 设置规则，如果实际概率与反事实概率的值均低于0.5，采用原label；如果实际概率比反事实概率高，采用反事实的label
Label <- c(0,1)
ReverseTreattace <- c()
for(i in 1:nrow(AllData[1:151,])){
	Tmp <- AllData[i,c(5,3,4)]
	if(Tmp[2]<=0.9 & Tmp[3]<=0.9){
		ReverseTreattace <- c(ReverseTreattace,0)#Tmp[1] %>% as.numeric()
	}
	if(Tmp[2]>0.9 & Tmp[3]>0.9 & Tmp[2]-Tmp[3]>=0.05){
		ReverseTreattace <- c(ReverseTreattace,setdiff(Label,0))
	}
	if(Tmp[2]>0.9 & Tmp[3]>0.9 & Tmp[3]-Tmp[2]>=0.05){
		ReverseTreattace <- c(ReverseTreattace,0)
	}
	if(Tmp[2]>=0.9 & Tmp[3]<0.9){
		ReverseTreattace <- c(ReverseTreattace,setdiff(Label,0))
	}
	if(Tmp[2]<0.9 & Tmp[3]>=0.9){
		ReverseTreattace <- c(ReverseTreattace,0)
	}
}

as.numeric(ReverseTreattace) %>% length()
table(GT=AllData[1:151,5],Reversetace=as.numeric(ReverseTreattace))
# 手术
Label <- c(0,1)

ReverseTreatLR <- c()
for(i in 152:418){
	Tmp <- AllData[i,c(5,3,4)]
	if(Tmp[2]<=0.4 & Tmp[3]<=0.4){
		ReverseTreatLR <- c(ReverseTreatLR ,Tmp[1])#
		print(paste0(i,"_leibie1"))
	}
	if(Tmp[2]>0.4 & Tmp[3]>0.4 & Tmp[2]-Tmp[3]>=0.1){
		ReverseTreatLR <- c(ReverseTreatLR ,setdiff(Label,Tmp[1]))
		print(paste0(i,"_leibie2"))
	}
	if(Tmp[2]>0.4 & Tmp[3]>0.4 & Tmp[2]-Tmp[3]<0.1){
		ReverseTreatLR <- c(ReverseTreatLR ,Tmp[1])
		print(paste0(i,"_leibie3"))
	}
	if(Tmp[2]>=0.4 & Tmp[3]<0.4){
		ReverseTreatLR <- c(ReverseTreatLR ,setdiff(Label,Tmp[1]))
		print(paste0(i,"_leibie4"))
	}
	if(Tmp[2]<0.4 & Tmp[3]>=0.4){
		ReverseTreatLR <- c(ReverseTreatLR ,Tmp[1])
		print(paste0(i,"_leibie5"))
	}
}
as.numeric(ReverseTreatLR ) %>% length()
table(GT=AllData[152:418,5],Reversetace=as.numeric(ReverseTreatLR))


AllData <- cbind(AllData,ReversTreatments=c(as.numeric(ReverseTreattace),as.numeric(ReverseTreatLR)))
table(原始治疗方案=AllData[,5],反事实治疗方案=AllData[,6])



##### 再加上早期。中期和晚期患者的信息
Clincial <- read.csv("1-Clinical_featuresAddCNLCAll.csv")
str(Clincial)
ClincalAll <- Clincial[,c(5,1,2,3,9,11,26,10,15,16,25)]
ClincalAll[which(ClincalAll[,3]==2),3] <- 0
ClincalAll[which(ClincalAll[,3]==3),3] <- 0
ClincalAll[which(ClincalAll[,3]==1),3] <- 1
str(ClincalAll)
table(ClincalAll[,3])
### 顺序调整
Tace <- read.csv("7.TAce治疗患者的反事实数据.csv")[1:5]
LR <- read.csv("7.手术治疗患者的反事实数据.csv")[c(1:3,5,6)]
TarData <- rbind(Tace,LR)
names(TarData )[1] <- "Status"
str(TarData)
table(TarData[,1]==AllData[,2])

# Tace数据的顺序是变过了的
ReverseLRA3POrder <- c()
for(i in 1:nrow(TarData)){
	Tar <- TarData[i,]
	for(j in 1:nrow(ClincalAll)){
		Tmp <- ClincalAll[j,-c(1:2)]
		if(sum(Tar[c(5,4)])==sum(Tmp[c(5,4)]) & sum(Tar[c(2,4,5)])==sum(Tmp[c(2,4,5)])){#sum(Tar)==sum(Tmp) &  & sum(Tar[c(2,4)])==sum(Tmp[c(2,4)])
			ReverseLRA3POrder <- c(ReverseLRA3POrder,j)
		}
	}
}

length(ReverseLRA3POrder)

ClincalAll <- ClincalAll[ReverseLRA3POrder,]
table(ClincalAll[,3]==TarData[,1])
table(AllData[,2]==TarData[,1])
AllData <- cbind(AllData,CNLC=ClincalAll)

table(AllData[,"CNLC.Portal.Vein.Tumor.Thrombus"] == ClincalAll[,"Portal.Vein.Tumor.Thrombus"])

AllData <- AllData[c(152:418,1:151),c(1,8:10,15,16,14,17,5,6,7)]

write.csv(AllData,"9.治疗效应结果.csv",row.names=F)

