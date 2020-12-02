#data-preparation
setwd("/Users/Aleksei/Desktop/Barcelona/Data/ttest/sezerman")       #set correct wd
set.seed(42)
library(caret)

SZ19original<-readMat("19_subset.mat")                              #read file (library(R.matlab))
dfSZ18orig<-as.data.frame(SZ18original)

dfSZ19orig$V1<-as.factor(dfSZ19orig$V1)
inTraindfSZ19<-createDataPartition(dfSZ19orig$V1,p=0.8)[[1]]
testingdfSZ19<-dfSZ19orig[-inTraindfSZ19,]
trainingdfSZ19<-dfSZ19orig[inTraindfSZ19,]
dfRFfitSZ19<-train(V1~.,data=trainingdfSZ19,method = "rf",importance=TRUE)
varImpPlot(dfRFfitSZ19$finalModel)
predictdfSZ19<-predict(dfRFfitSZ19,newdata = testingdfSZ19)
confusionMatrix(testingdfSZ19$V1,predictdfSZ19)
ConMatdfSZ19<-confusionMatrix(testingdfSZ19$V1,predictdfSZ19)
save(dfRFfitSZ19,file="dfRFfitSZ19.rda")
save(ConMatdfSZ19,file="ConMatdfSZ19.rda")

#with CV
fitControl<-trainControl(method="repeatedcv",number=5,repeats=5)
dfCVRFfitSZ19<-train(V1~.,data=trainingdfSZ19,method = "rf",importance=TRUE,trControl=fitControl)
varImpPlot(dfCVRFfitSZ19$finalModel)
predictdfCVSZ19<-predict(dfCVRFfitSZ19,newdata = testingdfSZ19)
confusionMatrix(testingdfSZ19$V1,predictdfCVSZ19)

# Plot variable importance
png(filename = "dfRFfitSZ19.png", width = 1200, height = 800)
varImpPlot(dfRFfitSZ19$finalModel)
dev.off()

#New Era (randomForest package) ==========================================
#With randomForest library:
library(randomForest)
library(R.matlab)
set.seed(42)

SZ18orig<-readMat("18_subset.mat")                                         #read file
dfSZ18<-as.data.frame(SZ18orig)                                            #prepare
dfSZ18$targets<-as.factor(dfSZ18$targets)
dfSZ18.rf<-randomForest(targets~.,data=dfSZ18,strata=TRUE,importance=TRUE) #training model

#Comparing Votes
CompVotesSZ18<-cbind(dfSZ18$targets,dfSZ18.rf$votes[,1:7])                 #merge columns to have identificator
CompVotesSZ18<-as.data.frame(CompVotesSZ18)
AvCompVotesSZ18class1<-colMeans(CompVotesSZ18[CompVotesSZ18$V1==1,])       #comparing votes for 1st class
AvCompVotesSZ18class1<-AvCompVotesSZ18class1[-1]

#Comparing Votes (std)
AvCompVotesSZ18class1_std<-colSds(data.matrix(CompVotesSZ18[CompVotesSZ18$V1==1,]))       #comparing votes for 1st class
AvCompVotesSZ18class1_std<-AvCompVotesSZ18class1_std[-1]
AvCompVotesSZ18class1_std<-round(AvCompVotesSZ18class1_std, digits = 2)

#Standard deviation
library(matrixStats)
compvotesSZ18<-cbind(dfSZ18$targets,dfSZ18.rf$votes[,1:7])
compvotesSZ18<-as.data.frame(compvotesSZ18)
stdevclass1<-compvotesSZ18[compvotesSZ18$V1==1,]
stdevclass1<-stdevclass1[-1]
StDevVotesClass1<-rowSds(as.matrix.data.frame(stdevclass1))
StDevVotesClass1<-data.frame(index=c(1:length(StDevVotesClass1)),StDev=StDevVotesClass1)

#StDev ploting
AdjStDevClass1<-StDevVotesClass1$StDev
maxStDev<-max(StDevVotesClass1$StDev)
AdjStDevClass1<-AdjStDevClass1/maxStDev*100
png(filename = "Consistency Class 1.png", width = 800, height = 500)
plot(AdjStDevClass1, main="Consistency Class1", xlab="Index of a sample", ylab="Consistency in %")
dev.off()

#In .eps
setEPS()
postscript("Consistency, Class 1, AA.eps")
plot(AdjStDevClass1aa, main=NULL)
dev.off()


#Fail samples analysis for class 5
class1aa18fail<-subset(stdevclass1aa,stdevclass1aa[,1]<(1/7))

#Extracting names
library(seqinr)
mG<-read.fasta("MetabotropicGlutamate.fasta")
attributes(mG[[40]])$name

#Checking time
ptm <- proc.time() #Start the clock
dfSZ18.rf<-randomForest(targets~.,data=dfSZ18,strata=TRUE,importance=TRUE) #training model
proc.time() - ptm

#Loop for mG low section
for (i in c(25:110)) {
  +     print(attributes(mG[[i]])$name)
  + }

#Prob. garbage
failclass1AA14<-ifelse(rownames(class2aa14fail)==rownames(class2sz14fail)&rownames(class2dv14fail)==rownames(class2sz14fail)&rownames(class2aa14fail)==rownames(class2dv14fail),rownames(class2aa14fail),0)
failclass1AA14<-subset(compvotesAA14$V1==1,compvotesAA14$`1`>compvotesAA14$`2`,compvotesAA14$`1`>compvotesAA14$`3`,compvotesAA14$`1`>compvotesAA14$`4`,compvotesAA14$`1`>compvotesAA14$`5`,compvotesAA14$`1`>compvotesAA14$`6`,compvotesAA14$`1`>compvotesAA14$`7`)