library(stringr)
library(randomForest)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(survival)
library(survivalROC)
library(survminer)
library(pROC)
library(data.table)
a <- "LUAD"
a <- commandArgs(T)
dds <- fread('dds_cox.csv',head=T,sep=",")
dds <- fread('dds_clinical.csv')
dds2 <- dds[,-2:-10]


sample <- c()
for (i in 1:nrow(dds)){
  sample <- c(sample,str_sub(dds$ID_original[i],14L,15L)=='01')}

hsi <- fread("/Users/hexintong/Documents/WD/data/LUAD.csv",head =T)
setkey(hsi,samples)
setkey(dds2,ID)
dds3 <- hsi[dds2]
dds4 <- na.omit(dds3)

i <- 3
mean_i <- mean(dds[,i])
for (rows in 1:dim(dds)[2])
{if (dds[rows,i] > mean_i){
    dds[rows,i] <- 1}
  else{dds[rows,i] <- 0
    }}
sample <- as.factor(sample)
dds$sample <- sample
dds$sample  <- as.factor(sample)
ddt <- dds[,-c(2:10)]
sub <- sample(nrow(ddt),2/3*nrow(ddt))
train_data <- ddt[sub,]
test_data <- ddt[-sub,]
rf <- randomForest(as.factor(sample)  ~ .,
                   data = train_data,
                   ntree = 500,
                   mtry = 3,
                   importance = TRUE ,
                   proximity = TRUE)
#varImpPlot(rf, main = "variable importance")
rf_imp <- rf$importance
rf_imp2 <- as.data.frame(rf_imp)
rf_imp3 <- rfrf_imp <- rf$importance
rf_imp2 <- as.data.frame(rf_imp)
rf_imp3 <- rf_imp2[order(rf_imp2$MeanDecreaseAccuracy,decreasing =T),]
#rf_imp4 <- rf_imp3[1:200,]
rf_imp4 <- rf_imp3[1:200,]
#annotation
rf_imp4$Ensembl_ID <- rownames(rf_imp4)
Ensembl_ID <- rf_imp4$Ensembl_ID
# remove Ensembl version
Ensembl_ID <- str_replace(Ensembl_ID,pattern = ".[0-9]+$",
                          replacement = "")
#save file
fwrite(rf_imp4,'Rf.csv')
file <- paste('/Users/hexintong/Documents/WD/run/',a,'.csv',sep='')
write.table(rf_imp4,file)
#checkout
pre_ran <- predict(rf,newdata=test_data)
obs_p_ran = data.frame(prob=pre_ran,obs=test_data$sample)
table(test_data$sample,pre_ran,dnn=c(TRUE,FALSE))

#pred<-predict(rf,newdata=test_data)  
#pred_out_1<-predict(object=rf,newdata=test_data,type="prob")  #输出概率
#table <- table(pred,test_data$sample)  
#sum(diag(table))/sum(table)  #预测准确率
#plot(margin(rf,test_data$sample),main=观测值被判断正确的概率图)
#sort(rf,decreasing
dds2 <-dds[,1:10]
dds3 <-dds[colnames (dds) %in% rownames(rf_imp4)]
dds4 <-cbind(dds2,dds3)
dds4$status <-gsub('Dead',1,dds4$status)
dds4$status <-gsub('Alive',0,dds4$status)
dds4$status <- as.numeric(dds4$status)
cox <- coxph(Surv(dds4$OS,dds4$status)~.,data=dds4[,11:210])
cox.summary <- summary(cox)
# dds4 <- dds4
# cox <- coxph(Surv(dds4$OS,dds4$status)~.,data=dds4[,11:110])
# risk <-  predict(cox,type="risk",data=dds4[,11:110])
# calculate marker
marker0 <- cox$linear.predictor

# .C makes it a curve
SROC<-survivalROC.C(Stime=dds4$OS,
                  status=dds4$status,
                  marker = marker0,
                  predict.time =1825
)
pdf('survivalROC.pdf')
plot(SROC$FP,SROC$TP,type="l",
     col="red",
     xlim=c(0,1),ylim=c(0,1),
     xlab=paste("FP","\n","AUC=",round(SROC$AUC,3)),
     ylab="TP",
     main="5-year survival ROC")
abline(0,1,col="gray",lty=2)
legend("bottomright",c("GENES"),col="red",lty=c(1,1))
dev.off()
pdf('survivalplot.pdf')
# survival
n <- 45
dds5 <- dds4
dds5[,n] <- dds4[,n] > mean(dds4[,n])
dds5[,n] <- gsub('TRUE',"HIGH",dds5[,n])
dds5[,n] <- gsub('FALSE',"LOW",dds5[,n])
fit <- survfit(Surv(OS, status) ~ dds5[,n], data = dds5)
ggsurvplot(fit, pval = TRUE)
dev.off()

