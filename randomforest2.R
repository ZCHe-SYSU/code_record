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

# debug Code
setwd("/Users/hexintong/Documents/WD/run/LUAD")
a <- "LUAD"
i <- 3

# reading project ID
a <- commandArgs(T)
# reading annotated expression data
dds <- fread('dds_cox.csv',head=T,sep=",")
# deleting clinical data
dds2 <- dds[,-2:-10]

# reading image learing result data
hsi <- fread("/Users/hexintong/Documents/WD/data/LUAD.csv",head =T)

# merge image learing result data and expression data
setkey(hsi,samples)
setkey(dds2,ID)
dds3 <- hsi[dds2]
dds4 <- na.omit(dds3)
dds4 <- dds4[complete.cases(dds4),]
dds <- as.data.frame(dds4)

# create a df to store accuracy
columns = c("i","accuracy") 
df_accuracy = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(df_accuracy) = columns


for (i in 2:47){
  # dividing each clusters into high/low
  mean_i <- mean(dds[,i],na.rm=T)
  dds$sample <- dds[,i] > mean_i
  sample<- dds$sample
  # deleting other clusters data
  ddt <- dds[,-1:-47]
  # 2/3 are training data and 1/3 are testing data
  sub <- sample(nrow(ddt),2/3*nrow(ddt))
  train_data <- ddt[sub,]
  test_data <- ddt[-sub,]
  # runing a random forest
  rf <- randomForest(as.factor(train_data$sample)  ~ .,
                     data = train_data,
                     ntree = 500,
                     mtry = 3,
                     importance = TRUE ,
                     proximity = TRUE)
  # output predicted probability
  pred<-predict(rf,newdata=test_data)  
  pred_out_1<-predict(object=rf,newdata=test_data,type="prob")
  table <- table(pred,test_data$sample)  
  # prediction accuracy
  accuracy <-sum(diag(table))/sum(table)
  # store accuracy in df
  df_accuracy[nrow(df_accuracy) + 1,] = c(i,accuracy)
}

# order rf learning result by accuracy
df_accuracy_order <- df_accuracy[order(df_accuracy$accuracy,decreasing = T),]

# select 5 clusters with highest accuracy 
df_accuracy_order_5 <- df_accuracy_order[1:5,]

i = 39
for (i in df_accuracy_order_5$i){
  # dividing each clusters into high/low
  mean_i <- mean(dds[,i],na.rm=T)
  dds$sample <- dds[,i] > mean_i
  sample<- dds$sample
  # deleting other clusters data
  ddt <- dds[,-1:-47]
  # 2/3 are training data and 1/3 are testing data
  sub <- sample(nrow(ddt),2/3*nrow(ddt))
  train_data <- ddt[sub,]
  test_data <- ddt[-sub,]
  # run randomforest again
  rf <- randomForest(as.factor(train_data$sample)  ~ .,
                     data = train_data,
                     ntree = 500,
                     mtry = 3,
                     importance = TRUE ,
                     proximity = TRUE)
  # output importance
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
  rf_imp4$Ensembl_ID <- str_replace(Ensembl_ID,pattern = ".[0-9]+$",
                            replacement = "")
  #save file
  file <- paste('/Users/hexintong/Documents/WD/run/',a,'/rf_imp','_',i,'.csv',sep='')
  fwrite(rf_imp4,file)
  #checkout
  pre_ran <- predict(rf,newdata=test_data)
  obs_p_ran = data.frame(prob=pre_ran,obs=test_data$sample)
  table(test_data$sample,pre_ran,dnn=c(TRUE,FALSE))
}



