#load clinical data
library(data.table)
library(stringr)
library(DESeq2)
library(survival)
cox <-fread('clinical.csv')
#calculate PFS & OS
cox[is.na(cox)]  <- 0
for ( i in 1:nrow(cox)){
	cox$OS[i] <- max(cox$days_to_death[i],cox$days_to_last_followup[i])}
for ( i in 1:nrow(cox)){
  if(cox$days_to_new_tumor[i]!=0){cox$PFS[i] <- cox$days_to_new_tumor[i]}else{cox$PFS[i] <- cox$days_to_last_followup[i]}
  }
load('dds.Rdata')
dds.norm <- counts(dds,normalized =T,replaced =F)
dds.df <-as.data.frame(dds.norm)
#read diff
res.diff <-read.table('res.diff.csv',head=T)
row <- rownames(res.diff)
dds.subset <-subset(dds.df,rownames(dds.df) %in% row)
dds.t <- t(dds.subset)
#write normalized expression data within clinical data
dds.dt <- as.data.table(dds.t,keep.rownames='ID_original')
dds.dt$ID <-str_sub(dds.dt$ID_original,1L,12L)
setkey(dds.dt,ID)
setkey(cox,ID)
dds.dt2  <- cox[dds.dt]
dds.all <- as.data.frame(dds.dt2)
dds <- subset(dds.all,str_sub(dds.all$ID_original,14L,15L)=='01')
dds$PFS <-as.numeric(dds$PFS)
dds$status <-gsub("Dead",1,dds$status)
dds$status <-gsub("Alive",0,dds$status)
dds$status <-as.numeric(dds$status)
write.table(dds,'dds_clinical.csv',sep=',')
cox <- coxph(Surv(dds$PFS,dds$status)~dds[,11],data=dds)
colnames.dds <-colnames(dds)
cox.coef <- data.frame(ID=0,coef=0,exp=0,Pr=0)
for (i in 11:(ncol(dds)-1)){
        cox.summary <- summary(coxph(Surv(dds$PFS,dds$status)~dds[,i],data=dds))
        cox.sm <- cox.summary$coefficients[,5]
        cox.coef[i,1] <- colnames.dds[i]
        cox.coef[i,2] <- cox.summary$coefficients[,1]
        cox.coef[i,3] <- cox.summary$coefficients[,2]
        cox.coef[i,4] <- cox.sm}
cox.minus <- cox.coef[-1,]
cox.pr <- subset(cox.minus,cox.minus[,4] < 0.01)
cox.na <- na.omit(cox.pr)
#num <- c(1,7,8)
num <- c(1:10)
for ( i in 11:ncol(dds)){
	if (colnames(dds)[i] %in% cox.na[,1]){
	num <- c(num,i)}}
dds2 <- dds.all[,num]
dds2$status <-gsub("Dead",1,dds2$status)
dds2$status <-gsub("Alive",0,dds2$status)
dds2$status <-as.numeric(dds2$status)
fwrite(dds2,'dds_cox.csv')
# # cox univiriant analysis
# library("survival")
# library("survminer")
# res.cox <- coxph(Surv(dds2$PFS,dds2$status)~ ENSG00000001084.13,data=dds2)
# res.cox <- survfit(Surv(dds2$PFS,dds2$status)~ ENSG00000001084.13,data=dds2)
# # Plot the baseline survival function
# ggsurvplot(res.cox, color = "#2E9FDF",
#            ggtheme = theme_minimal())
