#data preprocessing
library(DESeq2)
library(data.table)
#b <- commandArgs(T)
count <- fread("counts.csv")
#only normal and primary tumor sample
num <- c()
for(i in 1:length(colnames(count))){
        if(grepl("11A|01A",unlist(strsplit(colnames(count),"-")[i])[4])){
        num <- c(num,i)}}
count2 <- cbind(count[,'ID'],count[,..num])
fwrite(count2,"count2_dt.csv")
count3 <- count2
#clear 0 data
count3$sum <- rowSums(count2[,2:ncol(count2)])
count4 <- subset(count3,count3$sum>0)
count5 <- count4[,ncol(count4):=NULL]
fwrite(count5,"count3_dt.csv")
#set condition
a <-c()
condition <-c()
for ( i in 2:length(colnames(count5))) {
        if(grepl("01A",unlist(strsplit(colnames(count5),"-")[i])[4])){
        a <- c(a,"tumor")
        }else{
                a <- c(a,"normal")
        }
}
condition <- factor(a,levels=c("normal","tumor"))
count6 <-as.data.frame(count5)
rownames(count6) <-count6$ID
dds.count <- count6[,-1]
count7 <- as.data.frame(t(count6[,-1]))    #change table to frame
count8 <- cbind(condition,count7)
write.table(count8,"countfinal_df.csv")
# 
dds.before <- DESeqDataSetFromMatrix(dds.count,colData=DataFrame(condition),design= ~ condition)
dds <- DESeq(dds.before)
save(dds,file='dds.Rdata')
#get result
res <- results(dds)
#set axis
datax <- res$log2FoldChange
datay <- (-1*log10(res$padj))
#set significant point(two side)
sign_point1 = (datax >= 1) & (datay>= 1.30103)
sign_point2 = (datax <=-1) & (datay>= 1.30103)
#set color
col_point = rep("#BCBABE", length(datax))
col_point[sign_point1] = rgb(1,0,0,0.8)
col_point[sign_point2] = rgb(0,1,0,0.8)
#create pdf
pdf('diff_exp.pdf')
#plot
plot(x=res$log2FoldChange, y=(-1*log10(res$padj)),
     xlim=c(-8,8),ylim=c(0,24),
     col=col_point,pch=16
     )
#set abline
abline(h=-1*log10(0.05),lwd=3,lty=3,col="#4C5B61")
abline(v=log2(2) ,lwd=3,lty=3,col="#4C5B61")
abline(v=log2(1/2) ,lwd=3,lty=3,col="#4C5B61")
dev.off()

#calculate res
res_padj <- subset(res,padj<0.001)
res_diff <- subset(res_padj,log2FoldChange>1|log2FoldChange <(-1))
res_df <- as.data.frame(res_diff)
res_df$ID <- rownames(res_df)
#signif 
fwrite(res_df,file='res_diff.csv')
#clinical data
c0 <- fread('Clinical/vital_status',head=F,col.names=c("ID","status"))
c4 <- fread('Clinical/days_to_death',head=F,col.names=c("ID","days_to_death"))
c5 <- fread('Clinical/days_to_last_followup',head=F,col.names=c("ID","days_to_last_followup"))
c1 <- fread('Clinical/gender',head=F,col.names=c("ID","gender"))
c3 <- fread('Clinical/pathologic_stage',head=F,col.names=c("ID","stage"))
c6 <- fread('Clinical/new',head=F,col.names=c("ID","days_to_new_tumor"))
setkey(c0,"ID")
setkey(c1,"ID")
setkey(c3,"ID")
setkey(c4,"ID")
setkey(c5,"ID")
setkey(c6,"ID")
c0 <-c6[c0]
c0 <-c5[c0]
c0 <-c4[c0]
c0 <-c3[c0]
c0 <-c1[c0]
clinic.dp <-as.data.frame(c0)
num <- c()
for (i in 2:nrow(clinic.dp)){
	if (clinic.dp[i,'ID']==clinic.dp[i-1,'ID']){
		num  <- c(num,i)}}
clinic <- clinic.dp[-num,]
fwrite(clinic,"clinical.csv")

