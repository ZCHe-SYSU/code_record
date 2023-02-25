#go analysis
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(dplyr)
library(stringr)
library(data.table)

# Read differential expression data
res_diff <- fread('res_diff.csv') 
# Remove the version number of ensembl
res_diff$ID  <- res_diff$ID %>% str_replace(.,pattern = ".[0-9]+$",replacement = "")
Symbol_ID <- bitr(res_diff$ID, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")
# set key value for tables to merge
setkey(res_diff,"ID")
Symbol_ID <- Symbol_ID %>% as.data.table
setkey(Symbol_ID,"ENSEMBL")
# annotation
res_diff_annotated <-res_diff[Symbol_ID]

# Extracting information
res_cp <- data.frame(ID=res_diff_annotated$ID,SYMBOL=res_diff_annotated$SYMBOL,logFC=res_diff_annotated$log2FoldChange)
res_cp_logFC <- res_cp$logFC
names(res_cp_logFC) <- res_cp$ID
geneList <- sort(res_cp_logFC,decreasing = T)

# running go analysis
go_result <- gseGO(geneList     = geneList,
                   keyType= "ENSEMBL",
                   ont = "BP",
                   OrgDb = org.Hs.eg.db,
                   minGSSize    = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 1,
                   nPermSimple = 100000,
                   eps = 0)
# set readable
go <- setReadable(go_result, OrgDb = org.Hs.eg.db)

# store dot plot
pdf('diff_dotplot.pdf')
dotplot(go)
dev.off()


# store cnet plot
pdf('diff_cnetplot.pdf')
p1 <- cnetplot(go,showCategory=3,node_label="category")
p2 <- cnetplot(go,showCategory=3,node_label="gene") 
p3 <- cnetplot(go,showCategory=3,node_label="all") 
p4 <- cnetplot(go,showCategory=3,node_label="none") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
dev.off()

# gene_entrezid <- bitr(res_merged$gene_id, fromType = "ENSEMBL",
#                       toType = "ENTREZID",
#                       OrgDb = org.Hs.eg.db)
# res_cp_merge <- merge(gene.entrezid,res.cp,by.x='ENSEMBL',by.y='ID')
res_cp_logFC2 <- res_diff_annotated$log2FoldChange
names(res_cp_logFC2) <- res_diff_annotated$ENTREZID
geneList <- sort(res_cp_logFC2,decreasing = T)
kegg <- gseKEGG(geneList     = geneList,
                organism     = 'hsa',
                nPerm        = 1000,
                minGSSize    = 10,
                maxGSSize = 500,
                pvalueCutoff = 1,
                verbose      = FALSE)
pdf('diff_keggdotplot.pdf')
dotplot(kegg)
dev.off()

#
res_diff <- fread('/Users/hexintong/Documents/WD/run/LUAD/rf_imp_7.csv') 
res_diff$ID  <- res_diff$Ensembl_ID 
Symbol_ID <- bitr(res_diff$Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")
# set key value for tables to merge
setkey(res_diff,"ID")
Symbol_ID <- Symbol_ID %>% as.data.table
setkey(Symbol_ID,"ENSEMBL")
# annotation
res_diff_annotated <-res_diff[Symbol_ID]
# res.merged <-read.table('resmerged2.csv',sep=',',head =T)
res.cp <- data.frame(ID=res_diff_annotated$ID,SYMBOL=res_diff_annotated$SYMBOL,logFC=res_diff_annotated$MeanDecreaseAccuracy)
res.cp.logFC <- res.cp$logFC
names(res.cp.logFC) <- res.cp$ID
geneList <- sort(res.cp.logFC,decreasing = T)
go_result <- gseGO(geneList     = geneList,
                   keyType= "ENSEMBL",
                   ont = "BP",
                   OrgDb = org.Hs.eg.db,
                   minGSSize    = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 1)
go <- setReadable(go_result, OrgDb = org.Hs.eg.db)
pdf('rf_dotplot.pdf')
dotplot(go)
dev.off()
#pdf('rf_barplot.pdf')
#barplot(go, showCategory=15)
#dev.off()
#pdf('rf_enrichMap.pdf')
#enrichMap(go)
#dev.off()
pdf('rf_goGraph.pdf')
dev.off()
pdf('rf_cnetplot.pdf')
p1 <- cnetplot(go,showCategory=3,node_label="category")
p2 <- cnetplot(go,showCategory=3,node_label="gene")   
p3 <- cnetplot(go,showCategory=3,node_label="all")   
p4 <- cnetplot(go,showCategory=3,node_label="none")   
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
dev.off()

gene.entrezid <- bitr(res_diff$Ensembl_ID, fromType = "ENSEMBL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
res.cp.merge <- merge(gene.entrezid,res.cp,by.x='ENSEMBL',by.y='ID')
res.cp.logFC2 <- res.cp.merge$logFC
names(res.cp.logFC2) <- res.cp.merge$ENTREZID
geneList <- sort(res.cp.logFC2,decreasing = T)
kegg <- gseKEGG(geneList     = geneList,
                organism     = 'hsa',
                nPerm        = 1000,
                minGSSize    = 10,
                maxGSSize = 500,
                pvalueCutoff = 1,
                verbose      = FALSE)
pdf('rf_keggdotplot.pdf')
dotplot(kegg)
dev.off()

