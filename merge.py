#merge RNA-seq and clinical
import pandas as pd
from gtfparse import read_gtf
df  = read_gtf("/home/zc/Homo_sapiens.GRCh38.103.chr.gtf")
df2 = df[["gene_id","gene_name","gene_biotype"]]
df3 = pd.read_csv('restomerge.csv', index_col = 0,sep=' ' )
df4 = pd.merge(df2,df3,left_on='gene_id',right_index = True)
df5 = df4.drop_duplicates(['gene_id'])
df5.to_csv('resmerged.csv')
df32 = pd.read_csv('restomerge2.csv',index_col =0,sep =' ')
df42 = pd.merge(df2,df32,left_on = "gene_id",right_index = True)
df52 = df42.drop_duplicates(['gene_id'])
df52.to_csv('resmerged2.csv')