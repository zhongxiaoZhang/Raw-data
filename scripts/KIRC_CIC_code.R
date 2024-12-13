
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)

my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","black"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+#
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+#
    theme_bw()+#
    theme(text = element_text(family = 'Times',size = 15),
          legend.title = element_blank(),#
          legend.position = leg.pos,
    )+
    ylab(ylab)+#
    xlab(xlab)+#
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +#
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)#
  return(p)
}
#00.data.pre##########
##TCGA#################
tcga.cli<-read.delim('origin_datas/TCGA/Merge_KIRC_clinical.txt',sep='\t',header = T)
colnames(tcga.cli)[1:20]
tcga.cli$A7_Grade
tcga.cli=data.frame(Samples=tcga.cli$A0_Samples,Age=tcga.cli$A17_Age,Gender=tcga.cli$A18_Sex,
                    T.stage=tcga.cli$A3_T,N.stage=tcga.cli$A4_N,M.stage=tcga.cli$A5_M,
                    Stage=tcga.cli$A6_Stage,Grade=tcga.cli$A7_Grade,
                    Status=tcga.cli$A2_Event,OS.time=tcga.cli$A1_OS)
tcga.cli$Samples=paste0(tcga.cli$Samples,'-01')
rownames(tcga.cli)=tcga.cli$Samples
head(tcga.cli)


table(tcga.cli$T.stage)
tcga.cli$T.stage=gsub('[abc]','',tcga.cli$T.stage)


table(tcga.cli$N.stage)
tcga.cli$N.stage[tcga.cli$N.stage=='NX']<-NA

table(tcga.cli$M.stage)
tcga.cli$M.stage[tcga.cli$M.stage=='MX']<-NA

table(tcga.cli$Stage)
tcga.cli$Stage[tcga.cli$Stage=='']<-NA
tcga.cli$Stage=gsub('Stage ','',tcga.cli$Stage)

table(tcga.cli$Grade)
tcga.cli$Grade[tcga.cli$Grade %in% c('GX','Not Available')]<-NA

tcga.cli$OS.time
tcga.cli=tcga.cli[tcga.cli$OS.time>30,]
tcga.cli$Status
tcga.cli$OS=ifelse(tcga.cli$Status=='Alive',0,1)
table(tcga.cli$OS)
dim(tcga.cli)


#
tcga.data<-read.delim('origin_datas/TCGA/KIRC_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga.data[1:4,1:4]
table(substr(colnames(tcga.data),14,15))


sample_T=colnames(tcga.data)[which(as.numeric(substr(colnames(tcga.data),14,15))==1)]#
sample_T=intersect(sample_T,tcga.cli$Samples)
sample_N=colnames(tcga.data)[which(as.numeric(substr(colnames(tcga.data),14,15))==11)]#
length(sample_N)
tcga_type=data.frame(Samples=c(sample_T,sample_N),type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$type)
# Normal  Tumor 
# 72    513

genecode=read.delim('GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]


range(tcga.data)
tcga.data=log2(tcga.data[intersect(rownames(tcga.data),mrna_genecode$SYMBOL),tcga_type$Samples]+1)
range(tcga.data)
tcga.exp=tcga.data[,sample_T]
dim(tcga.exp)
# 19503   513
tcga.cli=tcga.cli[sample_T,]
dim(tcga.cli)
head(tcga.cli)

##ICGC########
icgc.cli<-read.delim('origin_datas/ICGC/Merge_clinical.txt',sep='\t',header = T)
colnames(icgc.cli)[1:30]
table(icgc.cli$specimen_type)
icgc.cli=icgc.cli[which(icgc.cli$specimen_type=='Primary tumour - solid tissue'),]
icgc.cli$donor_tumour_stage_at_diagnosis
icgc.cli$tumour_stage
icgc.cli=data.frame(Samples=icgc.cli$icgc_sample_id,
                    Age=icgc.cli$donor_age_at_diagnosis,
                    Gender=icgc.cli$donor_sex,
                    Stage=icgc.cli$donor_tumour_stage_at_diagnosis,
                    OS.time=icgc.cli$donor_survival_time,
                    Status=icgc.cli$donor_vital_status)
rownames(icgc.cli)=icgc.cli$Samples
head(icgc.cli)
icgc.cli$OS.time
table(icgc.cli$Status)
icgc.cli$OS=ifelse(icgc.cli$Status=='alive',0,1)

icgc.cli$T.stage=substr(icgc.cli$Stage,1,2)
icgc.cli$N.stage=substr(icgc.cli$Stage,3,4)
icgc.cli$M.stage=substr(icgc.cli$Stage,5,6)
icgc.cli$Stage=rep(' ',nrow(icgc.cli))
icgc.cli$Stage[which(icgc.cli$M.stage=='MX'|icgc.cli$N.stage=='NX')]='Unknow'
icgc.cli$Stage[which(icgc.cli$M.stage=='M1'|icgc.cli$T.stage=='T4')]='IV'
icgc.cli$Stage[which(icgc.cli$T.stage=='T3' & icgc.cli$M.stage=='M0')]='III'
icgc.cli$Stage[which(icgc.cli$T.stage=='T1' & icgc.cli$N.stage=='N1')]='III'
icgc.cli$Stage[which(icgc.cli$T.stage=='T2' & icgc.cli$N.stage=='N1')]='III'
icgc.cli$Stage[which(icgc.cli$Stage==' ' & icgc.cli$T.stage=='T1')]='I'
icgc.cli$Stage[which(icgc.cli$Stage==' ' & icgc.cli$T.stage=='T2')]='II'
head(icgc.cli)



icgc.data=read.delim('origin_datas/ICGC/ICGC_renal_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
icgc.data[1:4,1:4]
dim(icgc.data)
range(icgc.data)
icgc.data=log2(icgc.data+1)
icgc.exp=icgc.data[,intersect(colnames(icgc.data),rownames(icgc.cli))]
dim(icgc.exp)


icgc_type=data.frame(Samples=c(colnames(icgc.exp),setdiff(colnames(icgc.data),rownames(icgc.cli))),
                     type=rep(c('Tumor','Normal'),c(91,45)))
rownames(icgc_type)=icgc_type$Samples
table(icgc_type$type)






#01.########
dir.create('results/01.DEGs')
tcga.limma=mg_limma_DEG(exp = tcga.data[,tcga_type$Samples],group = tcga_type$type,ulab = 'Tumor',dlab = 'Normal')
tcga.limma$Summary
tcga.degs=tcga.limma$DEG[tcga.limma$DEG$adj.P.Val<0.05 & abs(tcga.limma$DEG$logFC)>1,]
write.csv(tcga.degs,'results/01.DEGs/TCGA_DEGs.csv')
tcga.up.degs=rownames(tcga.degs[tcga.degs$logFC>0,])
tcga.dn.degs=rownames(tcga.degs[tcga.degs$logFC<0,])
fig1a=my_volcano(dat = tcga.limma,p_cutoff = 0.05,fc_cutoff = 1,
                 col = c('#E6AB02','#8DD3C7','grey'),leg.pos = 'top')+ggtitle('TCGA')
fig1a

icgc.limma=mg_limma_DEG(exp = icgc.data[,icgc_type$Samples],group = icgc_type$type,ulab = 'Tumor',dlab = 'Normal')
icgc.limma$Summary
icgc.degs=icgc.limma$DEG[icgc.limma$DEG$adj.P.Val<0.05 & abs(icgc.limma$DEG$logFC)>1,]
write.csv(icgc.degs,'results/01.DEGs/ICGC_DEGs.csv')
icgc.up.degs=rownames(icgc.degs[icgc.degs$logFC>0,])
icgc.dn.degs=rownames(icgc.degs[icgc.degs$logFC<0,])
fig1b=my_volcano(dat = icgc.limma,p_cutoff = 0.05,fc_cutoff = 1,
                 col = c('#E6AB02','#8DD3C7','grey'),leg.pos = 'top')+ggtitle('ICGC')
fig1b

up.degs=intersect(tcga.up.degs,icgc.up.degs)
dn.degs=intersect(tcga.dn.degs,icgc.dn.degs)
length(up.degs);length(dn.degs)
all.dgs=c(up.degs,dn.degs)
length(all.dgs)
#1968


CIC.genes=read.xlsx('origin_datas/CIC_related_genes_PMID35974300.xlsx')
CIC.genes=CIC.genes$genes
length(CIC.genes)
#101



library(eulerr)
v=list(tcga.up.degs,icgc.up.degs)
names(v)=c('TCGA up DEGs','ICGC up DEGs')
fig1c=plot(venn(v),labels = list(col = "gray20", font = 2), 
               edges = list(col="gray60", lex=1),
               fills = list(fill = c("#FB8072","#80B1D3" ,"#B3DE69"), alpha = 0.6),
               quantities = list(cex=1, col='gray20'))
fig1c


DE.CIC.genes=intersect(CIC.genes,up.degs)
DE.CIC.genes
length(DE.CIC.genes)
#12


DEGs.enrichKEGG = mg_clusterProfiler(genes = up.degs)
DEGs.enrichKEGG.res=DEGs.enrichKEGG$KEGG@result
table(DEGs.enrichKEGG.res$p.adjust<0.05)
DEGs.enrichKEGG.res=DEGs.enrichKEGG.res[DEGs.enrichKEGG.res$p.adjust<0.05,]
write.xlsx(DEGs.enrichKEGG.res,'results/01.DEGs/upDEGs.enrichKEGG.res.xlsx',overwrite = T)
head(DEGs.enrichKEGG.res)
DEGs.enrichKEGG.res=DEGs.enrichKEGG.res %>% slice_max(n =10, order_by = Count)
DEGs.enrichKEGG.res$Description <- factor(DEGs.enrichKEGG.res$Description, 
                                                levels = DEGs.enrichKEGG.res$Description[order(DEGs.enrichKEGG.res$Count,decreasing = F)])
# 
fig1d=ggplot(data = DEGs.enrichKEGG.res, aes(x = Count, y = Description, fill=-log10(p.adjust))) +
  geom_bar(width = 0.5,stat = 'identity') +
  theme_classic() +
  scale_x_continuous(expand = c(0,0.5)) +
  scale_fill_gradient(low = "#80B1D3", high = "#E6AB02")+
  theme(axis.text.y = element_blank()) +
  geom_text(data = DEGs.enrichKEGG.res,aes(x = .5, y = Description, label = Description),
            family = 'Times',size = 5,hjust = 0) +
  theme(text = element_text(family = 'Times',size = 15,face = 'bold'),
        legend.position = 'right')+ggtitle('TOP10 KEGG enrichment analysis')
fig1d

# 
degs.genrichGO.res=rbind(DEGs.enrichKEGG$GO_BP@result,DEGs.enrichKEGG$GO_CC@result,DEGs.enrichKEGG$GO_MF@result)
degs.genrichGO.res$ONTOLOGY=rep(c('BP','CC','MF'),c(nrow(DEGs.enrichKEGG$GO_BP@result),nrow(DEGs.enrichKEGG$GO_CC@result),nrow(DEGs.enrichKEGG$GO_MF@result)))
head(degs.genrichGO.res)

top5 <- degs.genrichGO.res %>%  group_by(ONTOLOGY) %>%  arrange(pvalue) %>%  slice_head(n = 5)
df_top5 <- rbind(subset(top5, ONTOLOGY=="BP"), subset(top5, ONTOLOGY=="CC"), subset(top5, ONTOLOGY=="MF"))
df_top5$ONTOLOGY <- factor(df_top5$ONTOLOGY, levels=c('BP','CC', 'MF'))
# df_top5$Description <- factor(df_top5$Description, levels = rev(df_top5$Description))
df_top5$Description <- factor(df_top5$Description, 
                              levels = df_top5$Description[order(df_top5$ONTOLOGY,decreasing = T)])

fig1e=ggplot(df_top5, aes(x = Description, y = Count, fill = ONTOLOGY, group = Description)) +
  geom_bar(stat="identity", position="dodge", colour=aes(ONTOLOGY))+
  scale_fill_manual(values =c("#297CA0", "#E9EA77",'pink'))+
  xlab('')+ggtitle('TOP5 GO enrichment analysis')+
  theme_classic()+theme(text = element_text(family = 'Times',size = 14,face = 'bold'),
                        axis.text.x = element_text(size = 14,angle = 45,hjust = 1),
                        legend.position = 'left')

fig1e


pdf('results/01.DEGs/Fig1.pdf',height = 17,width = 12)
mg_merge_plot(mg_merge_plot(fig1a,fig1b,labels = c('A','B')),
              mg_merge_plot(fig1c,fig1d,labels = c('C','D')),
              fig1e,nrow=3,labels = c('','','E'),heights = c(.8,.8,1))
dev.off()

#02.TCGA#########
dir.create('results/02.ML')
##RFE#####
library(caret)
tcga.ML.dat=t(tcga.data[DE.CIC.genes,])
tcga.ML.dat=cbind.data.frame(type=tcga_type$type,tcga.ML.dat[tcga_type$Samples,])
tcga.ML.dat$type=as.factor(tcga.ML.dat$type)

mat=as.matrix(tcga.ML.dat[,DE.CIC.genes])
group=as.factor(tcga.ML.dat$type)


set.seed(123)
rfeControl = rfeControl(functions = caretFuncs, method ="cv", number= 10, verbose = FALSE)
tcga.rf1 = rfe(x = mat,y = group,sizes=c(1:20),rfeControl = rfeControl, method ="svmLinear")
result_svm = tcga.rf1$result


# summarize the results
print(tcga.rf1)
# list the chosen features
predictors(tcga.rf1)#2
#[1] "CDKN2A" "VIM"    "GZMB"   "TGFB1"  "CTSS"   "CDC20" 

pdf('results/02.ML/rfe_plot.pdf',height = 5,width = 6)
plot(tcga.rf1, type=c("g", "o"), main = 'TCGA RFE',xlab='Gene Number')
dev.off()

##LASSO##########
#lasso
library(glmnet)
set.seed(321)
fit1=glmnet(x = mat,y = group,family = "binomial",nlambda=100, alpha=1) 
cv.fit<-cv.glmnet(x = mat,y = group,family = "binomial",nlambda=100, alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
cv.fit$lambda.min
#0.001517271

pdf('results/02.ML/LASSO.pdf',height = 5,width = 9,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()

coefficient = coef(cv.fit,s=cv.fit$lambda.min)
Active.Index = which(as.numeric(coefficient)!=0)
active.coefficients = as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox = rownames(coefficient)[Active.Index]
sig_gene_multi_cox = sig_gene_multi_cox[-1]
active.coefficients=active.coefficients[-1]
data.frame(gene=sig_gene_multi_cox,coef=active.coefficients)
write.table(data.frame(gene=sig_gene_multi_cox,coef=active.coefficients),
'results/02.ML/LASSO_result.txt',
sep = '\t', col.names = NA, quote = FALSE)



hub.genes=intersect(tcga.rf1$optVariables,sig_gene_multi_cox)
hub.genes

v=list(tcga.rf1$optVariables,sig_gene_multi_cox)
names(v)=c('RFE','LASSO')
fig2c=plot(venn(v),labels = list(col = "gray20", font = 2), 
               edges = list(col="gray60", lex=1),
               fills = list(fill = c("#297CA0", "#E9EA77"), alpha = 0.6),
               quantities = list(cex=.8, col='gray20'))
fig2c


tcga.expr.df=cbind.data.frame(t(tcga.data[hub.genes,tcga_type$Samples]),type=tcga_type$type)
tcga.expr.df=melt(tcga.expr.df)
head(tcga.expr.df)

fig2d=ggviolin(tcga.expr.df, x = "variable", y = "value", fill = "type",add = "boxplot")+
  ggpubr::stat_compare_means(aes(group=type), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values = pal_simpsons()(9)[c(2,5)])+ylab('Gene Expression')+xlab('')+
  theme(text = element_text(family = 'Times',size=15))+ggtitle('TCGA')
fig2d

icgc.expr.df=cbind.data.frame(t(icgc.data[hub.genes,icgc_type$Samples]),type=icgc_type$type)
icgc.expr.df=melt(icgc.expr.df)
head(icgc.expr.df)


fig2e=ggviolin(icgc.expr.df, x = "variable", y = "value", fill = "type",add = "boxplot",width = 1.2)+
  ggpubr::stat_compare_means(aes(group=type), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values = pal_simpsons()(9)[c(2,5)])+ylab('Gene Expression')+xlab('')+
  theme(text = element_text(family = 'Times',size=15))+ggtitle('ICGC')
fig2e

pdf('results/02.ML/Fig2CDE.pdf',height = 5,width = 15)
mg_merge_plot(fig2c,fig2d,fig2e,ncol=3,labels = LETTERS[3:5])
dev.off()

#03.ROC曲线###########
dir.create('results/03.Diagnosis.model')
library(pROC)
TCGA.roc <- list()
for (i in 1:length(hub.genes)){
  roc1 <- roc(tcga.ML.dat$type, tcga.ML.dat[,hub.genes[i]])
  TCGA.roc[[i]]=roc1
  names(TCGA.roc)[i] <- paste0(hub.genes[i],' AUC=',round(roc1$auc[1],2))
}

fig3a=ggroc(TCGA.roc)+
  geom_segment(aes(x = 1, y = 0, xend =0, yend = 1), color="darkgrey", linetype="dashed")+
  ggtitle(label = 'TCGA')+
  scale_color_manual(values =pal_tron()(7)[c(1,2,3,4,7)])+
  theme_bw()+theme(panel.grid = element_blank(),legend.position =c(0.75,0.3),
                   text = element_text(family = 'Times',size = 14))
fig3a
ggsave('results/03.Diagnosis.model/Fig3a.pdf',fig3a,height = 5,width = 5)



icgc.ML.dat=t(icgc.data[DE.CIC.genes,])
icgc.ML.dat=cbind.data.frame(type=icgc_type$type,icgc.ML.dat[icgc_type$Samples,])
icgc.ML.dat$type=as.factor(icgc.ML.dat$type)

library(pROC)
icgc.roc <- list()
for (i in 1:length(hub.genes)){
  roc1 <- roc(icgc.ML.dat$type, icgc.ML.dat[,hub.genes[i]])
  icgc.roc[[i]]=roc1
  names(icgc.roc)[i] <- paste0(hub.genes[i],' AUC=',round(roc1$auc[1],2))
}

fig3b=ggroc(icgc.roc)+
  geom_segment(aes(x = 1, y = 0, xend =0, yend = 1), color="darkgrey", linetype="dashed")+
  ggtitle(label = 'ICGC')+
  scale_color_manual(values =pal_tron()(7)[c(1,2,3,4,7)])+
  theme_bw()+theme(panel.grid = element_blank(),legend.position =c(0.75,0.3),
                   text = element_text(family = 'Times',size = 14))
fig3b
ggsave('results/03.Diagnosis.model/Fig3b.pdf',fig3b,height = 5,width = 5)





#####
library(e1071)
paste0(hub.genes,collapse = '+')
tcga.svm.model=svm(type ~CDKN2A+VIM+TGFB1+CTSS+CDC20,data = tcga.ML.dat)
summary(tcga.svm.model)
tcga.svm.model<-predict(tcga.svm.model)
tcga.svm.model <- as.ordered(tcga.svm.model)

modelroc_TCGA <- roc(tcga.ML.dat$type,tcga.svm.model)
modelroc_TCGA
pdf('results/03.Diagnosis.model/TCGA_Diagnosis_AUC.pdf',height = 5,width = 5)
plot(modelroc_TCGA, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.5, 0.2), grid.col=c("green", "red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()

train_table<-table(tcga.svm.model,tcga.ML.dat$type)
#classAgreement
classAgreement(train_table)
#调用confusionMatrix
train_result_matrix=confusionMatrix(train_table)
draw_confusion_matrix <- function(cm,disease='Case',control='Control',group.col=c('#3F97D0','#F7AD50')) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col=group.col[1])
  text(195, 435, control, cex=1.2)
  rect(250, 430, 340, 370, col=group.col[2])
  text(295, 435, disease, cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=group.col[2])
  rect(250, 305, 340, 365, col=group.col[1])
  text(140, 400, control, cex=1.2, srt=90)
  text(140, 335, disease, cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  

pdf('results/03.Diagnosis.model/TCGA_confusion_mat.pdf',height = 5,width = 5)
draw_confusion_matrix(cm = train_result_matrix,disease = 'Tumor',group.col=pal_simpsons()(9)[c(2,5)])
dev.off()







library(e1071)
paste0(hub.genes,collapse = '+')
icgc.svm.model=svm(type ~CDKN2A+VIM+TGFB1+CTSS+CDC20,data = icgc.ML.dat)
summary(icgc.svm.model)
icgc.svm.model<-predict(icgc.svm.model)
icgc.svm.model <- as.ordered(icgc.svm.model)

modelroc_icgc <- roc(icgc.ML.dat$type,icgc.svm.model)
modelroc_icgc
pdf('results/03.Diagnosis.model/ICGC_Diagnosis_AUC.pdf',height = 5,width = 5)
plot(modelroc_icgc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.5, 0.2), grid.col=c("green", "red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()

test_table<-table(icgc.svm.model,icgc.ML.dat$type)
#
classAgreement(test_table)
#
test_result_matrix=confusionMatrix(test_table)
pdf('results/03.Diagnosis.model/ICGC_confusion_mat.pdf',height = 5,width = 5)
draw_confusion_matrix(cm = test_result_matrix,disease = 'Tumor',group.col=pal_simpsons()(9)[c(2,5)])
dev.off()
 


#04.##########
dir.create('results/04.TME')
tcga.cibersort=read.delim('results/04.TME/TCGA_CIBERSORT_Results.txt',row.names = 1,check.names = F)
tcga.cibersort=tcga.cibersort[,1:22]


tcga.cibersort.df=data.frame(type=tcga_type$type,tcga.cibersort[tcga_type$Samples,])
tcga.cibersort.df=melt(tcga.cibersort.df)
head(tcga.cibersort.df)

tcga.cibersort.df %>%
  ggplot(aes(x=variable, y=value,fill=type)) +
  geom_boxplot()+
  scale_fill_manual(values = pal_simpsons()(9)[c(2,5)])+   #
  ggpubr::stat_compare_means(aes(group=type), label = "p.signif", method = 'wilcox.test')+
  labs(x="", y = "Fraction", fill = "type") +
  theme_classic()+
  theme(legend.position = 'top',text = element_text(family = 'Times',size = 14,color='black',face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) # 
ggsave('results/04.TME/TCGA_CIBERSORT_boxplot.pdf',height = 5,width = 12)


tcga.ciber.cor.df=cbind.data.frame(t(tcga.exp[hub.genes,tcga.cli$Samples]),tcga.cibersort[tcga.cli$Samples,])
cor_res <- Hmisc::rcorr(as.matrix(tcga.ciber.cor.df),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 1
cor_res$r[is.na(cor_res$r)] <- 0
library(corrplot)
pdf('results/04.TME/TCGA_cibersort_cor_plot.pdf',height = 12,width = 6)
corrplot(as.matrix(cor_res$r[colnames(tcga.cibersort),hub.genes]),
         p.mat = as.matrix(cor_res$P[colnames(tcga.cibersort),hub.genes]),
         mar = c(0,0,1,1),
         col=colorRampPalette(c('#8787FFFF', 'white','#FB8072'))(100),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[1],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.1,cl.cex = 0.8,
         addgrid.col = 'white',
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[6],
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()







icgc.cibersort=read.delim('results/04.TME/ICGC_CIBERSORT_Results.txt',row.names = 1,check.names = F)
icgc.cibersort=icgc.cibersort[,1:22]

icgc.cibersort.df=data.frame(type=icgc_type$type,icgc.cibersort[icgc_type$Samples,])
icgc.cibersort.df=melt(icgc.cibersort.df)
head(icgc.cibersort.df)

icgc.cibersort.df %>%
  ggplot(aes(x=variable, y=value,fill=type)) +
  geom_boxplot()+
  scale_fill_manual(values = pal_simpsons()(9)[c(2,5)])+   #
  ggpubr::stat_compare_means(aes(group=type), label = "p.signif", method = 'wilcox.test')+
  labs(x="", y = "Score", fill = "type") +
  theme_classic()+
  theme(legend.position = 'top',text = element_text(family = 'Times',size = 14,color='black',face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) # 
ggsave('results/04.TME/ICGC_CIBERSORT_boxplot.pdf',height = 5,width = 12)

icgc.ciber.cor.df=cbind.data.frame(t(icgc.exp[hub.genes,icgc.cli$Samples]),icgc.cibersort[icgc.cli$Samples,])
cor_res <- Hmisc::rcorr(as.matrix(icgc.ciber.cor.df),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 1
cor_res$r[is.na(cor_res$r)] <- 0
library(corrplot)
pdf('results/04.TME/ICGC_cibersort_cor_plot.pdf',height = 12,width = 6)
corrplot(as.matrix(cor_res$r[colnames(icgc.cibersort),hub.genes]),
         p.mat = as.matrix(cor_res$P[colnames(icgc.cibersort),hub.genes]),
         mar = c(0,0,1,1),
         col=colorRampPalette(c('#8787FFFF', 'white','#FB8072'))(100),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[1],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.1,cl.cex = 0.8,
         addgrid.col = 'white',
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[6],
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()




#05.###########
dir.create('results/05.Clinical_correlation')
##TCGA###########
tcga.cli.merge=data.frame(tcga.cli,t(tcga.exp[hub.genes,tcga.cli$Samples]))
head(tcga.cli.merge)

stage.p=list()
##########Stage
stage.p[[1]]=tcga.cli.merge%>% drop_na(Stage)%>%
  ggplot( aes(x=Stage, y=CDKN2A,fill=Stage)) +geom_boxplot()+
  scale_fill_manual(values = pal_futurama()(9))+   #
  ggpubr::stat_compare_means(aes(group=Stage), label = "p.format", method = 'kruskal.test')+
  labs(x="",y = "gene expression", subtitle = "CDKN2A", fill = "Stage") +
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0.5),
                   text = element_text(family = 'Times',size = 14,color='black',face = 'bold')) # 

stage.p[[2]]=tcga.cli.merge%>% drop_na(Stage)%>%
  ggplot( aes(x=Stage, y=VIM,fill=Stage)) +geom_boxplot()+
  scale_fill_manual(values = pal_futurama()(9))+   #
  ggpubr::stat_compare_means(aes(group=Stage), label = "p.format", method = 'kruskal.test')+
  labs(x="",y = "gene expression", subtitle = "VIM", fill = "Stage") +
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0.5),
                   text = element_text(family = 'Times',size = 14,color='black',face = 'bold')) # 

stage.p[[3]]=tcga.cli.merge%>% drop_na(Stage)%>%
  ggplot( aes(x=Stage, y=TGFB1,fill=Stage)) +geom_boxplot()+
  scale_fill_manual(values = pal_futurama()(9))+   #
  ggpubr::stat_compare_means(aes(group=Stage), label = "p.format", method = 'kruskal.test')+
  labs(x="",y = "gene expression", subtitle = "TGFB1", fill = "Stage") +
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0.5),
                   text = element_text(family = 'Times',size = 14,color='black',face = 'bold')) # 

stage.p[[4]]=tcga.cli.merge%>% drop_na(Stage)%>%
  ggplot( aes(x=Stage, y=CTSS,fill=Stage)) +geom_boxplot()+
  scale_fill_manual(values = pal_futurama()(9))+   #
  ggpubr::stat_compare_means(aes(group=Stage), label = "p.format", method = 'kruskal.test')+
  labs(x="", y = "gene expression",subtitle = "CTSS", fill = "Stage") +
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0.5),
                   text = element_text(family = 'Times',size = 14,color='black',face = 'bold')) # 


stage.p[[5]]=tcga.cli.merge%>% drop_na(Stage)%>%
  ggplot( aes(x=Stage, y=CDC20,fill=Stage)) +geom_boxplot()+
  scale_fill_manual(values = pal_futurama()(9))+   #
  ggpubr::stat_compare_means(aes(group=Stage), label = "p.format", method = 'kruskal.test')+
  labs(x="",y = "gene expression", subtitle = "CDC20", fill = "Stage") +
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0.5),
                   text = element_text(family = 'Times',size = 14,color='black',face = 'bold')) # 


stage.p.merge=mg_merge_plot(stage.p,ncol=5,labels = LETTERS[1:5],common.legend = T)
stage.p.merge

grade.p=list()
########Grade
grade.p[[1]]=tcga.cli.merge%>% drop_na(Grade)%>%
  ggplot( aes(x=Grade, y=CDKN2A,fill=Grade)) +geom_boxplot()+
  scale_fill_manual(values = pal_futurama()(9))+   #
  ggpubr::stat_compare_means(aes(group=Grade), label = "p.format", method = 'kruskal.test')+
  labs(x="", y = "gene expression",subtitle = "CDKN2A", fill = "Grade") +
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0.5),
                   text = element_text(family = 'Times',size = 14,color='black',face = 'bold')) # 

grade.p[[2]]=tcga.cli.merge%>% drop_na(Grade)%>%
  ggplot( aes(x=Grade, y=VIM,fill=Grade)) +geom_boxplot()+
  scale_fill_manual(values = pal_futurama()(9))+   #
  ggpubr::stat_compare_means(aes(group=Grade), label = "p.format", method = 'kruskal.test')+
  labs(x="",y = "gene expression", subtitle= "VIM", fill = "Grade") +
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0.5),
                   text = element_text(family = 'Times',size = 14,color='black',face = 'bold')) # 

grade.p[[3]]=tcga.cli.merge%>% drop_na(Grade)%>%
  ggplot( aes(x=Grade, y=TGFB1,fill=Grade)) +geom_boxplot()+
  scale_fill_manual(values = pal_futurama()(9))+   #
  ggpubr::stat_compare_means(aes(group=Grade), label = "p.format", method = 'kruskal.test')+
  labs(x="",y = "gene expression", subtitle = "TGFB1", fill = "Grade") +
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0.5),
                   text = element_text(family = 'Times',size = 14,color='black',face = 'bold')) # 


grade.p[[4]]=tcga.cli.merge%>% drop_na(Grade)%>%
  ggplot( aes(x=Grade, y=CTSS,fill=Grade)) +geom_boxplot()+
  scale_fill_manual(values = pal_futurama()(9))+   #
  ggpubr::stat_compare_means(aes(group=Grade), label = "p.format", method = 'kruskal.test')+
  labs(x="", y = "gene expression", fill = "Grade",subtitle = 'CTSS') +
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0.5),
                   text = element_text(family = 'Times',size = 14,color='black',face = 'bold')) # 


grade.p[[5]]=tcga.cli.merge%>% drop_na(Grade)%>%
  ggplot( aes(x=Grade, y=CDC20,fill=Grade)) +geom_boxplot()+
  scale_fill_manual(values = pal_futurama()(9))+   #
  ggpubr::stat_compare_means(aes(group=Grade), label = "p.format", method = 'kruskal.test')+
  labs(x="",y = "gene expression", subtitle = "CDC20", fill = "Grade") +
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0.5),
                   text = element_text(family = 'Times',size = 14,color='black',face = 'bold')) # 


grade.p.merge=mg_merge_plot(grade.p,ncol=5,common.legend = T)
grade.p.merge

pdf('results/05.Clinical_correlation/Fig5.pdf',height = 8,width = 16)
mg_merge_plot(stage.p.merge,grade.p.merge,nrow=2)
dev.off()

cox_batch(tcga.exp[hub.genes,tcga.cli$Samples],time = tcga.cli$OS.time,event = tcga.cli$OS)


head(tcga.cli.merge)

hubgene.km=list()
for (i in 1:5) {
  # cutoff<-surv_cutpoint(tcga.cli.merge,time="OS.time",event="OS",variables=hub.genes[i])
  # tcga.cli.merge$group=ifelse(tcga.cli.merge[,hub.genes[i]]>cutoff$cutpoint$cutpoint,'H','L')
  tcga.cli.merge$group=ifelse(tcga.cli.merge[,hub.genes[i]]>median(tcga.cli.merge[,hub.genes[i]]),'H','L')
  hubgene.km[[i]]=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ group,data = tcga.cli.merge),data=tcga.cli.merge,
             conf.int = F,pval = T,fun = "pct",risk.table =T, size = 1,surv.median.line = 'hv',
             linetype = 'solid', legend = 'top', legend.title =hub.genes[i])
  hubgene.km[[i]]=mg_merge_plot(hubgene.km[[i]]$plot,hubgene.km[[i]]$table,nrow=2,heights = c(2.5,1))
}

hubgene.km.merge=mg_merge_plot(hubgene.km,ncol=5)
hubgene.km.merge



#06.#######
dir.create('results/06.scRNA')
library(Seurat)
dir= 'origin_datas/GEO/GSE224630_RAW/'
dir_name=list.files(dir)
datalist=list()
for (i in 1:length(dir_name)){
  # i=1
  dir.10x = paste0(dir,dir_name[i])
  list.files(dir.10x)
  my.data <- Read10X(data.dir = dir.10x)
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], min.cells = 3, min.features = 200)
  datalist[[i]]$Samples=dir_name[i]
  # datalist[[i]]$tissue=str_split_fixed(dir_name[i],'_',3)[3]
  rm(my.data)
}
names(datalist)=dir_name

##########
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 
  datalist[[i]] <- sce
  rm(sce)
}
sce <- merge(x=datalist[[1]], y=datalist[2:length(datalist)])
#

raw_meta=sce@meta.data
raw_count <- table(raw_meta$Samples)
raw_count
sum(raw_count)#  14981
pearplot_befor<-VlnPlot(sce,group.by ='Samples',
                        features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        pt.size = 0,
                        ncol = 3)
pearplot_befor
ggsave('results/06.scRNA/pearplot_befor.pdf',pearplot_befor,height = 6,width = 15)
#
sce=subset(sce, subset=nFeature_RNA>200 & nFeature_RNA<5000 & percent.mt<10)

clean_meta=sce@meta.data
clean_count <- table(clean_meta$Samples)
clean_count
sum(clean_count)#11675
pearplot_after <- VlnPlot(sce,group.by ='Samples',
                          features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                          pt.size = 0,
                          ncol = 3)
pearplot_after
ggsave('results/06.scRNA/pearplot_after.pdf',pearplot_after,height = 6,width = 15)
rm(datalist)
########
# sce <- NormalizeData(sce)
# sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
# sce <- ScaleData(sce, features = rownames(sce))

sce = SCTransform(sce, vars.to.regress = "percent.mt", verbose = F)
#
sce <- RunPCA(sce, features = VariableFeatures(sce))
colnames(sce@meta.data)
##
library(harmony)
sce = RunHarmony(sce, group.by.vars="Samples", max.iter.harmony=50, lambda=0.5,assay.use = "SCT")
#
pca.plot=ElbowPlot(sce,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
pca.plot
ggsave('results/06.scRNA/pca.plot.pdf',pca.plot,height = 6,width = 6)

###
sce <- RunTSNE(sce, dims=1:25, reduction="harmony")
DimPlot(sce,group.by='Samples',reduction="tsne",label =F,pt.size = 0.2)
# sce <- RunTSNE(sce, dims=1:15, reduction="harmony")
library(clustree)
sce <- FindNeighbors(sce, dims = 1:25, reduction="harmony")
#
sce <- FindClusters(object = sce,resolution = .1)
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
length(table(sce@meta.data$seurat_clusters))
DimPlot(sce,group.by='seurat_clusters',reduction="tsne",label =T,pt.size = 0.2)
saveRDS(sce,file = 'results/06.scRNA/sce.rds')



###################
#B/Plasma cell  1
VlnPlot(sce,features = c('CD79A','MZB1','CD27','MS4A1'),pt.size = 0,group.by = 'seurat_clusters')
#T cell 4
VlnPlot(sce,features = c('CD3E','CCL5','CD69','GZMA','NKG7','CD8A','CD3D'),pt.size = 0,group.by = 'seurat_clusters')
# #Fibroblast   2
VlnPlot(sce,features = c('COL1A1','COL3A1','COL1A2'),pt.size = 0,group.by = 'seurat_clusters')
#Endothelial cell  3,5
VlnPlot(sce,features = c('EMCN','VWF','PLVAP'),pt.size = 0,group.by = 'seurat_clusters')
#Epithelial cell  0,6
VlnPlot(sce,features = c('CA9','EPCAM','KRT8'),pt.size = 0,group.by = 'seurat_clusters')
VlnPlot(sce,features = c('TOP2A','MKI67'),pt.size = 0,group.by = 'seurat_clusters')



#
Logfc = 0.25
#
Minpct = 0.25
colnames(sce@meta.data)
Idents(sce)<-'seurat_clusters'

sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
write.csv(sce.markers, "results/06.scRNA/All_cluster_type.csv")


####
marker <- data.frame(cluster = 0:6,cell = 0:6)
marker[marker$cluster %in% c(0),2] <- 'Cancer cells'
marker[marker$cluster %in% c(1),2] <- 'B cells'
marker[marker$cluster %in% c(2),2] <- 'Fibroblast'
marker[marker$cluster %in% c(3,5),2] <- 'Endothelial cells'
marker[marker$cluster %in% c(4),2] <- 'Natural killer T (NKT) cells'
marker[marker$cluster %in% c(6),2] <- 'Profiling Cancer cells'

marker
sce@meta.data$cell_type <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})

# library(cols4all)
colors = c("#E7298A","#FDB462","#FB8072","#80B1D3","#B3DE69","#A6761D")


seurat_clusters_umap=DimPlot(sce,group.by='seurat_clusters',reduction="tsne",label =F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12),legend.position = 'none')+
  ggtitle('seurat_clusters')
seurat_clusters_umap=LabelClusters(plot = seurat_clusters_umap, id = 'seurat_clusters', family = 'Times', size = 5, colour = "black", repel = T)
seurat_clusters_umap
ggsave("results/06.scRNA/UMAP_seurat_clusters.pdf", seurat_clusters_umap, width=7, height=6)


cell_type_umap=DimPlot(sce,group.by='cell_type',reduction="tsne",label =F,pt.size = 0.2,cols =colors)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12),legend.position = 'none')+
  ggtitle('Cell type')
cell_type_umap=LabelClusters(plot = cell_type_umap, id = 'cell_type', family = 'Times', size = 5, colour = "black", repel = T)
cell_type_umap
ggsave("results/06.scRNA/UMAP_cell_type.pdf", cell_type_umap, width=7, height=6)

genes = c('CD79A','MS4A1',
          'MKI67','TOP2A','CA9','EPCAM','KRT8',
          'EMCN','VWF','PLVAP',
          'COL1A1','COL3A1','COL1A2',
          # 'IL1B','C1QB','LYZ','G0S2',
          'GNLY','NKG7','GZMA')

dotplot_gene_marker=DotPlot(sce, features=genes,group.by = 'cell_type')+coord_flip()+
  scale_color_gradientn(colors=c( "dodgerblue", "white", "orange", "firebrick1"))+theme_light()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=14, face="bold"), 
        axis.text.y=element_text(face="bold", size=14), axis.title.x=element_blank(), axis.title.y=element_blank())
dotplot_gene_marker
ggsave("results/06.scRNA/dotplot_gene_marker.pdf", dotplot_gene_marker, width=7, height=6)


library(lessR)
colors = c("#B3DE69","#1B9E77",'#BC80BD',"#FB8072","#E7298A","#FCCDE5","#A6761D","#FDB462","#80B1D3")
# pdf("07.scRNA/Proportion_cell_type.pdf", width=7, height=6)
PieChart(cell_type, data=sce.celltype, hole=0.5, main="Cell Proportion", main_cex=1.3, fill=colors)
dev.off()



##hubgene############
DotPlot(sce, features=hub.genes,group.by = 'cell_type')+
  scale_color_gradientn(colors=c( "dodgerblue", "white", "orange", "firebrick1"))+theme_light()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=14, face="bold"), 
        axis.text.y=element_text(face="bold", size=14), axis.title.x=element_blank(), axis.title.y=element_blank())

fig6d=VlnPlot(sce, features = hub.genes,pt.size = 0,group.by = 'cell_type',log = T)+NoLegend()
pdf("results/06.scRNA/hub_vlnplot.pdf", width=7, height=6)
fig6d
dev.off()

save.image(file = 'project.RData')


fig6=mg_merge_plot(dotplot_gene_marker,dotplot_gene_marker,dotplot_gene_marker,fig6d,nrow=2,ncol=2,labels = LETTERS[1:4])
ggsave('results/06.scRNA/Fig6.pdf',fig6,height = 12,width = 12)
