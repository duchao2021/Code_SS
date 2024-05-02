library(GEOquery)
library(limma)
library(ggplot2)
library(dplyr)
library(forcats)
library(MASS)
library(pROC)
library(ROCR)#
library(ggfortify)
library(pheatmap)
library(glmnet)
library(DESeq2)
library(enrichplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(RColorBrewer)
library(momr)
library(patchwork)
library(ggpubr)
library(VennDiagram) 
library(ggridges)
library(GOplot)
library(grid)
library(cowplot)
library(ggbiplot)

#Violin
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}+scale_fill_manual(values = c("#7697CB","#D293C1"))

#List of GSE to use in the anlaysis
GEOs <- c("GSE95233","GSE57065","GSE54514")
#Load all GSEs and save to a local files for later easy load
for(i in GEOs){
  if(!file.exists(paste0("~/Desktop/Sepsis septic shock/2022 publication/Code/",i,".annot")) & !file.exists(paste0(i,".gset"))){
    gset <- getGEO(i,GSEMatrix =TRUE, getGPL = T,AnnotGPL = T)[[1]]
    annot <- fData(gset)
    save(annot,file = paste0(i,".annot"))
    save(gset,file = paste0(i,".gset"))
    write.csv(as.data.frame(pData(gset)),paste0(i," pData.csv"),row.names = T) }}
#Load gset and annot files to the environment
for(i in GEOs){
  load(paste0("~/Desktop/Sepsis septic shock/2022 publication/Code/",i,".gset"))
  assign(paste0("gset",i), gset)
  load(paste0("~/Desktop/Sepsis septic shock/2022 publication/Code/",i,".annot"))
  assign(paste0("annot",i),annot) }
#ADJUST GENE NAMES AND REMOVE DUPLICATED GENES
for(i in GEOs){
  annot <- get(paste0("annot",i))
  annot <- annot[annot$`Gene symbol` != "",]
  annot <- annot[grep("///",annot$`Gene symbol`,invert = T),]
  data <- exprs(get(paste0("gset",i)))[as.character(annot$ID),]
  rownames(data) <- annot$`Gene symbol`
  data <- data[!duplicated(rownames(data)),]
  data <- data[order(rownames(data)),]
  assign(paste0("data",i),data) }
#Free some space by removing the raw annot and temporary files
rm(list = c("data","annot","gset",paste0("annot",GEOs)))
#Gene expression profile
dat1=data.frame(dataGSE95233)
dat2=data.frame(dataGSE57065)
dat3=data.frame(dataGSE54514)
#Gene expression shared by the three datasets
dat1=dat1[rownames(dat1)%in%rownames(dat2),]
dat1=dat1[rownames(dat1)%in%rownames(dat3),]
dat2=dat2[rownames(dat2)%in%rownames(dat1),]
dat3=dat3[rownames(dat3)%in%rownames(dat1),]
#check normalization
boxplot(data.frame(dat1),las=2)
boxplot(data.frame(dat2),las=2)
boxplot(data.frame(dat3),las=2)
#Prepare Group information for each dataset
#GSE95233
pData1=read.csv(file.choose(),stringsAsFactors = F)
pData1=pData1[which(pData1$characteristics_ch1.2=="time point: NA"|
                      pData1$time.point.ch1=="D01"),]
pData1$Group[pData1$characteristics_ch1.2=="time point: NA"]=c("Control")
pData1$Group[pData1$characteristics_ch1.2=="time point: D01"]=c("Patient")
table(pData1$Group)
Group1=data.frame(pData1$Group)
colnames(Group1)[1]=c("Group")
Group1$id=pData1$X
data1=dat1[,colnames(dat1)%in%pData1$X,]
#GSE57065
pData2=read.csv(file.choose(),stringsAsFactors = F)
pData2=pData2[which(pData2$characteristics_ch1.2=="sapsii: NA"|
                      pData2$collection.time.ch1=="24 hr"),]
pData2$Group[pData2$characteristics_ch1.2=="sapsii: NA"]=c("Control")
pData2$Group[pData2$collection.time.ch1=="24 hr"]=c("Patient")
table(pData2$Group)
Group2=data.frame(pData2$Group)
colnames(Group2)[1]=c("Group")
Group2$id=pData2$X
data2=dat2[,colnames(dat2)%in%pData2$X,]
#GSE54514
pData3=read.csv(file.choose(),stringsAsFactors = F)
pData3=pData3[which(pData3$group_day.ch1=="HC_D1"|
                      pData3$group_day.ch1=="NS_D1"|
                      pData3$group_day.ch1=="S_D1"),]
pData3$Group[pData3$group_day.ch1=="HC_D1"]=c("Control")
pData3$Group[pData3$group_day.ch1=="NS_D1"|pData3$group_day.ch1=="S_D1"]=c("Patient")
table(pData3$Group)
Group3=data.frame(pData3$Group)
colnames(Group3)[1]=c("Group")
Group3$id=pData3$X
data3=dat3[,colnames(dat3)%in%pData3$X,]

CIBERSORTx
CIBERSORTx_1=read.csv(file.choose(),stringsAsFactors = F)
CIBERSORTx_2=read.csv(file.choose(),stringsAsFactors = F)
CIBERSORTx_3=read.csv(file.choose(),stringsAsFactors = F)
CIBERSORTx_4=read.csv(file.choose(),stringsAsFactors = F)
CIBER_1=CIBERSORTx_1[,(2:3)]
CIBER_1$Tcell=CIBERSORTx_1$T.cells.CD8+CIBERSORTx_1$T.cells.CD4.naive+CIBERSORTx_1$T.cells.CD4.memory.resting+
  CIBERSORTx_1$T.cells.CD4.memory.activated+CIBERSORTx_1$T.cells.follicular.helper+CIBERSORTx_1$T.cells.gamma.delta+CIBERSORTx_1$T.cells.regulatory..Tregs.
CIBER_1$CD4T=CIBERSORTx_1$T.cells.CD4.naive+CIBERSORTx_1$T.cells.CD4.memory.resting+CIBERSORTx_1$T.cells.CD4.memory.activated
CIBER_1$CD8T=CIBERSORTx_1$T.cells.CD8
CIBER_1$Bcell=CIBERSORTx_1$B.cells.naive+CIBERSORTx_1$B.cells.memory
CIBER_1$NKcell=CIBERSORTx_1$NK.cells.resting+CIBERSORTx_1$NK.cells.activated
CIBER_1$Neutro=CIBERSORTx_1$Neutrophils
CIBER_1$Mono=CIBERSORTx_1$Monocytes

CIBER_2=CIBERSORTx_2[,(2:3)]
CIBER_2$Tcell=CIBERSORTx_2$T.cells.CD8+CIBERSORTx_2$T.cells.CD4.naive+CIBERSORTx_2$T.cells.CD4.memory.resting+
  CIBERSORTx_2$T.cells.CD4.memory.activated+CIBERSORTx_2$T.cells.follicular.helper+CIBERSORTx_2$T.cells.gamma.delta+CIBERSORTx_2$T.cells.regulatory..Tregs.
CIBER_2$CD4T=CIBERSORTx_2$T.cells.CD4.naive+CIBERSORTx_2$T.cells.CD4.memory.resting+CIBERSORTx_2$T.cells.CD4.memory.activated
CIBER_2$CD8T=CIBERSORTx_2$T.cells.CD8
CIBER_2$Bcell=CIBERSORTx_2$B.cells.naive+CIBERSORTx_2$B.cells.memory
CIBER_2$NKcell=CIBERSORTx_2$NK.cells.resting+CIBERSORTx_2$NK.cells.activated
CIBER_2$Neutro=CIBERSORTx_2$Neutrophils
CIBER_2$Mono=CIBERSORTx_2$Monocytes

CIBER_3=CIBERSORTx_3[,(2:3)]
CIBER_3$Tcell=CIBERSORTx_3$T.cells.CD8+CIBERSORTx_3$T.cells.CD4.naive+CIBERSORTx_3$T.cells.CD4.memory.resting+
  CIBERSORTx_3$T.cells.CD4.memory.activated+CIBERSORTx_3$T.cells.follicular.helper+CIBERSORTx_3$T.cells.gamma.delta+CIBERSORTx_3$T.cells.regulatory..Tregs.
CIBER_3$CD4T=CIBERSORTx_3$T.cells.CD4.naive+CIBERSORTx_3$T.cells.CD4.memory.resting+CIBERSORTx_3$T.cells.CD4.memory.activated
CIBER_3$CD8T=CIBERSORTx_3$T.cells.CD8
CIBER_3$Bcell=CIBERSORTx_3$B.cells.naive+CIBERSORTx_3$B.cells.memory
CIBER_3$NKcell=CIBERSORTx_3$NK.cells.resting+CIBERSORTx_3$NK.cells.activated
CIBER_3$Neutro=CIBERSORTx_3$Neutrophils
CIBER_3$Mono=CIBERSORTx_3$Monocytes

CIBER_4=CIBERSORTx_4[,(2:3)]
CIBER_4$Tcell=CIBERSORTx_4$T.cells.CD8+CIBERSORTx_4$T.cells.CD4.naive+CIBERSORTx_4$T.cells.CD4.memory.resting+
  CIBERSORTx_4$T.cells.CD4.memory.activated+CIBERSORTx_4$T.cells.follicular.helper+CIBERSORTx_4$T.cells.gamma.delta+CIBERSORTx_4$T.cells.regulatory..Tregs.
CIBER_4$CD4T=CIBERSORTx_4$T.cells.CD4.naive+CIBERSORTx_4$T.cells.CD4.memory.resting+CIBERSORTx_4$T.cells.CD4.memory.activated
CIBER_4$CD8T=CIBERSORTx_4$T.cells.CD8
CIBER_4$Bcell=CIBERSORTx_4$B.cells.naive+CIBERSORTx_4$B.cells.memory
CIBER_4$NKcell=CIBERSORTx_4$NK.cells.resting+CIBERSORTx_4$NK.cells.activated
CIBER_4$Neutro=CIBERSORTx_4$Neutrophils
CIBER_4$Mono=CIBERSORTx_4$Monocytes

write.csv(CIBER_1,"GSE95233_CIBERSORT.csv")
write.csv(CIBER_2,"GSE57065_CIBERSORT.csv")
write.csv(CIBER_3,"GSE54514_CIBERSORT.csv")
write.csv(CIBER_4,"GSE154918_CIBERSORT.csv")

#DEG analysis of the three studies
group_list=list(Group1,Group2,Group3)
N=c(1,2,3)
List1=list(data1,data2,data3)
for (i in N) {
  design <- model.matrix(~0+factor(data.frame(group_list[i])$Group))
  colnames(design) <- levels(factor(data.frame(group_list[i])$Group))
  rownames(design) <- colnames(data.frame(List1[i]))
  fit <- lmFit(data.frame(List1[i]),design)
  cont.matrix <- makeContrasts("Patient-Control",levels = design)
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)
  tempOutput <- topTable(fit2,coef = 1,n=Inf,adjust="BH")
  nrDEG <- na.omit(tempOutput)
  nrDEG$gene=row.names(nrDEG)
  nrDEGs=nrDEG[nrDEG$adj.P.Val<0.05&abs(nrDEG$logFC)>0.6,]
  assign(paste0("nrDEG",i),nrDEG)
  assign(paste0("nrDEGs",i),nrDEGs)
}

save(data1,data2,data3,Group1,Group2,Group3,nrDEG1,nrDEG2,nrDEG3,nrDEG_shared,file = "Microarray_data.Rdata")

#Septic shock but not sepsis-related genes:DEGs shared by GSE95233 and GSE57065 but not GSE54514
#Venn plot
vcol=c("goldenrod1","skyblue","palegreen")
Venn=venn.diagram(list(" "=rownames(nrDEGs1),
                       " "=rownames(nrDEGs2),
                       " "=rownames(nrDEGs3)),
                  filename = NULL,lwd=1,cex=3.5,cat.cex=2.5,
                  fill=vcol,alpha=0.7,margin=0.1)
grid.draw(Venn)

nrDEG_shared=nrDEGs1[rownames(nrDEGs1)%in%rownames(nrDEGs2),]
nrDEG_shared=nrDEG_shared[!rownames(nrDEG_shared)%in%rownames(nrDEGs3),]
nrDEG_shared=data.frame(rownames(nrDEG_shared))
colnames(nrDEG_shared)=c("Gene")
nrDEG1_shared=nrDEG1[rownames(nrDEG1)%in%nrDEG_shared$Gene,]
nrDEG1_shared_down=nrDEG1_shared[nrDEG1_shared$logFC<0,]
nrDEG2_shared=nrDEG2[rownames(nrDEG2)%in%nrDEG_shared$Gene,]
nrDEG2_shared_down=nrDEG2_shared[nrDEG2_shared$logFC<0,]
nrDEG_shared_down=data.frame(rownames(nrDEG1_shared_down[rownames(nrDEG1_shared_down)%in%rownames(nrDEG2_shared_down),]))
colnames(nrDEG_shared_down)=c("Gene")
nrDEG_shared_up=data.frame(nrDEG_shared[!nrDEG_shared$Gene%in%nrDEG_shared_down$Gene,])
colnames(nrDEG_shared_up)=c("Gene")
nrDEG3_shared=nrDEG3[rownames(nrDEG3)%in%nrDEG_shared$Gene,]

#Enrichment analysis: fgsea
List2_1=list(nrDEG1_shared,nrDEG2_shared,nrDEG3_shared)
N
for (i in N) {
  df <- data.frame(data.frame(List2_1[i])$logFC)
  df$SYMBOL <- rownames(data.frame(List2_1[i]))
  colnames(df)[1] <- c("logFC")
  df_id <- bitr(df$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  df_all <- merge(df,df_id,by="SYMBOL",all=F)
  assign(paste0("df_all_",i),df_all)
  df_all_sort <- df_all[order(df_all$logFC,decreasing = T),]
  gene_fc <- df_all_sort$logFC
  names(gene_fc) <- df_all_sort$ENTREZID
  #GSEA analysis
  KEGG <- gseKEGG(gene_fc,organism = "hsa")
  assign(paste0("KEGG_GSEA_",i),KEGG)
  #GO enrichment
  GO <- gseGO(gene_fc,ont = "ALL",OrgDb = org.Hs.eg.db,keyType = "ENTREZID",
              pvalueCutoff = 0.05,pAdjustMethod = "BH", verbose = TRUE, seed = FALSE)#BP,MF,CC,ALL
  assign(paste0("GO_GSEA_",i),GO)
  #fgsea
  KEGG_f <- gseKEGG(gene_fc,organism = "hsa",keyType = "kegg",exponent = 1,minGSSize = 10,
                    maxGSSize = 500,eps = 1e-10,pvalueCutoff = 0.05,pAdjustMethod = "BH",
                    verbose = TRUE,use_internal_data = FALSE,seed = FALSE,by = "fgsea")
  assign(paste0("KEGG_fgsea_",i),KEGG_f)
  sortKEGG_f<-KEGG_f[order(KEGG_f$enrichmentScore, decreasing = F),]
  assign(paste0("sortKEGG_fgsea_",i),sortKEGG_f)
}
#NK cell-related pathway
gseaplot2(KEGG_fgsea_1, "hsa04650", color = "blue", rel_heights=c(1, .2, .6),
title = "          KEGG pathway: Natural killer cell mediated cytotoxicity
          ",base_size = 18)#1000*460
gseaplot2(KEGG_fgsea_2, "hsa04650", color = "blue", rel_heights=c(1, .2, .6),
          title = "          KEGG pathway: Natural killer cell mediated cytotoxicity
          ",base_size = 18)#1000*460
#T cell-related pathway
gseaplot2(KEGG_fgsea_1, "hsa04660", color = "blue", rel_heights=c(1, .2, .6),
          title = "",base_size = 18)#707*517
gseaplot2(KEGG_fgsea_2, "hsa04660", color = "blue", rel_heights=c(1, .2, .6),
          title = "",base_size = 18)

Gene_NK_1 <- data.frame(unlist(strsplit(sortKEGG_fgsea_1[15,]$core_enrichment,split = '/')))
colnames(Gene_NK_1)[1] <- c("ENTREZID")
Gene_NK_1 <- merge(Gene_NK_1,df_all_1,by="ENTREZID",all=F)
Gene_NK_1 <- Gene_NK_1[order(Gene_NK_1$logFC,decreasing = F),]
Gene_NK_2 <- data.frame(unlist(strsplit(sortKEGG_fgsea_2[19,]$core_enrichment,split = '/')))
colnames(Gene_NK_2)[1] <- c("ENTREZID")
Gene_NK_2 <- merge(Gene_NK_2,df_all_2,by="ENTREZID",all=F)
Gene_NK_2 <- Gene_NK_2[order(Gene_NK_2$logFC,decreasing = F),]

Gene_T_1 <- data.frame(unlist(strsplit(sortKEGG_fgsea_1[13,]$core_enrichment,split = '/')))
colnames(Gene_T_1)[1] <- c("ENTREZID")
Gene_T_1 <- merge(Gene_T_1,df_all_1,by="ENTREZID",all=F)
Gene_T_1 <- Gene_T_1[order(Gene_T_1$logFC,decreasing = F),]
Gene_T_2 <- data.frame(unlist(strsplit(sortKEGG_fgsea_2[14,]$core_enrichment,split = '/')))
colnames(Gene_T_2)[1] <- c("ENTREZID")
Gene_T_2 <- merge(Gene_T_2,df_all_2,by="ENTREZID",all=F)
Gene_T_2 <- Gene_T_2[order(Gene_T_2$logFC,decreasing = F),]

#fgsea pathways
fgRes1s <- read.csv(file.choose(),stringsAsFactors = F)
fgRes1s$Enrichment = ifelse(fgRes1s$NES > 0, "Up-regulated", "Down-regulated")
#fgRes1s_up=fgRes1s[which(fgRes1s$Enrichment=="Up-regulated"&fgRes1s$NES>1.722),]#keep top 10
#fgRes1s_down=fgRes1s[which(fgRes1s$Enrichment=="Down-regulated"&abs(fgRes1s$NES)>1.65),]#keep top 10
filtRes = fgRes1s
ggplot(filtRes, aes(reorder(Description, ID), NES)) +
  geom_segment( aes(reorder(Description, ID), xend=Description, y=0, yend=NES)) +
  geom_point( size=5, aes( fill = Enrichment),
              shape=21, stroke=2) +
  scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                               "Up-regulated" = "firebrick") ) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + #,title="KEGG (Septic shock vs Control)"
  theme_minimal()+
  theme(axis.text.x = element_text(face="bold", size=16,color = "black"),
        axis.text.y = element_text(face="bold", size=16,color = "black"),
        axis.title.x = element_text(face="bold", size=16,color = "black"),
        axis.title.y = element_text(face="bold", size=16,color = "black"))#color="blue",,hjust=0.2,lineheight=0.2;,title=element_text(size=16,face="bold")
#tiff size 860W * 400H
fgRes2s <- read.csv(file.choose(),stringsAsFactors = F)
fgRes2s$Enrichment = ifelse(fgRes2s$NES > 0, "Up-regulated", "Down-regulated")
#fgRes1s_up=fgRes1s[which(fgRes1s$Enrichment=="Up-regulated"&fgRes1s$NES>1.722),]#keep top 10
#fgRes1s_down=fgRes1s[which(fgRes1s$Enrichment=="Down-regulated"&abs(fgRes1s$NES)>1.65),]#keep top 10
filtRes = fgRes2s
ggplot(filtRes, aes(reorder(Description, ID), NES)) +
  geom_segment( aes(reorder(Description, ID), xend=Description, y=0, yend=NES)) +
  geom_point( size=5, aes( fill = Enrichment),
              shape=21, stroke=2) +
  scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                               "Up-regulated" = "firebrick") ) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + #,title="KEGG (Septic shock vs Control)"
  theme_minimal()+
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16))#color="blue",,hjust=0.2,lineheight=0.2;,title=element_text(size=16,face="bold")
#tiff size 860W * 400H





#Hlty = "gray90", Seps_P = "green3", Shock_P = "dodgerblue1" #darkslateblue gray1

#heatmap_GSE95233
Group1=Group1[order(Group1$Group,decreasing = F),]
ann1=Group1[,1]
ann1=data.frame(ann1)
rownames(ann1)=Group1$id
colnames(ann1)[1]=c("Group")
ann1$Group[ann1$Group=="Patient"]=c("Septic_shock")
rts1=data1[as.character(Gene_NK_1$SYMBOL),]
rts1=rts1[,as.character(rownames(ann1))]
ann_colors = list(Group = c(Control = "gray50", Septic_shock = "dodgerblue3"))
pheatmap(rts1, annotation=ann1, 
         color = colorRampPalette(c("purple4", "white", "firebrick1"))(50),
         cluster_cols =F,cluster_rows = F,
         fontsize = 20,
         fontsize_row=18,
         fontsize_col=5,
         #fontface = "bold",
         border_color = NA,
         cellwidth = 4.95,
         cellheight = 28,
         annotation_colors = ann_colors,
         show_colnames = F)
#962*596
#heatmap_GSE57065
Group2=Group2[order(Group2$Group,decreasing = F),]
ann2=Group2[,1]
ann2=data.frame(ann2)
rownames(ann2)=Group2$id
colnames(ann2)[1]=c("Group")
ann2$Group[ann2$Group=="Patient"]=c("Septic_shock")
rts2=data2[as.character(Gene_NK_2$SYMBOL),]
rts2=rts2[,as.character(rownames(ann2))]
ann_colors = list(Group = c(Control = "gray90", Septic_shock = "dodgerblue3"))
pheatmap(rts2, annotation=ann2, 
         color = colorRampPalette(c("purple4", "white", "firebrick1"))(50),
         cluster_cols =F,cluster_rows = F,
         fontsize = 20,
         fontsize_row=18,
         fontsize_col=5,
         #fontface = "bold",
         border_color = NA,
         cellwidth = 6.5,
         cellheight = 28,
         annotation_colors = ann_colors,
         show_colnames = F)
#962*596

Gene_NK <- Gene_NK_1[Gene_NK_1$SYMBOL%in%Gene_NK_2$SYMBOL,]$SYMBOL
Gene_NK_1_1.5 <- Gene_NK_1[abs(Gene_NK_1$logFC)>1.5,]
Gene_NK_2_1.5 <- Gene_NK_2[abs(Gene_NK_2$logFC)>1.5,]
Gene <- Gene_NK_1_1.5[Gene_NK_1_1.5$SYMBOL%in%Gene_NK_2_1.5$SYMBOL,]$SYMBOL
Gene <- data.frame(Gene)

exprSet=`GSE154918_GeneLevel_Normalized(CPM.and.TMM)_data`
exprSet_VDR=exprSet[exprSet$gene_symbol%in%"VDR",]
rownames(exprSet_VDR)=exprSet_VDR[,1]
exprSet_VDR=exprSet_VDR[,-1]
data4_1=exprSet_VDR[,colnames(exprSet_VDR)%in%pData4_H_S$id]
data4_2=exprSet_VDR[,colnames(exprSet_VDR)%in%pData4_H_SS$id]

exprSet_GZMB=exprSet[exprSet$gene_symbol%in%"GZMB",]
data4_VDR=data4[rownames(data4)%in%"VDR",]

VDR1=data1[rownames(data1)%in%"VDR",]
VDR1=data.frame(t(VDR1))
VDR1$id=rownames(VDR1)
VDR1=merge(VDR1,Group1,by='id')
VDR1$value=VDR1$VDR/mean(VDR1$VDR[VDR1$Group=="Control"])
VDR2=data2[rownames(data2)%in%"VDR",]
VDR2=data.frame(t(VDR2))
VDR2$id=rownames(VDR2)
VDR2=merge(VDR2,Group2,by='id')
VDR2$value=VDR2$VDR/mean(VDR2$VDR[VDR2$Group=="Control"])
VDR3=data3[rownames(data3)%in%"VDR",]
VDR3=data.frame(t(VDR3))
VDR3$id=rownames(VDR3)
VDR3=merge(VDR3,Group3,by='id')
VDR3$value=VDR3$VDR/mean(VDR3$VDR[VDR3$Group=="Control"])
VDR4_1=data4_1
VDR4_1=data.frame(t(VDR4_1))
VDR4_1$id=rownames(VDR4_1)
VDR4_1=merge(VDR4_1,pData4_H_S,by='id')
VDR4_1$value=VDR4_1$VDR/mean(VDR4_1$VDR[VDR4_1$Group=="Control"])
VDR4_2=data4_2[rownames(data4_2)%in%"VDR",]
VDR4_2=data.frame(t(VDR4_2))
VDR4_2$id=rownames(VDR4_2)
VDR4_2=merge(VDR4_2,pData4_H_SS,by='id')
VDR4_2$value=VDR4_2$VDR/mean(VDR4_2$VDR[VDR4_2$Group=="Control"])


#Prediction models based on the 6 genes
 #LDA
#List1
#group_list
#N
for (i in N) {
  data=data.frame(List1[i])
  data=data[rownames(data)%in%Gene$Gene,]
  data=data.frame(t(data))
  data$id=rownames(data)
  data=merge(data,data.frame(group_list[i]),by='id')
  rownames(data)=data$id
  data=data[,!colnames(data)%in%c("id")]
  nsample=nrow(data)
  ngenes=ncol(data[,!colnames(data)%in%c("Group")])
  data.df=data
  # split data into training and test sets
  set.seed(0)
  train.ix = sample(1 : nsample, round(nsample / 2))
  test.ix = setdiff(1 : nsample, train.ix)
  train.df = data.df[train.ix, ]
  test.df = data.df[test.ix, ]
  # check class proportions in training and test sets
  table(train.df$Group); table(test.df$Group)
  # run LDA
  ldamodel = lda(Group ~ ., data = train.df)
  # perform classification for train set
  pred.id0 = predict(ldamodel, train.df)$class
  test.id0 = train.df$Group
  mcr0 = mean(pred.id0 != test.id0)
  result0=data.frame(table(test.id0, pred.id0))
  assign(paste0("LDA_train_result_",i),result0)
  # perform classification for test set
  pred.id = predict(ldamodel, test.df)$class
  test.id = test.df$Group
  mcr = mean(pred.id != test.id)
  result=data.frame(table(test.id, pred.id))
  assign(paste0("LDAresult_",i),result)
  #ROC curve
  pred<-predict(ldamodel, test.df,type='prob')
  pred1<-as.data.frame(pred$posterior)
  #test.df$group[test.df$Group=="Control"]=0
  #test.df$group[test.df$Group=="Patient"]=1
  #pred2<-prediction(pred1[,2],test.df$group)
  #roc.perf<-performance(pred2,measure="tpr",x.measure="fpr")
  #auc.train<-performance(pred2,measure="auc")
  #auc.train<-auc.train@y.values
  #yourfilename=paste("ROCtesting",i,".pdf",sep="")
  #pdf(file=yourfilename,width = 5,height = 5)
  #plot(roc.perf)
  #abline(a=0,b=1)
  #text(x=.25,y=.65,paste("AUC = ",round(auc.train[[1]],3),sep=""))
  #dev.off()
  modelroc<-roc(test.df$Group,pred1[,2])
  yourfilename=paste("ROCLDA_",i,".tiff",sep="")
  #pdf(file=yourfilename,width = 5,height = 5)
  tiff(file=yourfilename,width = 1600,height = 1600,res = 300)
  plot(modelroc,print.auc=FALSE, print.auc.x=0.2,print.auc.y=0.2,cex.lab=2,cex.axis=2,cex.main=3,
       font.lab=2,#font=2,
       max.auc.polygon=T, grid=c(0.5, 0.5),grid.col=c("green", "red"),max.auc.polygon.col="white",
       main="",xlab="False Positive Rate",ylab="True Positive Rate",col="black",legacy.axes=TRUE,lwd=8)
  dev.off()
  ##LD1 of the three studies
  linear=lda(Group~.,data)
  p=predict(linear,data)
  posterior=data.frame(p$posterior)
  #assign(paste0("posterior_",i),posterior)
  LD1=p$x[,1]
  LD1=data.frame(LD1)
  LD1$Group=data$Group
  LD1$id=rownames(data)
  assign(paste0("LD1value_",i),LD1)
  roc=roc(LD1$Group,LD1$LD1)
  best=coords(roc,"best",ret = "all",transpose = F)
  assign(paste0("LD1_roc_",i),roc)
  assign(paste0("LD1_best_",i),best)
  yourfilename=paste("ROCLD1_",i,".tiff",sep="")
  #pdf(file=yourfilename,width = 5,height = 5)
  tiff(file=yourfilename,width = 1600,height = 1600,res = 300)
  plot(roc,print.auc=FALSE, print.auc.x=0.2,print.auc.y=0.2,cex.lab=2,cex.axis=2,cex.main=3,
       font.lab=2,#font=2,
       max.auc.polygon=T, grid=c(0.5, 0.5),grid.col=c("green", "red"),max.auc.polygon.col="white",
       main="",xlab="False Positive Rate",ylab="True Positive Rate",col="black",legacy.axes=TRUE,lwd=8)
  legend("bottomright",c(paste("AUC=",round(roc$auc,3),sep = ""),
                         paste("Sensitivity=",round(best$sensitivity,3),sep = ""),
                         paste("Specificity=",round(best$specificity,3),sep = ""),
                         paste("Accuracy=",round(best$accuracy,3),sep = "")),col = c("red"),
         #lwd = 1,
         bty = "n",cex = 2)
  dev.off()
}
aLD1value_1=LD1value_1
aLD1value_2=LD1value_2
aLD1value_3=LD1value_3
LD1value_1$Group[LD1value_1$Group=="Patient"]=c("Septic shock")
LD1value_2$Group[LD1value_2$Group=="Patient"]=c("Septic shock")
LD1value_3$Group[LD1value_3$Group=="Patient"]=c("Sepsis")
ggplot(LD1value_1,aes(x=LD1,fill=Group))+
  geom_histogram(binwidth = .2,alpha=.5)+
  geom_count()
#
theme_rbook <- function(base_size = 24, base_family = "", base_line_size = base_size/24,base_rect_size = base_size/24) #bold
{theme(text=element_text(family = "Arial"),axis.title = element_text(size = 24),axis.text.x = element_text(size = 24,colour = "black"),                              
       axis.text.y = element_text(size = 24,colour = "black"),plot.caption = element_text(size = 24, face = "italic"),            
       panel.background = element_rect(fill="white"),axis.line = element_line(size = 1, colour = "black"),
       axis.ticks.y = element_line(size = 1, colour = "black"),legend.title=element_text(size=24) , legend.text=element_text(size=24),
       legend.position = c(0.70,0.8),
       strip.background =element_rect(fill = "#cddcdd"),
       #panel.border = element_rect(colour = "black", fill=NA, size=1),
       strip.text = element_text(colour = "black"),legend.key=element_blank())}

 #Risk score model
List4=list(data.frame(t(data1)),data.frame(t(data2)),data.frame(t(data3)))
#List2
#group_list
#N
List2=List2_1
N0=c(1:6)
for (n in N) {
  risk=data.frame(List4[n])
  Sum=data.frame(group_list[n])
  rownames(Sum)=Sum$id
  for(i in Gene$Gene){
    GAPDH=risk[,colnames(risk)==c("GAPDH")]
    GAPDH=data.frame(GAPDH)
    rownames(GAPDH)=rownames(risk)
    data=data.frame(risk[,colnames(risk)==c(i)])
    rownames(data)=rownames(risk)
    colnames(data)[1]=c(i)
    data$relative=data[,1]/GAPDH$GAPDH
    Mean=mean(data$relative)
    Std=sd(data$relative)
    Weight=ifelse(data.frame(List2[n])[rownames(data.frame(List2[n]))==c(i),]$logFC<0,-1,1)
    data$riskscore=(data$relative-Mean)/Std*Weight
    assign(paste0("data",i),data)
    for (k in N0) {
      Sum[,N0[k]]=data[,3]
    }
    colnames(Sum)=Gene$Gene
    assign(paste0("Riskscore_Sum_",n),Sum)
  }
  RScore=Sum[,1]+Sum[,2]+Sum[,3]+Sum[,4]+Sum[,5]+Sum[,6]
  RScore=data.frame(RScore)
  RScore$id=rownames(risk)
  Group=data.frame(group_list[n])
  RScore=merge(RScore,Group,by='id')
  roc=roc(RScore$Group,RScore$RScore)
  best=coords(roc,"best",ret = "all",transpose = F)
  assign(paste0("RScore_roc_",n),roc)
  assign(paste0("RScore_best_",n),best)
  assign(paste0("RScore_",n),RScore)
  yourfilename=paste("ROCRisk_",n,".tiff",sep="")
  #pdf(file=yourfilename,width = 5,height = 5)
  tiff(file=yourfilename,width = 1600,height = 1600,res = 300)
  plot(roc,print.auc=FALSE, print.auc.x=0.2,print.auc.y=0.2,cex.lab=2,cex.axis=2,cex.main=3,
      font.lab=2,#font=2,
      max.auc.polygon=T, grid=c(0.5, 0.5),grid.col=c("green", "red"),max.auc.polygon.col="white",
      main="",xlab="False Positive Rate",ylab="True Positive Rate",col="black",legacy.axes=TRUE,lwd=8)
  
  dev.off()
}


#density plot of risk scores
RScore_1$Group[RScore_1$Group=="Patient"]=c("SS")
RScore_2$Group[RScore_2$Group=="Patient"]=c("SS")
RScore_3$Group[RScore_3$Group=="Patient"]=c("Sep")
RScore_1$Group[RScore_1$Group=="Control"]=c("HC")
RScore_2$Group[RScore_2$Group=="Control"]=c("HC")
RScore_3$Group[RScore_3$Group=="Control"]=c("HC")
List5=list(RScore_1[,-1],RScore_2[,-1],RScore_3[,-1])
#N
Dataset=c("GSE95233       ","GSE57065       ","GSE54514      ")
N=3
for (i in N) {
  data=data.frame(List5[i])
  yourfilename=paste("Fig4_Density_6genes_RS-",Dataset[i],".tiff",sep="")
  tiff(file=yourfilename,width = 2000,height = 1600,res = 300)
  print(ggplot(data)+
          geom_density(aes(x=RScore,fill=Group),alpha=0.5)+
          scale_x_continuous(limits=c(-20, 20), breaks=c(-20,-10,0,10,20))+
          scale_y_continuous(limits=c(0, 0.2), breaks=c(0,0.05,0.1,0.15,0.2))+
          labs(y="density",x="Risk score",fill="Group"
               #,title = Dataset[i]
               )+
          scale_fill_manual(
            values = c("gray90","green4"))+
          #geom_vline(xintercept = -2.035,lty=2,lwd=1.5)+
          #geom_vline(xintercept = 0.703,lty=3,lwd=1.5)+
          #theme(plot.title = element_text(size=33,face="bold",hjust = 0.5)) +
          theme(legend.position="none")+
          theme_rbook())
  dev.off()
  }

N=c(1,2)
for (i in N) {
  data=data.frame(List5[i])
  yourfilename=paste("Fig4_Density_6genes_RS-",Dataset[i],".tiff",sep="")
  tiff(file=yourfilename,width = 2000,height = 1600,res = 300)
  print(ggplot(data)+
          geom_density(aes(x=RScore,fill=Group),alpha=0.5)+
          scale_x_continuous(limits=c(-20, 20), breaks=c(-20,-10,0,10,20))+
          scale_y_continuous(limits=c(0, 0.2), breaks=c(0,0.05,0.1,0.15,0.2))+
          labs(y="density",x="Risk score",fill="Group"
               #,title = Dataset[i]
          )+
          scale_fill_manual(
            values = c("gray90","dodgerblue3"))+
          #geom_vline(xintercept = -2.035,lty=2,lwd=1.5)+
          #geom_vline(xintercept = 0.703,lty=3,lwd=1.5)+
          #theme(plot.title = element_text(size=33,face="bold",hjust = 0.5)) +
          theme(legend.position="none")+
          theme_rbook())
  dev.off()
}

#GSE154918
data4=GSE154918.gene.expression
rownames(data4)=data4[,1]
data4=data4[,-1]
pData4=data.frame(GSE154918_full_metadata$status.ch1)
pData4$id=GSE154918_full_metadata$X
colnames(pData4)[1]=c("Group")
pData4=pData4[which(pData4$Group=="Hlty"|pData4$Group=="Shock_P"|pData4$Group=="Seps_P"),]
pData4_H_S=pData4[which(pData4$Group=="Hlty"|pData4$Group=="Seps_P"),]#Sepsis vs Control
pData4_H_S$Group[pData4_H_S$Group=="Hlty"]=c("Control")
pData4_H_S$Group[pData4_H_S$Group=="Seps_P"]=c("Patient")
pData4_H_SS=pData4[which(pData4$Group=="Hlty"|pData4$Group=="Shock_P"),]#Septic shock vs Control
pData4_H_SS$Group[pData4_H_SS$Group=="Hlty"]=c("Control")
pData4_H_SS$Group[pData4_H_SS$Group=="Shock_P"]=c("Patient")
pData4_S_SS=pData4[which(pData4$Group=="Seps_P"|pData4$Group=="Shock_P"),]#Septic shock vs Sepsis
pData4_S_SS$Group[pData4_S_SS$Group=="Seps_P"]=c("Control")
pData4_S_SS$Group[pData4_S_SS$Group=="Shock_P"]=c("Patient")
pData4_full=GSE154918_full_metadata[GSE154918_full_metadata$X%in%pData4$id,]
table(pData4_full[which(pData4_full$status.ch1=="Hlty"),]$Sex.ch1)
table(pData4_full[which(pData4_full$status.ch1=="Seps_P"),]$Sex.ch1)
table(pData4_full[which(pData4_full$status.ch1=="Shock_P"),]$Sex.ch1)
data4=data4[,colnames(data4)%in%pData4$id]
data4_H_S=data4[,colnames(data4)%in%pData4_H_S$id]
data4_H_SS=data4[,colnames(data4)%in%pData4_H_SS$id]
data4_S_SS=data4[,colnames(data4)%in%pData4_S_SS$id]

List6=list(pData4_H_S,pData4_H_SS,pData4_S_SS)
List7=list(data4_H_S,data4_H_SS,data4_S_SS)

#DEG analysis
N
for (i in N) {
  (colData=data.frame(row.names = colnames(data.frame(List7[i])),group_list=data.frame(List6[i])$Group))
  dds=DESeqDataSetFromMatrix(countData = data.frame(List7[i]),
                             colData = colData,
                             design = ~ group_list)
  dds=DESeq(dds)
  res=results(dds,
              contrast = c("group_list","Patient","Control"))
  resOrdered=res[order(res$padj),]
  DEG=as.data.frame(resOrdered)
  DEG=na.omit(DEG)
  DEG=DEG[DEG$baseMean>20,]
  DEG$gene=rownames(DEG)#####final gene list####
  assign(paste0("nrDEG4_",i),DEG)
}
save(data1,data2,data3,data4_H_S,data4_H_SS,data4_S_SS,Group1,Group2,Group3,pData4_H_S,pData4_H_SS,pData4_S_SS,
     nrDEG1,nrDEG2,nrDEG3,nrDEG4_1,nrDEG4_2,nrDEG4_3,nrDEG_shared,file = "S_1st_part_raw data.Rdata")

#heatmap_GSE57065
data4_N=`GSE154918_GeneLevel_Normalized(CPM.and.TMM)_data`
data4_N6=data4_N[data4_N$gene_symbol%in%Gene$Gene,]
rownames(data4_N6)=data4_N6$gene_symbol
data4_N6=data4_N6[,-1]
data4_N6=log2(data4_N6)
Group4=pData4_full[,(1:2)]
Group4$Group=pData4_full$status.ch1
Group4=Group4[order(Group4$Group,decreasing = F),]
ann4=Group4[,3]
ann4=data.frame(ann4)
rownames(ann4)=Group4$X
colnames(ann4)[1]=c("Group")
#ann2$Group[ann2$Group=="Patient"]=c("Septic_shock")
rts4=data4_N6[as.character(Gene$Gene),]
rts4=rts4[,as.character(rownames(ann4))]
ann_colors = list(Group = c(Hlty = "gray60", Seps_P = "green3", Shock_P = "dodgerblue3"))
pheatmap(rts4, annotation=ann4, 
         color = colorRampPalette(c("purple4", "white", "firebrick1"))(50),
         cluster_cols =F,cluster_rows = F,
         fontsize = 20,
         fontsize_row=18,
         fontsize_col=5,
         #fontface = "bold",
         border_color = NA,
         cellwidth = 6.5,
         cellheight = 28,
         annotation_colors = ann_colors,
         show_colnames = F)
#962*596



#LDA model
for (i in N) {
  data=data.frame(List7[i])
  data=data[rownames(data)%in%Gene$Gene,]
  data=data.frame(t(data))
  data$id=rownames(data)
  data=merge(data,data.frame(List6[i]),by='id')
  rownames(data)=data$id
  data=data[,!colnames(data)%in%c("id")]
  nsample=nrow(data)
  ngenes=ncol(data[,!colnames(data)%in%c("Group")])
  data.df=data
  # split data into training and test sets
  set.seed(0)
  train.ix = sample(1 : nsample, round(nsample / 2))
  test.ix = setdiff(1 : nsample, train.ix)
  train.df = data.df[train.ix, ]
  test.df = data.df[test.ix, ]
  # check class proportions in training and test sets
  table(train.df$Group); table(test.df$Group)
  # run LDA
  ldamodel = lda(Group ~ ., data = train.df)
  # perform classification for test set
  #pred.id = predict(ldamodel, test.df)$class
  #test.id = test.df$Group
  #mcr = mean(pred.id != test.id)
  #result=data.frame(table(test.id, pred.id))
  pred.id = predict(ldamodel, train.df)$class
  train.id = train.df$Group
  mcr = mean(pred.id != train.id)
  result=data.frame(table(train.id, train.id))
  assign(paste0("LDAresult_4_",i),result)
  ##LD1 of the three studies
  linear=lda(Group~.,data)
  p=predict(linear,data)
  LD1=p$x[,1]
  LD1=data.frame(LD1)
  LD1$Group=data$Group
  LD1$id=rownames(data)
  assign(paste0("LD1value_4_",i),LD1)
  roc=roc(LD1$Group,LD1$LD1)
  best=coords(roc,"best",ret = "all",transpose = F)
  assign(paste0("LD1_roc_4_",i),roc)
  assign(paste0("LD1_best_4_",i),best)
  yourfilename=paste("ROCLD1_4_",i,".tiff",sep="")
  #pdf(file=yourfilename,width = 5,height = 5)
  tiff(file=yourfilename,width = 1600,height = 1600,res = 300)
  plot(roc,print.auc=FALSE, print.auc.x=0.2,print.auc.y=0.2,cex.lab=2,cex.axis=2,cex.main=3,
       font.lab=2,#font=2,
       max.auc.polygon=T, grid=c(0.5, 0.5),grid.col=c("green", "red"),max.auc.polygon.col="white",
       main="",xlab="False Positive Rate",ylab="True Positive Rate",col="black",legacy.axes=TRUE,lwd=8)
  dev.off()
}
#ROC curve of LD1 values
tiff("Fig4_6genes_ROC_LD1value_GSE154918.tiff",width = 1600,height = 1600,res = 300)
plot(LD1_roc_4_1,print.auc=FALSE, print.auc.x=0.2,print.auc.y=0.2,cex.lab=2,cex.axis=2,cex.main=2,font.lab=2,
     max.auc.polygon=T, grid=c(0.5, 0.5),grid.col=c("green", "red"),max.auc.polygon.col="white",
     main="ROC curves of LD1 values",xlab="False Positive Rate",ylab="True Positive Rate",col="purple",legacy.axes=TRUE,lwd=8)
plot(LD1_roc_4_2,add=T,col="palegreen",lwd=8,lty=3)
plot(LD1_roc_4_3,add=T,col="red",lwd=8,lty=4)
legend("bottomright",c(paste("AUC1=",round(LD1_roc_4_1$auc,3),sep = ""),
                       paste("AUC2=",round(LD1_roc_4_2$auc,3),sep = ""),
                       paste("AUC3=",round(LD1_roc_4_3$auc,3),sep = "")),
       lty=c(1,3,4),col = c("purple","palegreen","red"),lwd=4,bty = "n",cex = 2)
dev.off()

plot(roc,print.auc=FALSE, print.auc.x=0.2,print.auc.y=0.2,cex.lab=2,cex.axis=2,cex.main=3,
     font.lab=2,#font=2,
     max.auc.polygon=T, grid=c(0.5, 0.5),grid.col=c("green", "red"),max.auc.polygon.col="white",
     main="",xlab="False Positive Rate",ylab="True Positive Rate",col="black",legacy.axes=TRUE,lwd=8)
legend("bottomright",c(paste("AUC=",round(roc$auc,3),sep = ""),
                       paste("Sensitivity=",round(best$sensitivity,3),sep = ""),
                       paste("Specificity=",round(best$specificity,3),sep = ""),
                       paste("Accuracy=",round(best$accuracy,3),sep = "")),col = c("red"),
       #lwd = 1,
       bty = "n",cex = 2)




#ROC
AUC2=data.frame(c(round(LD1_roc_4_1$auc,3),round(LD1_roc_4_2$auc,3),round(LD1_roc_4_3$auc,3)))
colnames(AUC2)[1]=c("AUC")
Table_LD1_2=rbind(LD1_best_4_1[,(1:8)],LD1_best_4_2[,(1:8)],LD1_best_4_3[,(1:8)])
Table_LD1_2$AUC=AUC2$AUC
rownames(Table_LD1_2)=c("S_vs_C","SS_vs_C","SS_vs_S")
Table_LD1_2$No=c(1:3)
#density plot of LD1 values
aLD1value_4_1=LD1value_4_1
aLD1value_4_2=LD1value_4_2
aLD1value_4_3=LD1value_4_3

LD1value_4_1$Group[LD1value_4_1$Group=="Patient"]=c("Sepsis")
LD1value_4_2$Group[LD1value_4_2$Group=="Patient"]=c("Septic shock")
LD1value_4_3$Group[LD1value_4_3$Group=="Patient"]=c("Septic shock")
LD1value_4_3$Group[LD1value_4_3$Group=="Control"]=c("Sepsis")
List8=list(LD1value_4_1[,-3],LD1value_4_2[,-3],LD1value_4_3[,-3])
#Control=c("Control","Control","Sepsis")
#Patient=c("Sepsis","Septic shock","Septic shock")
#Dataset=c("Control vs Sepsis         ","Control vs Septic shock         ","Sepsis vs Septic shock         ")
#N
Table_LD1_pred2=Table_LD1_2[,(5:10)]
colnames(Table_LD1_pred2)[5]=c("Accuracy1")
Table_LD1_pred2$tn2=Table_LD1_pred2$tn
Table_LD1_pred2$tp2=Table_LD1_pred2$tp
Table_LD1_pred2$fn2=Table_LD1_pred2$fn
Table_LD1_pred2$fp2=Table_LD1_pred2$fp
Table_LD1_pred2$Accuracy2=Table_LD1_pred2$Accuracy1
for (i in N) {
  data=data.frame(List8[i])
  threshold1=-1.404
  data$pred1[data$LD1<threshold1]=c("Control")
  data$pred1[data$LD1>threshold1]=c("Patient")
  result1=data.frame(table(data$Group, data$pred1))
  Table_LD1_pred2[i,1]=result1[1,3]
  Table_LD1_pred2[i,2]=result1[4,3]
  Table_LD1_pred2[i,3]=result1[2,3]
  Table_LD1_pred2[i,4]=result1[3,3]
  Table_LD1_pred2[i,5]=(result1[1,3]+result1[4,3])/(result1[1,3]+result1[2,3]+result1[3,3]+result1[4,3])
  threshold2=-0.639
  data$pred2[data$LD1<threshold2]=c("Control")
  data$pred2[data$LD1>threshold2]=c("Patient")
  result2=data.frame(table(data$Group, data$pred2))
  Table_LD1_pred2[i,7]=result2[1,3]
  Table_LD1_pred2[i,8]=result2[4,3]
  Table_LD1_pred2[i,9]=result2[2,3]
  Table_LD1_pred2[i,10]=result2[3,3]
  Table_LD1_pred2[i,11]=(result2[1,3]+result2[4,3])/(result2[1,3]+result2[2,3]+result2[3,3]+result2[4,3])
}

N
for (i in N) {
  data=data.frame(List8[i])
  yourfilename=paste("LD1_density_4_",i,".tiff",sep="")
  tiff(file=yourfilename,width = 1600,height = 1600,res = 300)
  p=ggplot(data,aes(x=LD1,fill=Group))+geom_histogram(aes(y=..density..),binwidth = .2,alpha=.5)+
    geom_density(alpha=.2)+geom_boxplot(width=.2, position=position_nudge(x=0,y=1.0))+
    scale_x_continuous(limits=c(-6, 6), breaks=c(-6,-4,-2,0,2,4,6))+
    geom_vline(xintercept = -1.404,lty=2,lwd=1.5)+
    geom_vline(xintercept = -0.639,lty=3,lwd=1.5)+
    scale_fill_manual(values = c("turquoise4","tomato"))+
    scale_color_manual(values = c("turquoise4","tomato"))+
    scale_y_continuous(limits=c(0, 3), breaks=c(0,1,2))+theme_rbook()
  print(p)
  dev.off()
}


#Risk score model
List9=list(data.frame(t(data4_H_S)),data.frame(t(data4_H_SS)),data.frame(t(data4_S_SS)))
List10=list(nrDEG4_1,nrDEG4_2,nrDEG4_3)
#List6
N=c(1:3)
#N0
for (n in N) {
  risk=data.frame(List9[n])
  Sum=data.frame(List6[n])
  rownames(Sum)=Sum$id
  for(i in Gene$Gene){
    GAPDH=risk[,colnames(risk)==c("GAPDH")]
    GAPDH=data.frame(GAPDH)
    rownames(GAPDH)=rownames(risk)
    data=data.frame(risk[,colnames(risk)==c(i)])
    rownames(data)=rownames(risk)
    colnames(data)[1]=c(i)
    data$relative=data[,1]/GAPDH$GAPDH
    Mean=mean(data$relative)
    Std=sd(data$relative)
    Weight=ifelse(data.frame(List10[n])[rownames(data.frame(List10[n]))==c(i),]$log2FoldChange<0,-1,1)
    data$riskscore=(data$relative-Mean)/Std*Weight
    assign(paste0("data",i),data)
    for (k in N0) {
      Sum[,N0[k]]=data[,3]
    }
    colnames(Sum)=Gene$Gene
    assign(paste0("Riskscore_Sum_4_",n),Sum)
  }
  RScore=Sum[,1]+Sum[,2]+Sum[,3]+Sum[,4]+Sum[,5]+Sum[,6]
  RScore=data.frame(RScore)
  RScore$id=rownames(risk)
  Group=data.frame(List6[n])
  RScore=merge(RScore,Group,by='id')
  roc=roc(RScore$Group,RScore$RScore)
  best=coords(roc,"best",ret = "all",transpose = F)
  assign(paste0("RScore_roc_4_",n),roc)
  assign(paste0("RScore_best_4_",n),best)
  assign(paste0("RScore_4_",n),RScore)
  yourfilename=paste("ROCRisk_4_",n,".tiff",sep="")
  tiff(file=yourfilename,width = 1600,height = 1600,res = 300)
  plot(roc,print.auc=FALSE, print.auc.x=0.2,print.auc.y=0.2,cex.lab=2,cex.axis=2,cex.main=3,
       font.lab=2,#font=2,
       max.auc.polygon=T, grid=c(0.5, 0.5),grid.col=c("green", "red"),max.auc.polygon.col="white",
       main="",xlab="False Positive Rate",ylab="True Positive Rate",col="black",legacy.axes=TRUE,lwd=8)
  
  dev.off()
}
AUC2=data.frame(c(round(RScore_roc_4_1$auc,3),round(RScore_roc_4_2$auc,3),round(RScore_roc_4_3$auc,3)))
colnames(AUC2)[1]=c("AUC")
Table_RScore2=rbind(RScore_best_4_1[,(1:8)],RScore_best_4_2[,(1:8)],RScore_best_4_3[,(1:8)])
Table_RScore2$AUC=AUC2$AUC
rownames(Table_RScore2)=c("S_vs_C","SS_vs_C","SS_vs_S")

List10_1=list(RScore_4_1,RScore_4_2,RScore_4_3)
Control=c("Control","Control","Sepsis")
Patient=c("Sepsis","Septic shock","Septic shock")
Dataset=c("Control vs Sepsis         ","Control vs Septic shock         ","Sepsis vs Septic shock         ")
Table_RScore_pred2=Table_RScore2[,(5:9)]
colnames(Table_RScore_pred2)[5]=c("Accuracy1")
Table_RScore_pred2$tn2=Table_RScore_pred2$tn
Table_RScore_pred2$tp2=Table_RScore_pred2$tp
Table_RScore_pred2$fn2=Table_RScore_pred2$fn
Table_RScore_pred2$fp2=Table_RScore_pred2$fp
Table_RScore_pred2$Accuracy2=Table_RScore_pred2$Accuracy1
for (i in N) {
  data=data.frame(List10_1[i])
  cutoff1=-2.035
  data$pred1[data$RScore<cutoff1]=c("Control")
  data$pred1[data$RScore>cutoff1]=c("Patient")
  result1=data.frame(table(data$Group, data$pred1))
  Table_RScore_pred2[i,1]=result1[1,3]
  Table_RScore_pred2[i,2]=result1[4,3]
  Table_RScore_pred2[i,3]=result1[2,3]
  Table_RScore_pred2[i,4]=result1[3,3]
  Table_RScore_pred2[i,5]=(result1[1,3]+result1[4,3])/(result1[1,3]+result1[2,3]+result1[3,3]+result1[4,3])
  cutoff2=0.703
  data$pred2[data$RScore<cutoff2]=c("Control")
  data$pred2[data$RScore>cutoff2]=c("Patient")
  result2=data.frame(table(data$Group, data$pred2))
  Table_RScore_pred2[i,6]=result2[1,3]
  Table_RScore_pred2[i,7]=result2[4,3]
  Table_RScore_pred2[i,8]=result2[2,3]
  Table_RScore_pred2[i,9]=result2[3,3]
  Table_RScore_pred2[i,10]=(result2[1,3]+result2[4,3])/(result2[1,3]+result2[2,3]+result2[3,3]+result2[4,3])
}


#density plot of risk scores
aRScore_4_1=RScore_4_1
aRScore_4_2=RScore_4_2
aRScore_4_3=RScore_4_3
RScore_4_1$Group[RScore_4_1$Group=="Patient"]=c("Sepsis")
RScore_4_2$Group[RScore_4_2$Group=="Patient"]=c("Septic shock")
RScore_4_3$Group[RScore_4_3$Group=="Control"]=c("Sepsis")
RScore_4_3$Group[RScore_4_3$Group=="Patient"]=c("Septic shock")
RScore_4_3=RScore_4_3[order(RScore_4_3$Group,decreasing = F),]
List11=list(RScore_4_1[,-1],RScore_4_2[,-1],RScore_4_3[,-1])
#N
for (i in N) {
  data=data.frame(List11[i])
  yourfilename=paste("AFig4_Density_6genes_RS_4-",Dataset[i],".tiff",sep="")
  tiff(file=yourfilename,width = 1600,height = 1600,res = 300)
  print(ggplot(data)+
          geom_density(aes(x=RScore,fill=Group),alpha=0.5)+
          scale_x_continuous(limits=c(-20, 20), breaks=c(-20,-10,0,10,20))+
          scale_y_continuous(limits=c(0, 0.3), breaks=c(0,0.1,0.2,0.3))+
          labs(y="density",x="Risk score",fill="Group"
               #,title = Dataset[i]
          )+
          scale_fill_manual(
            values = c("turquoise4","tomato"))+
          geom_vline(xintercept = -2.035,lty=2,lwd=1.5)+
          geom_vline(xintercept = 0.703,lty=3,lwd=1.5)+
          #theme(plot.title = element_text(size=33,face="bold",hjust = 0.5)) +
          theme(legend.position="none")+
          theme_rbook()+
          guides(fill=FALSE))
  dev.off()
}
NEWNEW

#density plot of risk scores
RScore_4_1$Group[RScore_4_1$Group=="Sepsis"]=c("Sep")
RScore_4_2$Group[RScore_4_2$Group=="Septic shock"]=c("SS")
RScore_4_3$Group[RScore_4_3$Group=="Septic shock"]=c("SS")
RScore_4_1$Group[RScore_4_1$Group=="Control"]=c("HC")
RScore_4_2$Group[RScore_4_2$Group=="Control"]=c("HC")
RScore_4_3$Group[RScore_4_3$Group=="Sepsis"]=c("Sep")
List5_4=list(RScore_4_1[,-1],RScore_4_2[,-1],RScore_4_3[,-1])
#N
N=1
for (i in N) {
  data=data.frame(List5_4[i])
  yourfilename=paste("AFig4_Density_6genes_RS_4-",i,".tiff",sep="")
  tiff(file=yourfilename,width = 2000,height = 1600,res = 300)
  print(ggplot(data)+
          geom_density(aes(x=RScore,fill=Group),alpha=0.5)+
          scale_x_continuous(limits=c(-20, 20), breaks=c(-20,-10,0,10,20))+
          scale_y_continuous(limits=c(0, 0.2), breaks=c(0,0.05,0.1,0.15,0.2))+
          labs(y="density",x="Risk score",fill="Group"
               #,title = Dataset[i]
          )+
          scale_fill_manual(
            values = c("gray80","green4"))+
          #geom_vline(xintercept = -2.035,lty=2,lwd=1.5)+
          #geom_vline(xintercept = 0.703,lty=3,lwd=1.5)+
          #theme(plot.title = element_text(size=33,face="bold",hjust = 0.5)) +
          theme(legend.position="none")+
          theme_rbook()+
          guides(fill=FALSE))
  dev.off()
}
N=2
for (i in N) {
  data=data.frame(List5_4[i])
  yourfilename=paste("Fig4_Density_6genes_RS_4-",i,".tiff",sep="")
  tiff(file=yourfilename,width = 2000,height = 1600,res = 300)
  print(ggplot(data)+
          geom_density(aes(x=RScore,fill=Group),alpha=0.5)+
          scale_x_continuous(limits=c(-20, 20), breaks=c(-20,-10,0,10,20))+
          scale_y_continuous(limits=c(0, 0.2), breaks=c(0,0.05,0.1,0.15,0.2))+
          labs(y="density",x="Risk score",fill="Group"
               #,title = Dataset[i]
          )+
          scale_fill_manual(
            values = c("gray80","dodgerblue3"))+
          #geom_vline(xintercept = -2.035,lty=2,lwd=1.5)+
          #geom_vline(xintercept = 0.703,lty=3,lwd=1.5)+
          #theme(plot.title = element_text(size=33,face="bold",hjust = 0.5)) +
          theme(legend.position="none")+
          theme_rbook()+
          guides(fill=FALSE))
  dev.off()
}
N=3
for (i in N) {
  data=data.frame(List5_4[i])
  yourfilename=paste("Fig4_Density_6genes_RS_4-",i,".tiff",sep="")
  tiff(file=yourfilename,width = 2000,height = 1600,res = 300)
  print(ggplot(data)+
          geom_density(aes(x=RScore,fill=Group),alpha=0.5)+
          scale_x_continuous(limits=c(-20, 20), breaks=c(-20,-10,0,10,20))+
          scale_y_continuous(limits=c(0, 0.2), breaks=c(0,0.05,0.1,0.15,0.2))+
          labs(y="density",x="Risk score",fill="Group"
               #,title = Dataset[i]
          )+
          scale_fill_manual(
            values = c("green4","dodgerblue3"))+
          #geom_vline(xintercept = -2.035,lty=2,lwd=1.5)+
          #geom_vline(xintercept = 0.703,lty=3,lwd=1.5)+
          #theme(plot.title = element_text(size=33,face="bold",hjust = 0.5)) +
          theme(legend.position="none")+
          theme_rbook()+
          guides(fill=FALSE))
  dev.off()
}


#Table_6genes accuracy AUC_1st part
AUC=data.frame(c(round(LD1_roc_1$auc,3),round(LD1_roc_2$auc,3),round(LD1_roc_3$auc,3),round(LD1_roc_4_1$auc,3),round(LD1_roc_4_2$auc,3),
                 round(LD1_roc_4_3$auc,3)))
colnames(AUC)[1]=c("AUC")
Table_LD1=rbind(LD1_best_1[,(1:8)],LD1_best_2[,(1:8)],LD1_best_3[,(1:8)],LD1_best_4_1[,(1:8)],LD1_best_4_2[,(1:8)],
                LD1_best_4_3[,(1:8)])
Table_LD1$AUC=AUC$AUC
rownames(Table_LD1)=c("GSE95233","GSE57065","GSE54514","GSE154918Sepsis","GSE154918Septic","GSE154918SS")
Table_LD1$No=c(1:6)


#Table_6genes accuracy AUC
AUC2=data.frame(c(round(RScore_roc_1$auc,3),round(RScore_roc_2$auc,3),round(RScore_roc_3$auc,3),round(RScore_roc_4_1$auc,3),round(RScore_roc_4_2$auc,3),
                  round(RScore_roc_4_3$auc,3)))
colnames(AUC2)[1]=c("AUC")
Table_RScore=rbind(RScore_best_1[,(1:8)],RScore_best_2[,(1:8)],RScore_best_3[,(1:8)],RScore_best_4_1[,(1:8)],RScore_best_4_2[,(1:8)],
                   RScore_best_4_3[,(1:8)])
Table_RScore$AUC=AUC2$AUC
Table_RScore=Table_RScore[-11,]
rownames(Table_RScore)=c("GSE95233","GSE57065","GSE54514","GSE154918Sepsis","GSE154918Septic","GSE154918SS")
Table_RScore$No=c(1:6)


#Multiple nonlinear regression model#6genes
#varibles:1st LD1 value;2nd Risk score;3rd: group information(paitent or control)
List18=list(aLD1value_1,aLD1value_2,aLD1value_3,aLD1value_4_1,aLD1value_4_2,aLD1value_4_3)
List19=list(aRScore_1,aRScore_2,aRScore_3,aRScore_4_1,aRScore_4_2,aRScore_4_3)
N4=c(1:6)
for (i in N4) {
  lm=merge(data.frame(List18[i]),data.frame(List19[i]),by='id')
  lm$Group[lm$Group.x=="Control"]=0
  lm$Group[lm$Group.x=="Patient"]=1
  rownames(lm)=lm$id
  lm=lm[,-5]
  lm=lm[,-3]
  lm=lm[,-1]
  assign(paste0("lm",i),lm)
}
#Lasso 
#based on combined data: GSE95233 GSE57065
lm12=rbind(lm1,lm2)
nsample=nrow(lm12)
ngenes=ncol(lm12[,!colnames(lm12)%in%c("Group")])
data.df=lm12
# split data into training and test sets
set.seed(1111)
train.ix1 = sample(1 : nsample, round(nsample*0.7))
test.ix1 = setdiff(1 : nsample, train.ix1)
train.df1 = data.df[train.ix1, ]
test.df1 = data.df[test.ix1, ]
# check class proportions in training and test sets
table(train.df1$Group); table(test.df1$Group)
set.seed(1234)
x1=as.matrix(train.df1[,(1:2)])
outcome1=train.df1$Group
#fit1=cv.glmnet(x1,outcome1,family="binomial",type.measure = "deviance")
fit1=cv.glmnet(x1,outcome1,family="binomial",type.measure = "deviance",alpha=1,nfolds = 20)
plot(fit1,cex.axis=1.5,cex.lab=1.5,cex.main=2)#688*416
coef(fit1,s=fit1$lambda.min)
#coef(fit1,s=fit1$lambda.1se)
#predict other data

List20=list(train.df1,test.df1,lm12,lm1,lm2,lm3,lm4,lm5,lm6)
#N1
Table_lm=Table_RScore
rownames(Table_lm)=c("training","testing","GSE95233GSE57065","GSE95233","GSE57065","GSE54514","GSE154918Sepsis","GSE154918Septic","GSE154918SS")
Table_lm$No=c(1:9)
N1=c(1:9)
for (i in N1) {
  outcome=data.frame(List20[i])$Group
  pred=predict(fit1, newx=as.matrix(data.frame(List20[i])[,1:2]) ,s=c(fit1$lambda.min,fit1$lambda.1se))#
  re=as.data.frame(cbind(data.frame(List20[i])$Group,pred))
  colnames(re)=c("Group","prob_min","prob_1se")
  reL=re
  reL$Group[reL$Group==0]=c("Control")
  reL$Group[reL$Group==1]=c("Septic shock")
  reL=reL[order(reL$Group,decreasing = F),]
  re$Group=as.factor(re$Group)
  #ROC
  roc=roc(re$Group,re$prob_min)
  best=coords(roc,"best",ret = "all",transpose = F)
  Table_lm[i,(1:8)]=best[1,(1:8)]
  assign(paste0("lm_roc_",i),roc)
  assign(paste0("lm_best_",i),best)
  #assign(paste0("accuracy_LR_",i),best$accuracy)
  pred_min=prediction(re[,2],re[,1])
  auc_min=performance(pred_min,"auc")@y.values[[1]]
  perf_min=performance(pred_min,"tpr","fpr")
  assign(paste0("re",i),re)
  assign(paste0("reL",i),reL)
  assign(paste0("perf_min",i),perf_min)
  assign(paste0("pred_min",i),pred_min)
  Table_lm[i,9]=auc_min
  #assign(paste0("auc_min",i),auc_min)
}

List21=list(aLD1value_1,aLD1value_2,aLD1value_3,aLD1value_4_1,aLD1value_4_2,
            aLD1value_4_3)
List22=list(aRScore_1,aRScore_2,aRScore_3,aRScore_4_1,aRScore_4_2,aRScore_4_3)
List23=list(reL4,reL5,reL6,reL7,reL8,reL9)

N4
table1=table
rownames(table1)=c("GSE95233","GSE57065","GSE54514","GSE154918Sepsis","GSE154918Septic","GSE154918SS")

table2=table
rownames(table2)=c("GSE95233","GSE57065","GSE54514","GSE154918Sepsis","GSE154918Septic","GSE154918SS")

for (i in N4) {
  LDA=data.frame(List21[i])
  THR1=-1.404
  THR2=-0.0639
  LDA$Group1[LDA$LD1<THR1]=c("Control")
  LDA$Group1[LDA$LD1>THR1]=c("Patient")
  LDA$Group2[LDA$LD1<THR2]=c("Control")
  LDA$Group2[LDA$LD1>THR2]=c("Patient")
  result_LDA1=data.frame(table(LDA$Group, LDA$Group1))
  result_LDA2=data.frame(table(LDA$Group, LDA$Group2))
  table1[i,1]=result_LDA1[1,3]
  table1[i,2]=result_LDA1[2,3]
  table1[i,3]=result_LDA1[3,3]
  table1[i,4]=result_LDA1[4,3]
  table1[i,5]=(result_LDA1[1,3]+result_LDA1[4,3])/(result_LDA1[1,3]+result_LDA1[2,3]+result_LDA1[3,3]+result_LDA1[4,3])
  table1[i+14,1]=result_LDA2[1,3]
  table1[i+14,2]=result_LDA2[2,3]
  table1[i+14,3]=result_LDA2[3,3]
  table1[i+14,4]=result_LDA2[4,3]
  table1[i+14,5]=(result_LDA2[1,3]+result_LDA2[4,3])/(result_LDA2[1,3]+result_LDA2[2,3]+result_LDA2[3,3]+result_LDA2[4,3])
  
  RScore=data.frame(List22[i])
  #cutoff1,cutoff2
  RScore$Group1[RScore$RScore<cutoff1]=c("Control")
  RScore$Group1[RScore$RScore>cutoff1]=c("Patient")
  RScore$Group2[RScore$RScore<cutoff2]=c("Control")
  RScore$Group2[RScore$RScore>cutoff2]=c("Patient")
  result_RScore1=data.frame(table(RScore$Group, RScore$Group1))
  result_RScore2=data.frame(table(RScore$Group, RScore$Group2))
  table2[i,1]=result_RScore1[1,3]
  table2[i,2]=result_RScore1[2,3]
  table2[i,3]=result_RScore1[3,3]
  table2[i,4]=result_RScore1[4,3]
  table2[i,5]=(result_RScore1[1,3]+result_RScore1[4,3])/(result_RScore1[1,3]+result_RScore1[2,3]+result_RScore1[3,3]+result_RScore1[4,3])
  table2[i+14,1]=result_RScore2[1,3]
  table2[i+14,2]=result_RScore2[2,3]
  table2[i+14,3]=result_RScore2[3,3]
  table2[i+14,4]=result_RScore2[4,3]
  table2[i+14,5]=(result_RScore2[1,3]+result_RScore2[4,3])/(result_RScore2[1,3]+result_RScore2[2,3]+result_RScore2[3,3]+result_RScore2[4,3])
  
}
write.csv(table1,"table1predictforfig5.csv")
write.csv(table2,"table2predictforfig5.csv")
write.csv(Table_lm,"tablelmpredictforfig5.csv")

#"gray80","dodgerblue3" green4

#Volin plot Septic shock vs Control
List15=list(reL2,reL3,reL4,reL5,reL8)
title1=c("testing data","GSE95233+GSE57065","GSE95233","GSE57065","GSE154918")
N2=c(1:5)
for (i in N2) {
  data=data.frame(List15[i])
  yourfilename=paste("Fig6_LR_1st_6genes_septicshock-",title1[i],".tiff",sep="")
  tiff(file=yourfilename,width = 1500,height = 1600,res = 300)
  print(ggviolin(data,x="Group",y="prob_min",fill = "Group",palette = c("gray80","dodgerblue3"),
                 add = "boxplot",add.params = list(fill="white"),alpha=0.5)+guides(fill="none")+
          scale_y_continuous(limits=c(-40, 60), breaks=c(-40,0,40))+
          geom_signif(comparisons = list(c("Control", "Septic shock")),y_position = c(45, 48),
                      map_signif_level=T,textsize=6,test=wilcox.test,step_increase=0.2,vjust = 0,hjust=0.5)+
          labs(y="Probability",x="",fill="Group",title = title1[i])+
          theme(plot.title = element_text(size=24,hjust = 0.5))+#face="bold",
          theme_rbook())
  dev.off()
}
#Volin plot Sepsis vs Control
reL6$Group[reL6$Group=="Septic shock"]=c("Sepsis")
reL7$Group[reL7$Group=="Septic shock"]=c("Sepsis")
reL15$Group[reL15$Group=="Septic shock"]=c("Sepsis")
List16=list(reL6,reL7)
title2=c("GSE54514","GSE154918")
#N
for (i in N) {
  data=data.frame(List16[i])
  yourfilename=paste("Fig6_LR_1st_6genes_sepsis-",title2[i],".tiff",sep="")
  tiff(file=yourfilename,width = 1500,height = 1600,res = 300)
  print(ggviolin(data,x="Group",y="prob_min",fill = "Group",palette = c("gray80","green4"),
                 add = "boxplot",add.params = list(fill="white"),alpha=0.5)+guides(fill="none")+
          scale_y_continuous(limits=c(-40, 60), breaks=c(-40,0,40))+
          geom_signif(comparisons = list(c("Control", "Sepsis")),y_position = c(45, 48),
                      map_signif_level=T,textsize=6,test=wilcox.test,step_increase=0.2,vjust = 0,hjust=0.5)+
          labs(y="Probability",x="",fill="Group",title = title2[i])+
          theme(plot.title = element_text(size=24))+#,just = 0.5
          theme_rbook())
  dev.off()
}
#Volin plot Septic shock vs Sepsis
reL9$Group[reL9$Group=="Control"]=c("Sepsis")
reL17$Group[reL17$Group=="Control"]=c("Sepsis")
List17=list(reL9,reL17)
title3=c("GSE154918")
N3=c(1)
for (i in N3) {
  data=data.frame(List17[i])
  yourfilename=paste("Fig6_LR_1st_6genes_sepsis_septicshock-",title3[i],".tiff",sep="")
  tiff(file=yourfilename,width = 1500,height = 1600,res = 300)
  print(ggviolin(data,x="Group",y="prob_min",fill = "Group",palette = c("green4","dodgerblue3"),
                 add = "boxplot",add.params = list(fill="white"),alpha=0.5)+guides(fill="none")+
          scale_y_continuous(limits=c(-40, 60), breaks=c(-40,0,40))+
          geom_signif(comparisons = list(c("Sepsis", "Septic shock")),y_position = c(45, 48),
                      map_signif_level=T,textsize=6,test=wilcox.test,step_increase=0.2,vjust = 0,hjust=0.5)+
          labs(y="Probability",x="",fill="Group",title = title3[i])+
          theme(plot.title = element_text(size=24,hjust = 0.5))+
          theme_rbook())
  dev.off()
}




#Immune cells
CIBERSORTx_1=CIBERSORTx_GSE95233_Results
CIBERSORTx_1=CIBERSORTx_1[CIBERSORTx_1$Mixture%in%pData1$X,]
CIBERSORTx_2=CIBERSORTx_GSE57065_Results
CIBERSORTx_2=CIBERSORTx_2[CIBERSORTx_2$Mixture%in%pData2$X,]
CIBERSORTx_3=CIBERSORTx_GSE54514_Results
CIBERSORTx_3=CIBERSORTx_3[CIBERSORTx_3$Mixture%in%pData3$X,]
List14=list(CIBERSORTx_1,CIBERSORTx_2,CIBERSORTx_3)
group_list
for (i in N) {
  data=data.frame(List14[i])
  data$NK=data$NK.cells.activated+data$NK.cells.resting
  data2=data[,-(1:26)]
  data2=data.frame(data2)
  colnames(data2)=c("NK")
  data2$id=data$Mixture
  data2=merge(data2,data.frame(group_list[i]),by='id')
  assign(paste0("CIBERSORT_NK",i),data2)
}
#Immune cells
#CIBERSORTx_Job43_ResultsGSE154918
CIBERSORTx_Job43_ResultsGSE154918$NK=CIBERSORTx_Job43_ResultsGSE154918$NK.cells.activated+CIBERSORTx_Job43_ResultsGSE154918$NK.cells.resting
#GSE95233
CIBERSORT_NK_adult_1=CIBERSORT_NK1
CIBERSORT_NK_adult_1$Group[CIBERSORT_NK_adult_1$Group=="Patient"]=c("Septic_shock")
CIBERSORT_NK_adult_1$Dataset[CIBERSORT_NK_adult_1$Group=="Control"|CIBERSORT_NK_adult_1$Group=="Septic_shock"]=c("GSE95233")
CIBERSORT_NK_adult_1=CIBERSORT_NK_adult_1[,-1]
#GSE57065
CIBERSORT_NK_adult_2=CIBERSORT_NK2
CIBERSORT_NK_adult_2$Group[CIBERSORT_NK_adult_2$Group=="Patient"]=c("Septic_shock")
CIBERSORT_NK_adult_2$Dataset[CIBERSORT_NK_adult_2$Group=="Control"|CIBERSORT_NK_adult_2$Group=="Septic_shock"]=c("GSE57065")
CIBERSORT_NK_adult_2=CIBERSORT_NK_adult_2[,-1]
#GSE54514
CIBERSORT_NK_adult_3=CIBERSORT_NK3
CIBERSORT_NK_adult_3$Group[CIBERSORT_NK_adult_3$Group=="Patient"]=c("Sepsis")
CIBERSORT_NK_adult_3$Dataset[CIBERSORT_NK_adult_3$Group=="Control"|CIBERSORT_NK_adult_3$Group=="Sepsis"]=c("GSE54514")
CIBERSORT_NK_adult_3=CIBERSORT_NK_adult_3[,-1]
#GSE154918
CIBERSORT_NK_adult_4=data.frame(CIBERSORTx_Job43_ResultsGSE154918$NK)
CIBERSORT_NK_adult_4$id=CIBERSORTx_Job43_ResultsGSE154918$Mixture
CIBERSORT_NK_adult_4=merge(CIBERSORT_NK_adult_4,pData4,by='id')
write.csv(CIBERSORT_NK1,"GSE95233_NK.csv")
write.csv(CIBERSORT_NK2,"GSE57065_NK.csv")
write.csv(CIBERSORT_NK3,"GSE54514_NK.csv")
write.csv(CIBERSORT_NK_adult_4,"GSE154918_NK.csv")

CIBERSORT_NK_adult_4$Dataset[CIBERSORT_NK_adult_4$Group=="Hlty"|CIBERSORT_NK_adult_4$Group=="Shock_P"|CIBERSORT_NK_adult_4$Group=="Seps_P"]=c("GSE154918")
CIBERSORT_NK_adult_4$Group[CIBERSORT_NK_adult_4$Group=="Hlty"]=c("Control")
CIBERSORT_NK_adult_4$Group[CIBERSORT_NK_adult_4$Group=="Shock_P"]=c("Septic_shock")
CIBERSORT_NK_adult_4$Group[CIBERSORT_NK_adult_4$Group=="Seps_P"]=c("Sepsis")
CIBERSORT_NK_adult_4=CIBERSORT_NK_adult_4[,-1]
colnames(CIBERSORT_NK_adult_4)[1]=c("NK")
CIBERSORT_NK_adult_4_1=CIBERSORT_NK_adult_4[which(CIBERSORT_NK_adult_4$Group=="Control"|CIBERSORT_NK_adult_4$Group=="Sepsis"),]
CIBERSORT_NK_adult_4_2=CIBERSORT_NK_adult_4[which(CIBERSORT_NK_adult_4$Group=="Control"|CIBERSORT_NK_adult_4$Group=="Septic_shock"),]
CIBERSORT_NK_adult_4_3=CIBERSORT_NK_adult_4[which(CIBERSORT_NK_adult_4$Group=="Sepsis"|CIBERSORT_NK_adult_4$Group=="Septic_shock"),]
CIBERSORT_NK_adult_septicshock=rbind(CIBERSORT_NK_adult_1,CIBERSORT_NK_adult_2,CIBERSORT_NK_adult_4_2)
CIBERSORT_NK_adult_sepsis=rbind(CIBERSORT_NK_adult_3,CIBERSORT_NK_adult_4_1)
#Violin
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}+scale_fill_manual(values = c("#7697CB","#D293C1"))
#"gray80","dodgerblue3" green4
#split_violin_septic_shock
CIBERSORT_NK_adult_septicshock$Dataset=factor(CIBERSORT_NK_adult_septicshock$Dataset,levels = c('GSE95233','GSE57065','GSE154918'))
p1=ggplot(CIBERSORT_NK_adult_septicshock, aes(Dataset, NK, fill = Group)) + geom_split_violin()+theme_rbook()+
  scale_fill_manual(labels=c("Control","Septic shock"),values = c("gray80","dodgerblue3"))+
  geom_boxplot(width = 0.16, notch = FALSE, notchwidth = .3, outlier.shape = NA,coef=0) +
  labs(x=NULL,y="Estimated NK cell proportion")+
  scale_y_continuous(limits=c(0, 0.3), breaks=c(0,0.1,0.2,0.3))
#split_violin_sepsis
CIBERSORT_NK_adult_sepsis$Dataset=factor(CIBERSORT_NK_adult_sepsis$Dataset,levels = c('GSE54514','GSE154918'))
p2=ggplot(CIBERSORT_NK_adult_sepsis, aes(Dataset, NK, fill = Group)) + geom_split_violin()+theme_rbook()+
  scale_fill_manual(labels=c("Control","Sepsis"),values = c("gray80","green4"))+
  geom_boxplot(width = 0.16, notch = FALSE, notchwidth = .3, outlier.shape = NA,coef=0) +
  labs(x=NULL,y="Estimated NK cell proportion")+
  scale_y_continuous(limits=c(0, 0.3), breaks=c(0,0.1,0.2,0.3))
#split_violin_sepsis vs septic_shock
p3=ggplot(CIBERSORT_NK_adult_4_3, aes(Dataset, NK, fill = Group)) + geom_split_violin()+theme_rbook()+
  scale_fill_manual(labels=c("Sepsis","Septic shock"),values = c("green4","dodgerblue3"))+
  geom_boxplot(width = 0.16, notch = FALSE, notchwidth = .3, outlier.shape = NA,coef=0) +
  labs(x=NULL,y="Estimated NK cell proportion")+
  scale_y_continuous(limits=c(0, 0.3), breaks=c(0,0.1,0.2,0.3))
layout="AAABBC"
p1+p2+p3+plot_layout(design = layout)


#Predict testing data and other datasets
 #GSE154918
CIBERSORT_NK_4=data.frame(CIBERSORTx_Job43_ResultsGSE154918$NK)
CIBERSORT_NK_4$id=CIBERSORTx_Job43_ResultsGSE154918$Mixture
CIBERSORT_NK_4=merge(CIBERSORT_NK_4,pData4,by='id')
colnames(CIBERSORT_NK_4)[2]=c("NK")
CIBERSORT_NK_4_1=CIBERSORT_NK_4[which(CIBERSORT_NK_4$Group=="Hlty"|CIBERSORT_NK_4$Group=="Seps_P"),]
CIBERSORT_NK_4_2=CIBERSORT_NK_4[which(CIBERSORT_NK_4$Group=="Hlty"|CIBERSORT_NK_4$Group=="Shock_P"),]
CIBERSORT_NK_4_3=CIBERSORT_NK_4[which(CIBERSORT_NK_4$Group=="Seps_P"|CIBERSORT_NK_4$Group=="Shock_P"),]
CIBERSORT_NK_4_1$Group[CIBERSORT_NK_4_1$Group=="Hlty"]=0
CIBERSORT_NK_4_1$Group[CIBERSORT_NK_4_1$Group=="Seps_P"]=1
CIBERSORT_NK_4_2$Group[CIBERSORT_NK_4_2$Group=="Hlty"]=0
CIBERSORT_NK_4_2$Group[CIBERSORT_NK_4_2$Group=="Shock_P"]=1
CIBERSORT_NK_4_3$Group[CIBERSORT_NK_4_3$Group=="Seps_P"]=0
CIBERSORT_NK_4_3$Group[CIBERSORT_NK_4_3$Group=="Shock_P"]=1



save(lm_2nd_1,lm_2nd_2,lm_2nd_3,lm_2nd_4,lm_2nd_5,lm_2nd_6,file = "SepsisSeptic_1stdata.Rdata")
save(re3,re4,re5,re6,re7,re8,re9,re10,re11,
     re12,re13,re14,re15,re16,re17,file = "SepsisSeptic_1stdata_pred_min.Rdata")

load("Lass_2nd_3x_ROCs.Rdata")
roc.test(Lass_1st_3x_ROC_9,Lass_2nd_3x_ROC_9,method="delong")
roc.test(Lass_1st_3x_ROC_17,Lass_2nd_3x_ROC_17,method="delong")


#PCR human
outcome=lmPCR$Group
pred=predict(fit1, newx=as.matrix(lmPCR[,1:2]) , s=c(fit1$lambda.min,fit1$lambda.1se) )
re=as.data.frame(cbind(lmPCR$Group,pred))
colnames(re)=c("Group","prob_min","prob_1se")
reL=re
reL$Group[reL$Group==0]=c("Control")
reL$Group[reL$Group==1]=c("Septic shock")
reL=reL[order(reL$Group,decreasing = F),]
re$Group=as.factor(re$Group)
#ROC
roc=roc(re$Group,re$prob_min)
best=coords(roc,"best",ret = "all",transpose = F)

#assign(paste0("accuracy_LR_",i),best$accuracy)
pred_min=prediction(re[,2],re[,1])
auc_min=performance(pred_min,"auc")@y.values[[1]]#0.991
perf_min=performance(pred_min,"tpr","fpr")

data=reL
yourfilename=paste("HumanPCR_LR-",".tiff",sep="")
tiff(file=yourfilename,width = 1600,height = 1600,res = 300)
print(ggviolin(data,x="Group",y="prob_min",fill = "Group",palette = c("turquoise4","tomato"),
                 add = "boxplot",add.params = list(fill="white"),alpha=0.5)+guides(fill="none")+
          scale_y_continuous(limits=c(-40, 60), breaks=c(-40,0,40))+
          geom_signif(comparisons = list(c("Control", "Septic shock")),y_position = c(45, 48),
                      map_signif_level=T,textsize=6,test=wilcox.test,step_increase=0.2,vjust = 0,hjust=0.5)+
          labs(y="Probability",x="",fill="Group")+
          theme(plot.title = element_text(size=24,face="bold",hjust = 0.5))+
          theme_rbook())
  dev.off()

  #####human
  List30=list(lmPCR1,lmPCR2,lmPCR3)
  Table_lm_H=Table_RScore[(1:3),]
  rownames(Table_lm_H)=c("C_S","C_SS","S_SS")
  N
  for (i in N) {
    outcome=data.frame(List30[i])$Group
    pred=predict(fit1, newx=as.matrix(data.frame(List30[i])[,1:2]) , s=c(fit1$lambda.min,fit1$lambda.1se) )
    re=as.data.frame(cbind(data.frame(List30[i])$Group,pred))
    colnames(re)=c("Group","prob_min","prob_1se")
    reL=re
    reL$Group[reL$Group==0]=c("Control")
    reL$Group[reL$Group==1]=c("Septic shock")
    reL=reL[order(reL$Group,decreasing = F),]
    re$Group=as.factor(re$Group)
    #ROC
    roc=roc(re$Group,re$prob_min)
    best=coords(roc,"best",ret = "all",transpose = F)
    Table_lm_H[i,(1:8)]=best[1,(1:8)]
    assign(paste0("lm_H_roc_",i),roc)
    assign(paste0("lm_H_best_",i),best)
    #assign(paste0("accuracy_LR_",i),best$accuracy)
    pred_min=prediction(re[,2],re[,1])
    auc_min=performance(pred_min,"auc")@y.values[[1]]
    perf_min=performance(pred_min,"tpr","fpr")
    assign(paste0("re_H",i),re)
    assign(paste0("reL_H",i),reL)
    assign(paste0("perf_min_H",i),perf_min)
    assign(paste0("pred_min_H",i),pred_min)
    Table_lm_H[i,9]=auc_min
    #assign(paste0("auc_min",i),auc_min)
  }

write.csv(Table_lm,"Table_lm.csv")
write.csv(Table_lm_H,"Table_lm_H.csv")


#Lasso 
#based on combined data: GSE95233 GSE57065 GSE154918_sep lm5
lm125=rbind(lm1,lm2,lm5)
nsample=nrow(lm125)
ngenes=ncol(lm125[,!colnames(lm125)%in%c("Group")])
data.df=lm125
# split data into training and test sets
set.seed(1112)
train.ix1 = sample(1 : nsample, round(nsample*0.7))
test.ix1 = setdiff(1 : nsample, train.ix1)
train.df1 = data.df[train.ix1, ]
test.df1 = data.df[test.ix1, ]
# check class proportions in training and test sets
table(train.df1$Group); table(test.df1$Group)
set.seed(1234)
x1=as.matrix(train.df1[,(1:2)])
outcome1=train.df1$Group
#fit2=cv.glmnet(x1,outcome1,family="binomial",type.measure = "deviance")
fit2=cv.glmnet(x1,outcome1,family="binomial",type.measure = "deviance",alpha=1,nfolds = 20)
plot(fit2,cex.axis=1.5,cex.lab=1.5,cex.main=2)#688*416
coef(fit2,s=fit2$lambda.min)
#coef(fit1,s=fit1$lambda.1se)
#predict other data
save(fit1,fit2,file = "SepShockFindeR.Rdata")

List205=list(train.df1,test.df1,lm125,lm1,lm2,lm3,lm4,lm5,lm6,lmPCR1,lmPCR2,lmPCR3)
#N1
Table_lm2=Table_RScore
rownames(Table_lm2)=c("training","testing","Merged","GSE95233","GSE57065","GSE54514","GSE154918Sepsis","GSE154918Septic","GSE154918SS","human_sepsis","human_septic","humanSS")
Table_lm2$No=c(1:12)
N5=c(1:12)
for (i in N5) {
  outcome=data.frame(List205[i])$Group
  pred=predict(fit2, newx=as.matrix(data.frame(List205[i])[,1:2]) ,s=c(fit2$lambda.min,fit2$lambda.1se))#
  re=as.data.frame(cbind(data.frame(List205[i])$Group,pred))
  colnames(re)=c("Group","prob_min","prob_1se")
  reL=re
  reL$Group[reL$Group==0]=c("Control")
  reL$Group[reL$Group==1]=c("Septic shock")
  reL=reL[order(reL$Group,decreasing = F),]
  re$Group=as.factor(re$Group)
  #ROC
  roc=roc(re$Group,re$prob_min)
  best=coords(roc,"best",ret = "all",transpose = F)
  Table_lm2[i,(1:8)]=best[1,(1:8)]
  assign(paste0("lm_roc_",i),roc)
  assign(paste0("lm_best_",i),best)
  #assign(paste0("accuracy_LR_",i),best$accuracy)
  pred_min=prediction(re[,2],re[,1])
  auc_min=performance(pred_min,"auc")@y.values[[1]]
  perf_min=performance(pred_min,"tpr","fpr")
  assign(paste0("re",i),re)
  assign(paste0("reL",i),reL)
  assign(paste0("perf_min",i),perf_min)
  assign(paste0("pred_min",i),pred_min)
  Table_lm2[i,9]=auc_min
  #assign(paste0("auc_min",i),auc_min)
}

#Volin plot Septic shock vs Control
List15_2=list(reL2,reL3,reL8)
title1_2=c("testing data","Merged data","GSE154918")
N2_2=c(1:3)
for (i in N2_2) {
  data=data.frame(List15_2[i])
  yourfilename=paste("Fig5_SepxFindR_3merged_6genes_septicshock-",title1_2[i],".tiff",sep="")
  tiff(file=yourfilename,width = 1500,height = 1600,res = 300)
  print(ggviolin(data,x="Group",y="prob_min",fill = "Group",palette = c("gray80","dodgerblue3"),
                 add = "boxplot",add.params = list(fill="white"),alpha=0.5)+guides(fill="none")+
          scale_y_continuous(limits=c(-40, 60), breaks=c(-40,0,40))+
          geom_signif(comparisons = list(c("Control", "Septic shock")),y_position = c(45, 48),
                      map_signif_level=T,textsize=6,test=wilcox.test,step_increase=0.2,vjust = 0,hjust=0.5)+
          labs(y="Probability",x="",fill="Group",title = title1[i])+
          theme(plot.title = element_text(size=24,hjust = 0.5))+#face="bold",
          theme_rbook())
  dev.off()
}

#Volin plot Septic shock vs Sepsis
reL9$Group[reL9$Group=="Control"]=c("Sepsis")

title3=c("GSE154918")
N3=c(1)
for (i in N3) {
  data=reL9
  yourfilename=paste("Fig5_SepxFindR_3merged_6genes_sepsis_septicshock-",title3[i],".tiff",sep="")
  tiff(file=yourfilename,width = 1500,height = 1600,res = 300)
  print(ggviolin(data,x="Group",y="prob_min",fill = "Group",palette = c("green4","dodgerblue3"),
                 add = "boxplot",add.params = list(fill="white"),alpha=0.5)+guides(fill="none")+
          scale_y_continuous(limits=c(-40, 60), breaks=c(-40,0,40))+
          geom_signif(comparisons = list(c("Sepsis", "Septic shock")),y_position = c(45, 48),
                      map_signif_level=T,textsize=6,test=wilcox.test,step_increase=0.2,vjust = 0,hjust=0.5)+
          labs(y="Probability",x="",fill="Group",title = title3[i])+
          theme(plot.title = element_text(size=24,hjust = 0.5))+
          theme_rbook())
  dev.off()
}