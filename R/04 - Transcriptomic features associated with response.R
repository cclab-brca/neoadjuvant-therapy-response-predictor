# Manuscript:   Multi-omic machine learning predictor of breast cancer therapy response
# Author:       Stephen-John Sammut
# Description:  Transcriptomic features associated with response
# Section:      Results - Tumour proliferation and immune signatures
#=======================================================================================

rm (list=ls())

#load packages
library (data.table)
library (edgeR)
library (EnsDb.Hsapiens.v86)
library (fgsea)
library (ggbiplot)
library (ggplot2)
library (ggpmisc)
library (ggridges)
library (MASS)
library (org.Hs.eg.db)
library (Pigengene)
library (ReactomePA)
library (readxl)
library (sm)
library (stringr)
library (viridis)

# change the contents of variable baseDir to the root analysis folder 
baseDir <- "~/Neoadjuvant-predictor/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/R/LoadDirectoryStructure.R"))

figure_font_size=12

# load clinical metadata (Supplementary Table 1)
metadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 1))
# combine ER and HER2 status
metadata$ERHER2.status <- ifelse(metadata$ER.status=="POS","ER+ HER2-","ER- HER2-")
metadata$ERHER2.status <- ifelse(metadata$HER2.status=="POS","HER2+",metadata$ERHER2.status)
metadataFull <- metadata

# Perform analyses with cases that had an RCB assessment and received more than one cycle of therapy
# as detailed in Methods (Statistical testing)
metadata <- metadata[metadata$RCB.category!="NA",]
metadata <- metadata[metadata$Chemo.cycles>1 & metadata$aHER2.cycles>1,]

# load list of breast cancer driver genes from resources directory
driverGenes<- scan(paste0(resourcesDir,"breast-cancer-driver-genes.txt"), what=character(),skip = 1)

# load Gene Ensembl ID to Hugo ID dictionary
ensemblToHugo <- read.table(paste0(resourcesDir,"EnsemblID.to.Hugo.v87.tsv.gz"), header=T, stringsAsFactors = F,sep="\t")

# load RNA data (Supplementary Table 3)
rnadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 3))


#=========================================================================
# Differential gene expression and enrichment
# Figure 3a
#=========================================================================

p <- metadata

# load RNAseq raw counts (Methods)
transneo.counts <- data.frame(fread(paste0(dataDir,"transneo-diagnosis-RNAseq-rawcounts.tsv.gz"),header=T, sep="\t",stringsAsFactors = F),row.names = 1)
transneo.counts <- transneo.counts[,colnames(transneo.counts) %in% p$Donor.ID]

# update metadata - retain only samples that have RNAseq data
p <- p[p$Donor.ID %in% colnames(transneo.counts),]

# sanity check
stopifnot(sum(p$Donor.ID!=colnames(transneo.counts))==0)

# 149 tumours have RNAseq data, and associated RCB assessment + had adequate chemotherapy exposure
dim(p)

# construct differential expression matrix
y <- DGEList(transneo.counts)
minCPM      <- 1
minNoToKeep <- 10
keep        <- rowSums(cpm(y)>minCPM)>=minNoToKeep
y <- y[keep , , keep.lib.sizes=FALSE]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method = "TMM")

# construct design matrix
rcb    <- factor(p$pCR.RD, levels=c("RD","pCR"))
design <- model.matrix(~rcb, data = y$samples)
coef   <- ncol(design)
head(design,2)

# perform DE: pCR vs RD - as per edgeR documentation
y       <- estimateDisp(y, design = design, robust = TRUE)
fit     <- glmQLFit(y, design, robust=TRUE)
qlf     <- glmQLFTest(fit, coef=coef)
results <- as.data.frame(topTags(qlf,n = Inf))
is.de   <- decideTestsDGE(qlf, p.value=0.05)

# 2,071 genes are under-expressed, and 2,439 genes are over-expressed in tumours attaining pCR (FDR<0.05). 
summary(is.de)

# Annotate full gene list with their Hugo Gene name
annotatedResults <- results
annotatedResults[rownames(annotatedResults) %in% rownames(y)[is.de[,1]==1],"expression"] <- "overexpressed"
annotatedResults[rownames(annotatedResults) %in% rownames(y$counts)[is.de[,1]==-1],"expression"] <- "underexpressed"
annotatedResults[rownames(annotatedResults) %in% rownames(y$counts)[is.de[,1]==0],"expression"] <- "notDE"
annotatedResults <- merge(x=annotatedResults,y=ensemblToHugo, by.x=0,by.y=1, sort=F)

# Characterise expression landscape of driver genes and plot Figure 3a
driverExpression <- annotatedResults[annotatedResults$expression!="notDE",]
driverExpression <- driverExpression[driverExpression$Hugo %in% c( driverGenes),]
driverExpression <- driverExpression[abs(driverExpression$logFC)>0.5 & driverExpression$FDR<0.05,]
driverExpression <- driverExpression[,c(8,2,6,7)]
driverExpression[order(driverExpression$logFC,decreasing = T),]

fig3a <- 
  ggplot(driverExpression,aes(x=reorder(Hugo,logFC),y=(logFC),color=expression))+
  geom_point(aes(size= -log10(FDR)))+
  geom_hline(yintercept=0,color="black", linetype="dotted")+
  labs(x="",y="log FC pCR")+
  coord_flip()+
  scale_colour_manual(values = c("#FB6542","#375E97"))+
  scale_size_continuous(name=expression(italic(-log[10]~FDR)),breaks = c(2:4))+
  guides(color="none")+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face="italic"),
        plot.margin = unit(c(1.1,1.5,0,0.5), "lines"))
fig3a


#=========================================================================
# Differentially expressed pathways/gene sets in pCR vs RD
# Figure 3b, Extended Figure 5a
#=========================================================================

# a. MSigDB enrichment

# load MSigDB Hallmarks Gene set from resources folder and run enrichment
hallmarks <- paste0(resourcesDir,"MSigDB.HallmarksGeneSet.Rdata")
geneList  <- readRDS(hallmarks)

# edgeR camera
indexes    <- ids2indices(geneList, rownames(y$counts))
gst.camera <- camera.DGEList(y, index=indexes, design=design, coef=coef)
gst.camera <- gst.camera[gst.camera$FDR<0.01,]
# fgsea - get Normalised Enrichment Score
results.ord  <- results[ order(-results[,"logFC"]), ]
ranks        <- results.ord$logFC
names(ranks) <- rownames(results.ord)
g <- geneList[names(geneList) %in% rownames(gst.camera)]
q <- data.frame(fgsea(g, ranks, minSize=15, maxSize = 500, nperm=1000),stringsAsFactors = F)[,c(1:5)]
q <- q[q$padj<0.05,]

# create Figure 3b - MSigDB gene set enrichment
q$process <- gsub("HALLMARK_","",q$pathway)
q$process <- gsub("_"," ",q$process)
q$process <- str_to_sentence(q$process)
q$process <- gsub("E2f","E2F",q$process)
q$process <- gsub("Ifn","IFN",q$process)
q$process <- gsub("G2m","G2M",q$process)
q$process <- gsub("Il6 jak stat3","IL6 JAK STAT3",q$process)
q$process <- gsub("Myc","MYC",q$process)
q$process <- gsub("Mtorc1","MTORC1",q$process)
q$process <- gsub("Tnfa signaling via nfkb","TNFA signaling via NFKB",q$process)
q$process <- gsub("Il2 stat","IL2 STAT",q$process)
q <- q[order(q$NES),]
q$yax <- ifelse(q$NES > 0, -0.02, 0.02)
q$col <- ifelse(q$NES > 0, "blue","red")

fig3b <- 
  ggplot(q,aes(y=NES,x=reorder(process,NES),label = process))+
  geom_text(aes(y = yax,hjust = NES > 0),size=(figure_font_size)/(14/5))+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score (pCR)", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression: pCR",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(0.1,0,0.05,0), "lines"))
fig3b

fig3ab <- ggarrange(fig3a,fig3b,widths=c(1,1.6), labels = c("a","b"))


# b. Reactome enrichment
res <- results

# map Ensembl IDs to Entrez gene names
entrez <- read.table(paste0(resourcesDir,"EnsemblID.to.Entrez.tsv.gz"), sep="\t",header = T,stringsAsFactors = F)
entrez <- entrez[!is.na(entrez$EntrezGene.ID),]
res$symbol    <- mapIds(EnsDb.Hsapiens.v86, keys=row.names(res), column="GENENAME", keytype="GENEID", multiVals="first")
res$entrezENS <- mapIds(EnsDb.Hsapiens.v86, keys=row.names(res), column="ENTREZID", keytype="GENEID", multiVals="first")
res$entrezORG <- mapIds(org.Hs.eg.db, keys=row.names(res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res[is.na(res$entrezORG) & !is.na(res$entrezENS),"entrezORG"] <- res[is.na(res$entrezORG) & !is.na(res$entrezENS),"entrezENS"]

# remove duplicate entries
m      <- table(entrez[entrez$Gene.stable.ID %in% rownames(res[is.na(res$entrezORG),]),"Ensembl.ID"])
entrez <- entrez[entrez$Ensembl.ID %in% names(m[m==1]),]
res    <- merge(res,entrez,by.x=0,by.y="Ensembl.ID",all.x=T)
res[is.na(res$entrezORG) & !is.na(res$EntrezGene.ID),"entrezORG"] <- res[is.na(res$entrezORG) & !is.na(res$EntrezGene.ID),"EntrezGene.ID"]
res <- res[order(res$logFC, decreasing = T),]

geneList <- res[,2]
names(geneList) <- res[,9]
geneList <- geneList[!is.na(names((geneList)))]
gse <- gsePathway(geneList, nPerm=1000, minGSSize=120, pvalueCutoff=0.05, pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(gse)

# sort reactome enrichment by NES
e     <- res[order(res$NES),c(1,2,3,5,7)]
e$col <- ifelse(e$NES > 0, "blue","red")
e$yax <- ifelse(e$NES > 0, -0.02, 0.02)
e     <- e[abs(e$NES)>1.6,]
e$Description <- str_to_sentence(e$Description)
e$Description <- gsub("Dna","DNA",e$Description)
e$Description <- gsub("Esr","ESR",e$Description)
e$Description <- gsub("hiv","HIV",e$Description)
e$Description <- gsub("g1","G1",e$Description)
e$Description <- gsub("s phase","S phase",e$Description)
e$Description <- gsub("s trans","S trans",e$Description)
e$Description<-gsub("and","&",e$Description)

# Plot Extended Figure 5a: Reactome pathways associated with pCR vs RD
eFig5a <- 
  ggplot(e,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=(figure_font_size-1)/(14/5))+ 
  geom_bar(stat="identity",aes(fill=col),width = 0.9)+
  geom_hline(yintercept = 0)+
  coord_flip()+
  labs(y = "Normalised enrichment score", x = "",title="pCR vs residual disease")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression: pCR",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(0.1,0,0.2,0), "lines"))
eFig5a



#=========================================================================
# Pathways associated with increasing RCB class (Extended Figure 5b)
#=========================================================================

# remove any cases with pCR from within clinical metadata dataframe (p)
p.noPCR <- p[p$RCB.category!="pCR",]

# redo differential expression, this time to detect genes/pathways associated with increasing RCB score
transneo.counts.de2 <- data.frame(fread(paste0(dataDir,"transneo-diagnosis-RNAseq-rawcounts.tsv.gz"),
                                        header=T, sep="\t",stringsAsFactors = F),row.names = 1)
transneo.counts.de2 <- transneo.counts.de2[,colnames(transneo.counts.de2) %in% p.noPCR$Donor.ID]
p.noPCR <- p.noPCR[p.noPCR$Donor.ID %in% colnames(transneo.counts.de2),]
stopifnot(sum(p$Donor.ID!=colnames(transneo.counts))==0)
y.de2  <- DGEList(transneo.counts.de2)
minCPM <- 1
minNoToKeep <- 10
keep <- rowSums(cpm(y.de2)>minCPM)>=minNoToKeep
y.de2 <- y.de2[keep , , keep.lib.sizes=FALSE]
y.de2$samples$lib.size <- colSums(y.de2$counts)
y.de2 <- calcNormFactors(y.de2, method = "TMM")

# design matrix and perform DE
rcb    <- as.numeric(p.noPCR$RCB.score)
design <- model.matrix(~rcb, data = y.de2$samples)
coef   <- ncol(design)
y.de2  <- estimateDisp(y.de2, design = design, robust = TRUE)
fit    <- glmQLFit(y.de2, design, robust=TRUE)
qlf    <- glmQLFTest(fit, coef=coef)
results <- as.data.frame(topTags(qlf,n = Inf))
is.de   <- decideTestsDGE(qlf, p.value=0.05)


# Reactome gene set enrichment
res    <- results
entrez <- read.table(paste0(resourcesDir,"EnsemblID.to.Entrez.tsv.gz"), sep="\t",header = T,stringsAsFactors = F)
entrez <- entrez[!is.na(entrez$EntrezGene.ID),]
res$symbol    <- mapIds(EnsDb.Hsapiens.v86, keys=row.names(res), column="GENENAME", keytype="GENEID", multiVals="first")
res$entrezENS <- mapIds(EnsDb.Hsapiens.v86, keys=row.names(res), column="ENTREZID", keytype="GENEID", multiVals="first")
res$entrezORG <- mapIds(org.Hs.eg.db, keys=row.names(res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
res[is.na(res$entrezORG) & !is.na(res$entrezENS),"entrezORG"] <- res[is.na(res$entrezORG) & !is.na(res$entrezENS),"entrezENS"]

# remove duplicate entries
m      <- table(entrez[entrez$Gene.stable.ID %in% rownames(res[is.na(res$entrezORG),]),"Ensembl.ID"])
entrez <- entrez[entrez$Ensembl.ID %in% names(m[m==1]),]
res    <- merge(res,entrez,by.x=0,by.y="Ensembl.ID",all.x=T)
res[is.na(res$entrezORG) & !is.na(res$EntrezGene.ID),"entrezORG"]<-res[is.na(res$entrezORG) & !is.na(res$EntrezGene.ID),"EntrezGene.ID"]
res <- res[order(res$logFC, decreasing = T),]
geneList <- res[,2]
names(geneList)<-res[,9]
geneList<-geneList[!is.na(names((geneList)))]
gse <- gsePathway(geneList, minGSSize=100, pvalueCutoff=0.05, maxGSSize = 2000, pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(gse)

# sort reactome enrichment by NES
e     <- res[order(res$NES), c(1,2,3,5,7)]
e$col <- ifelse(e$NES > 0, "blue","red")
e$yax <- ifelse(e$NES > 0, -0.1, 0.1)
e     <- e[abs(e$NES)>1.5,]
e$Description<-gsub("and","&",e$Description)

eFig5b <- 
  ggplot(e,aes(y=NES,x=reorder(Description,NES), label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=(figure_font_size-1)/(14/5))+ 
  geom_bar(stat="identity",aes(fill=col),width = 0.9)+
  geom_hline(yintercept = 0)+
  coord_flip()+
  labs(y = "Normalised enrichment score", x = "",title="Increasing RCB class")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression: Increasing RD",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_manuscript()+
  theme(legend.position="bottom", 
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(0.1,0,0.2,0), "lines"))
eFig5b


#=========================================================================
# Association between tumour proliferation and response 
# Figure 3c, Extended Data Figure 6
#=========================================================================

rnadata <- merge(rnadata,p[,c("Donor.ID","RCB.category","pCR.RD","ER.status","HER2.status","ERHER2.status","Grade.pre.NAT")],
                 by="Donor.ID")
rnadata$RCB.category <- factor(rnadata$RCB.category,
                               levels=c("pCR","RCB-I","RCB-II","RCB-III"),
                               ordered=T)

# a. Associations with Genomic Grade Index

# GGI gsva score is correlated with Grade
# p=3e-13
summary(lm(rnadata$GGI.gsva~rnadata$Grade.pre.NAT))

fig3c_1 <- 
  ggplot(rnadata,aes(x=as.character(Grade.pre.NAT),y=GGI.gsva,fill=as.character(Grade.pre.NAT)))+
  geom_boxplot(outlier.size = 0.6, width=0.5)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_manual(values = rev(c("#FB6542","#375E97")))+
  stat_compare_means(ref.group=1,label = "p.format", hide.ns = T,size=5.5,label.y.npc = 0.9,color="#FB6542")+
  scale_y_continuous(breaks = c(-0.5,0,0.5))+
  labs(y="GGI score",x="Tumour Grade")+
  guides(fill="none")+
  theme(plot.margin = unit(c(1.5,0.1,0.5,1), "lines"),
        panel.grid.major.x = element_blank())
fig3c_1

# GGI score is monotonically associated with RCB
# p=2e-05
m      <- polr(factor(rnadata$RCB.category,ordered = T) ~ rnadata$GGI.gsva,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
pval[1]

z <- rnadata
z$RCB.category< - factor(z$RCB.category,levels=c("RCB-III","RCB-II","RCB-I","pCR"),labels=c("III","II","I","0"))
fig3c_2 <-
  ggplot(z,aes(x=GGI.gsva,y=(RCB.category), , fill = ..x..))+
  theme_manuscript(base_size = figure_font_size)+
  geom_density_ridges_gradient(scale=0.9, rel_min_height = 0.01)+
  labs(x="GGI score",y="RCB")+
  theme_ridges(grid=FALSE, center_axis_labels = T,font_size = figure_font_size)+
  scale_fill_viridis(option = "B") +
  scale_x_continuous(expand=c(0.01,0),breaks = c(-1,0,1))+
  scale_y_discrete(expand=c(0.01,0))+
  guides(fill="none")+
  theme(plot.margin = unit(c(1.5,0.1,0.21,1), "lines"),
        axis.text = element_text(size = figure_font_size),
        axis.line.x = element_line(color="black"))
fig3c_2


# GGI score only associated with response in HER2- tumours, p=4e-05
wilcox.test(
  rnadata[rnadata$HER2.status=="NEG" & rnadata$pCR.RD=="pCR","GGI.gsva"],
  rnadata[rnadata$HER2.status=="NEG" & rnadata$pCR.RD=="RD","GGI.gsva"])
# HER2+ p=0.99
wilcox.test(
  rnadata[rnadata$HER2.status=="POS" & rnadata$pCR.RD=="pCR","GGI.gsva"],
  rnadata[rnadata$HER2.status=="POS" & rnadata$pCR.RD=="RD","GGI.gsva"])

wilcox.test(
  rnadata[rnadata$RCB.category=="pCR","GGI.gsva"],
  rnadata[rnadata$RCB.category=="RCB-II","GGI.gsva"])

# b. Associations with Embryonic stem cell score

# ESC score is monotonically associated with the degree of RD, p=0.0001
m      <- polr(rnadata$RCB.category ~ rnadata$ESC.gsva,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
pval[1]

fig3c_3 <- 
  ggplot(rnadata,aes(x=RCB.category,y=ESC.gsva,fill=RCB.category))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_RCB()+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,label.y.npc = 0.9,color="#FB6542")+
  labs(y="Stem cell score",x="RCB category")+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  scale_y_continuous(limits = c(-0.8,0.8))+
  guides(fill="none")+
  theme(plot.margin = unit(c(1.5,0,0.5,1), "lines"),
        panel.grid.major.x = element_blank())
fig3c_3

fig3c <- ggarrange(fig3c_1,fig3c_2,fig3c_3,nrow = 1,widths=c(1,1,1))

eFig6a <- 
  ggplot(rnadata,aes(x=pCR.RD,y=GGI.gsva,fill=pCR.RD))+
  geom_boxplot(outlier.size = 0.6, width=0.6)+
  labs(x="Response",y="GSVA score",title=NULL)+
  facet_grid(~ERHER2.status)+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=4,label.y.npc = 0.94,color="#FB6542")+
  scale_fill_pCR_RD()+
  coord_cartesian(ylim=c(-0.9,0.9))+
  guides(fill="none")+
  theme_manuscript(base_size = figure_font_size)+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.5,0.5,1), "lines"))
eFig6a

table(rnadata$ERHER2.status)

# c. Associations with paclitaxel response metagene (Extended Figure 6b)

pacli <- rnadata[, grep("Donor.ID|Taxane",colnames(rnadata))]
pacli <- merge(pacli,p,by="Donor.ID")
pacli$HER2.status <- factor(pacli$HER2.status,levels=c("NEG","POS"),labels=c("HER2-","HER2+"))

eFig6b_1 <- 
  ggplot(pacli,aes(x=Taxane.MitosisScore,y=Taxane.CeramideScore,color=pCR.RD))+
  geom_point()+
  stat_ellipse(type = "t")+
  labs(x="Mitotic score",y="Ceramide score")+
  facet_grid(~HER2.status)+
  scale_colour_pCR_RD()+
  guides(color="none")+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 11),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))

eFig6b_2 <- 
  ggplot(data=pacli,aes(y=Taxane.FinalScore,x=RCB.category,fill=RCB.category))+
  geom_boxplot(outlier.size = 0.6, width=0.6)+
  facet_grid(~HER2.status,scales="free")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=4,label.y.npc = 0.9,color="#FB6542")+
  labs(y="Taxane score",x="RCB class")+
  scale_fill_RCB()+
  guides(fill="none")+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 11),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))

eFig6b <- 
  ggarrange(eFig6b_1,eFig6b_2, 
            labels = c("b",""),nrow = 2,
            font.label = list(size = figure_font_size, family="Helvetica"))
eFig6b

pdf(paste0(outputDir,"EFig6.pdf"),height=4.5,width=4,useDingbats = F)
eFig6b
dev.off()

wilcox.test(pacli[pacli$HER2.status=="HER2-" & pacli$pCR.RD=="pCR","Taxane.FinalScore",],
            pacli[pacli$HER2.status=="HER2-" & pacli$pCR.RD!="pCR","Taxane.FinalScore",])

#=========================================================================
# Immune associations with RCB (Figure 3d)
#=========================================================================

# a. Lymphocyte density
digPath <- read.table(paste0(dataDir,"transneo-diagnosis-DigPathology.tsv.gz"), header=T, stringsAsFactors = F,sep="\t")
digPath$median_lymph_KDE_knn_50 <- 10^(digPath$median_lymph_KDE_knn_50)

#reload metadata, to retain those with digital path+response
digPathMd <- metadataFull[metadataFull$RCB.category!="NA",]
digPathMd <- digPathMd[digPathMd$Chemo.cycles>1 & digPathMd$aHER2.cycles>1,]
im <- merge(digPath,digPathMd[,c("Donor.ID","RCB.category","pCR.RD","ERHER2.status")],by.x="Trial.ID",by.y="Donor.ID")

# Tumours that attained pCR had higher lymphocyte density p=0.0006
wilcox.test(im[im$pCR.RD=="pCR","median_lymph_KDE_knn_50"], im[im$pCR.RD!="pCR","median_lymph_KDE_knn_50"])

# Lymphocyte density monotonically associated with RCB category, p=0
im$RCB.category <- factor(im$RCB.category,ordered = T)
m      <- polr(im$RCB.category ~ im$median_lymph_KDE_knn_50,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
pval

fig3d_1 <- 
  ggplot(im,aes(x=RCB.category,y=median_lymph_KDE_knn_50,fill=RCB.category))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,label.y.npc = 0.9,color="#FB6542")+
  labs(y="Lymphocyte dens.",x="")+
  scale_fill_RCB()+
  guides(fill="none")+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  theme_manuscript(base_size = figure_font_size)+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        plot.margin = unit(c(2.5,0.1,0.5,1), "lines"))
fig3d_1


# b. CYT score

# CYT monotonically associated with RCB category, p=0.001
m      <- polr(rnadata$RCB.category ~ rnadata$CytScore.log2,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
pval[1]

fig3d_2 <-
  ggplot(rnadata,aes(x=RCB.category,y=CytScore.log2,fill=RCB.category))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_RCB()+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,label.y.npc = 0.9,color="#FB6542")+
  labs(y=expression(log[2]~CYT),x="RCB category")+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  guides(fill="none")+
  theme(plot.margin = unit(c(2.5,0.1,0.5,1), "lines"),
        panel.grid.major.x = element_blank())
fig3d_2


# LD correlated with the CYT (R2=0.4, p=1x10-15), 
cytIM <- merge(im[,c("Trial.ID","median_lymph_KDE_knn_50")],rnadata[,c("Donor.ID","CytScore.log2")],
               by.x="Trial.ID",by.y="Donor.ID")
c     <- cor.test(cytIM$CytScore.log2,cytIM$median_lymph_KDE_knn_50)
c$p.value
c$estimate^2


# c. Danaher T cell enrichment
danaher <- rnadata[,c("Donor.ID","ERHER2.status","RCB.category","pCR.RD",grep("Danaher",colnames(rnadata),value = T))]
fig3d_3 <- 
  ggplot(danaher,aes(x=RCB.category,y=Danaher.CD8.T.cells,fill=RCB.category))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_RCB()+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,label.y.npc = 0.9,color="#FB6542")+
  labs(y="CD8 T cells",x="")+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  guides(fill="none")+
  theme(plot.margin = unit(c(2.5,0,0.5,1), "lines"),
        panel.grid.major.x = element_blank())
fig3d_3

# T cell infiltrate monotonically associated with RCB category, p=0.0002
m      <- polr(rnadata$RCB.category ~ rnadata$Danaher.CD8.T.cells,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
pval[1]


fig3cd <- ggarrange(fig3c_1,fig3c_2,fig3c_3,fig3d_1,fig3d_2,fig3d_3,nrow=2,ncol=3,align = "v",heights=c(1,1.05))
fig3cd <- annotate_figure(fig3cd, fig.lab = "c   Proliferation",fig.lab.face = "bold",fig.lab.size = figure_font_size)

fig3abcd<-ggarrange(fig3ab,ggplot()+theme_void(),fig3cd,widths=c(1.2,0.04,1),nrow=1)

pdf(paste0(outputDir,"Fig3a-d.pdf"),height=12/2.54, width=41.5/2.54, useDingbats = F, onefile = T)
fig3abcd
dev.off()


#=========================================================================
# Other associations with immune metrics (Extended data figure 7)
#=========================================================================

# a Danaher immune estimation

danaher <- rnadata[,c("Donor.ID","ERHER2.status","RCB.category","pCR.RD",grep("Danaher",colnames(rnadata),value = T))]

# Plot Danaher PCA
pca <- danaher[,c(5:18)]
rownames(pca) <- danaher$Donor.ID
colnames(pca) <- gsub("Danaher\\.","",colnames(pca))
colnames(pca) <- gsub("\\."," ",colnames(pca))
counts  <- t(cpm(y,normalized.lib.sizes = T,log = T))
pca$ER <- counts[,"ENSG00000091831"]
pca <- prcomp(pca, center = TRUE,scale. = T) 

eFig7a <- 
  ggbiplot(pca,var.axes = T, choices = c(1,2),groups=danaher$pCR.RD,ellipse = F,
           varname.adjust = 1.2,point.size = 2, alpha=0.7,
           arrows.text.size = 2.8, 
           color = "black")+
  theme_minimal()+
  theme(axis.line=element_blank(),
        axis.ticks = element_line(size=0.2,color="black"),
        axis.text = element_text(color="black",size=11 ),
        axis.title = element_text(color="black",size=12 ),
        legend.position = "bottom",
        panel.grid = element_blank())+
  guides(fill="none",color="none")+
  scale_colour_pCR_RD()+
  labs(x=paste0("PC1 (",round(summary(pca)$importance[2,"PC1"]*100),"%)"),
       y=paste0("PC2 (",round(summary(pca)$importance[2,"PC2"]*100),"%)"))+
  coord_cartesian(ylim = c(-2.5,2), xlim=c(-2,2))+
  geom_segment(x=-2,xend=2,y=-2.72,yend=-2.72,size=0.2)+
  geom_segment(x=-2.2,xend=-2.2,y=-2,yend=2,size=0.2)
eFig7a

# plot Danaher box plots
d <- melt(danaher,id.vars=c(1:4))
d$variable<-factor(d$variable,levels=levels(d$variable),gsub("Danaher.","",levels(d$variable)))

d <- d[!d$variable %in% c("NK.CD56dim.cells","DC","CD45","Macrophages",
                          "Cytotoxic.cells","Exhausted.CD8","T.cells","Neutrophils","Th1.cells","Treg"),]

d$variable <- factor(d$variable, levels=levels(d$variable),labels=gsub("\\."," ", levels(d$variable)))

d$ERHER2.status<-factor(d$ERHER2.status,levels=c("ER- HER2-","ER+ HER2-","HER2+"),
                        labels=c("ER-\nHER2-","ER+\nHER2-","HER2+"))

eFig7b <-
  ggplot(d,aes(x=ERHER2.status,y=value,fill=pCR.RD))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Response")+
  stat_compare_means(aes(group=pCR.RD),label = "p.format", hide.ns = F,size=2.5,
                     color="#FB6542", label.y.npc = 0.94)+
  labs(y="Immune enrichment score",x="ER/HER2 status")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
eFig7b


# Mast cells were enriched in resistant tumours (enrichment score pCR: 2.1, RD: 3.4, p=0.0001; Wilcoxon rank sum test). 
mastPcr <- danaher[danaher$pCR.RD=="pCR","Danaher.Mast.cells"]
mastRd <- danaher[danaher$pCR.RD!="pCR","Danaher.Mast.cells"]
#pCR: 2.1
median(mastPcr)
#RD: 3.4
median(mastRd)
#p=0.0001
wilcox.test(mastPcr,mastRd)


# b. MCP counter immune estimation

mcp <- read.table(paste0(dataDir,"transneo-diagnosis-immune-MCPcounter.tsv.gz"), 
                  header=T, stringsAsFactors = F,sep="\t",row.names = 1)
mcp <- merge(p[,c("Donor.ID","ERHER2.status","pCR.RD","RCB.category")],mcp,by.x="Donor.ID",by.y=0)
mcp <- melt(mcp,id.vars = c("pCR.RD","RCB.category","ERHER2.status","Donor.ID"))
mcp <- mcp[!mcp$variable %in% c("T.cells","Fibroblasts","Endothelial.cells","Neutrophils"),]

mcp$ERHER2.status<-factor(mcp$ERHER2.status,levels=c("ER- HER2-","ER+ HER2-","HER2+"),
                          labels=c("ER-\nHER2-","ER+\nHER2-","HER2+"))

mcp$variable<-gsub("\\."," ",as.character(mcp$variable))

eFig7c <- 
  ggplot(mcp,aes(x=ERHER2.status,y=value,fill=pCR.RD))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free",nrow=1)+
  scale_fill_pCR_RD(name="Response")+
  labs(x="ER/HER2 status",y="MCPcounter enrichment score")+
  theme_manuscript(base_size = figure_font_size)+
  stat_compare_means(aes(group=pCR.RD),label = "p.format", hide.ns = T,size=3,color="#FB6542", label.y.npc = 0.95)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
eFig7c

cairo_pdf(paste0(outputDir,"EFig7c.pdf"),height=9.3/2.54, width=32/2.54)
eFig7c
dev.off()


# c. Immunophenoscore

ips <- read.table(paste0(dataDir,"transneo-diagnosis-immune-IPS-components.tsv.gz"), header=T, stringsAsFactors = F,sep="\t",row.names = 1)
ips <- ips[rownames(ips) %in% p$Donor.ID,]

colnames(ips) <- gsub("\\."," ",colnames(ips))
colnames(ips) <- gsub("HLA ","HLA-",colnames(ips))
colnames(ips) <- gsub("PD ","PD-",colnames(ips))
colnames(ips) <- gsub("CTLA ","CTLA-",colnames(ips))

components <- read.table(paste0(resourcesDir,"IPS_genes.txt"),sep="\t",header = T)[,c(2,3)]
components <- rbind(components[components$CLASS=="MHC",],
                    components[components$CLASS=="EC",],
                    components[components$CLASS=="CP",],
                    components[components$CLASS=="SC",])
components=components[!duplicated(components),]
rownames(components) <- components[,1]
components <- components[,-1,drop=F]
colnames(components)[1] <- "Class"
components$Class <- factor(components$Class,levels=c("MHC","EC","CP","SC"),
                           labels=c("MHC","Effector cells","Checkpoints","Suppressor cells"))

ips <- ips[,match(rownames(components),colnames(ips))]

ips[,components$Class %in% c("Checkpoints","Suppressor cells")] <- ips[,components$Class %in% c("Checkpoints","Suppressor cells")]*-1
ips[,"CD27"] <- ips[,"CD27"]*-1
ips[,"ICOS"] <- ips[,"ICOS"]*-1

annotation <- data.frame(Response = p$RCB.category)
rownames(annotation) <- p$Donor.ID

Response        <- c("#20A39E","#ffe671","#fdb462","#ef3b2c")
names(Response) <- c("pCR","RCB-I","RCB-II","RCB-III")
anno_colors     <- list(Response = Response)

eFig7d <- pheatmap.type(data.frame(t(ips)), components, type = colnames(components)[1],
                        doTranspose=F, conditions="Auto",
                        show_colnames = F,treeheight_row = 0,cellwidth = 4,cellheight = 12,
                        legend = T, treeheight_col = 10, annotation_legend = F,
                        main = "Immunoscore enrichment",
                        clustering_method = "ward.D2",
                        clustering_distance_cols = "euclidean",
                        annotation_col = annotation,
                        annotation_colors = anno_colors,
                        scale = "row",
                        color = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                        fontsize_col = 12,
                        fontsize_row = 12,
                        fontsize = 12,
                        cutree_cols = 3, gaps_row = c(10,14,24))
eFig7d <- eFig7d$pheat$gtable
eFig7d$grobs[[1]]$gp$fontsize=12
ggarrange(eFig7d)


# d. Association between Danaher immune infiltration and lymphoyctic density
#    Extended Data Figure 7e

cytIM <- merge(im[,c("Trial.ID","median_lymph_KDE_knn_50")],rnadata,by.x="Trial.ID",by.y="Donor.ID")
d <- melt(cytIM,measure.vars=grep("Danaher",colnames(cytIM)))
d <- d[!d$variable%in% c("Danaher.T.cells","Danaher.NK.CD56dim.cells","Danaher.Cytotoxic.cells","Danaher.CD45","Danaher.DC"),]

d$variable <- factor(d$variable,levels=levels(d$variable),labels=gsub("\\."," ", gsub("Danaher.","",levels(d$variable))))

eFig7e <- 
  ggplot(d[d$variable!="Treg",],aes(x=median_lymph_KDE_knn_50,y=value))+
  geom_point(size=0.8)+
  facet_wrap(~variable,scales="free_y",nrow=2)+
  labs(x="Lymphocyte density (digital pathology)",y="Abdundance (RNA-seq)")+
  geom_smooth(method="lm", formula = y~x, se = T)+
  scale_x_continuous(breaks=c(0,0.0005,0.001),labels = comma)+
  theme_manuscript(base_size = 11)+
  stat_poly_eq(aes(label = paste(..rr.label..)), color="red",
               label.x.npc = "left", label.y.npc = 0.98,
               formula = y~x, parse = TRUE, size = 3)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 12),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))
eFig7e


#=========================================================================
# Immunity and Proliferation interplay (Figure 3e)
#=========================================================================

pdf(paste0(outputDir,"Fig3e.pdf"),width=(5.38*2)/2.54,height=(5.1*2)/2.54, useDingbats = F, onefile = T)

par(mfrow=c(2,2), oma = c(2, 2, 0, 0),
    mar=c(2,1.5,1.5,1), font.lab=2,
    mgp = c(2.5, 1, 0), cex.axis=1.2, 
    cex.main=1.2, cex.lab=1.5)

for (q in c("pCR","RCB-I","RCB-II","RCB-III")){
  GGI   <- rnadata[rnadata$RCB.category %in% q,"GGI.gsva"]
  STAT1 <- rnadata[rnadata$RCB.category %in% q,"STAT1.gsva"]
  y     <- cbind(`GGI score`=GGI,`STAT1 score` = STAT1)
  sm.density(y,display="image",panel=F,ylim=c(-1,1),xlim=c(-1,1),xlab="",ylab="")
  title(q, line = 0.7)
  abline(v=mean(rnadata$GGI.gsva),lty="dashed")
  abline(h=mean(rnadata$STAT1.gsva),lty="dashed")
}
mtext("Immune score (STAT1)",side=2,line=1,outer=TRUE,cex=1,las=0,font=1)
mtext("Proliferation score (GGI)",side=1,line=1,outer=TRUE,cex=1,font=1)
dev.off()


# Validation of relationship between immune activation and proliferation
# Extended Data Figure 7f 

extValid <- read.table(paste0(resourcesDir,"GGI-STAT1-enrichment.tsv"),  sep="\t",header=T)

pdf(paste0(outputDir,"EFig7f.pdf"),width=(5.38*2)/2.54,height=(5.1*2)/2.54, useDingbats = F, onefile = T)
par(mfrow=c(2,2),oma = c(2, 2, 0, 0),
    mar=c(2,1.5,1.5,1), font.lab=2,
    mgp = c(2.5, 1, 0), cex.axis=1.2, 
    cex.main=1.2, cex.lab=1.5)

for (q in c("pCR","RCB-I","RCB-II","RCB-III")){
  GGI   <- extValid[extValid$class %in% q,"ggi"]
  STAT1 <- extValid[extValid$class %in% q,"stat1"]
  y     <- cbind(`GGI score`=GGI,`STAT1 score` = STAT1)
  sm.density(y,display="image",panel=F,ylim=c(-1,1),xlim=c(-1,1),xlab="",ylab="")
  title(q, line = 0.7)
  abline(v=-0,lty="dashed")
  abline(h=0,lty="dashed")
}
mtext("Immune score (STAT1)",side=2,line=1,outer=TRUE,cex=1,las=0,font=1)
mtext("Proliferation score (GGI)",side=1,line=1,outer=TRUE,cex=1,font=1)
dev.off()

