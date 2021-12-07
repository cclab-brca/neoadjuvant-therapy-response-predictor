# Manuscript:   Multi-omic machine learning predictor of breast cancer therapy response
# Author:       Stephen-John Sammut
# Description:  Transcriptomic features associated with response
# Section:      Immune dysfunction in resistant tumours
#=======================================================================================

rm (list=ls())

library (data.table)
library (edgeR)
library (fgsea)
library (reshape2)
library (readxl)
library (stringr)

# change the contents of variable baseDir to the root analysis folder 
baseDir <- "~/Neoadjuvant-predictor/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/R/LoadDirectoryStructure.R"))

figure_font_size=12

# load MSigDB Hallmarks gene sets
hallmarks <- paste0(resourcesDir,"MSigDB.HallmarksGeneSet.Rdata")
hallmarksGeneList <- readRDS(hallmarks)

# load clinical metadata (Supplementary Table 1)
metadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 1))
# combine ER and HER2 status
metadata$ERHER2.status<-ifelse(metadata$ER.status=="POS","ER+ HER2-","ER- HER2-")
metadata$ERHER2.status<-ifelse(metadata$HER2.status=="POS","HER2+",metadata$ERHER2.status)

# Perform analyses with cases that had an RCB assessment and received more than one cycle of therapy
# as detailed in Methods (Statistical testing)
metadata <- metadata[metadata$RCB.category!="NA",]
metadata <- metadata[metadata$Chemo.cycles>1 & metadata$aHER2.cycles>1,]

# load RNA data (Supplementary Table 3)
rnadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 3))
rnadata  <- merge(metadata,rnadata,by="Donor.ID")
nrow(rnadata)

#=========================================================================
# Identification of immune-high and proliferation-high tumours that do not
# respond to therapy (Figure 3f)
#=========================================================================

g <- rnadata

# set thresholds
ggiThreshold <- mean(g$GGI.gsva)
stat1Threshold <- mean(g$STAT1.gsva)

#26 cases have RD post therapy despite having high STAT1 and GGI 
# and receiving >1 cycle of therapy
table(g[g$STAT1.gsva>stat1Threshold & g$GGI.gsva>ggiThreshold,"pCR.RD"])


fig3f_1 <- 
  ggplot(g,aes(x=GGI.gsva,y=STAT1.gsva,color=pCR.RD))+
  scale_x_continuous(limits=c(-0.83,0.83))+
  annotate("rect", xmin=0, xmax=0.7,ymin=0, ymax=0.83,fill="steelblue",color="steelblue", linetype="dashed",alpha=0.2)+
  geom_point(size=1.2)+
  labs(x="GGI",y="STAT1")+
  scale_colour_pCR_RD()+
  theme_bw()+
  guides(color="none")+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")))
fig3f_1

transneo.counts <- data.frame(fread(paste0(dataDir,"transneo-diagnosis-RNAseq-rawcounts.tsv.gz"),
                                    header=T, sep="\t",stringsAsFactors = F),row.names = 1)

doDE <- function(metadata, counts){
  counts <- counts[,colnames(counts) %in% metadata$Donor.ID]
  # construct DE matrix
  y <- DGEList(counts)
  minCPM  <- 1
  minNoToKeep <- 2
  keep <- rowSums(cpm(y)>minCPM)>=minNoToKeep
  y <- y[keep , , keep.lib.sizes=FALSE]
  y$samples$lib.size <- colSums(y$counts)
  y <- calcNormFactors(y, method = "TMM")
  
  # construct design matrix
  rcb <- factor(metadata$pCR.RD, levels=rev(c("RD","pCR")))
  er <- as.factor(metadata$ER.status)
  her2 <- factor(metadata$HER2.status)
  # design matrix
  design <- model.matrix(~er+her2+rcb, data = y$samples)
  coef <- ncol(design)
  
  # perform DE
  y   <- estimateDisp(y, design = design, robust = TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  qlf <- glmQLFTest(fit, coef=coef)
  results <- as.data.frame(topTags(qlf,n = Inf))
  is.de <- decideTestsDGE(qlf, p.value=0.05)
  indexes <- ids2indices(hallmarksGeneList, rownames(y$counts))
  gst.camera <- camera.DGEList(y, index=indexes, design=design, coef=coef)
  q <- gst.camera[gst.camera$FDR<0.01,]
  q$process <- gsub("HALLMARK_","",rownames(q))
  q$process <- gsub("_"," ",q$process)
  q[q$process =="EPITHELIAL MESENCHYMAL TRANSITION","process"]<-"EMT"
  q$process<-str_to_sentence(q$process)
  q$process<-gsub("E2f","E2F",q$process)
  q$process<-gsub("Interferon","IFN",q$process)
  q$process<-gsub("Emt","EMT",q$process)
  q$process<-gsub("G2m","G2M",q$process)
  q$process<-gsub("Myc","MYC",q$process)
  q$process<-gsub("Mtorc1","MTORC1",q$process)
  q$process<-gsub("Tnfa signaling via nfkb","TNFA signaling via NFKB",q$process)
  
  return (q)
}

# Do a DE on samples that have high prolif and high immune activation
g.hihi <- g[g$GGI.gsva>ggiThreshold & g$STAT1.gsva>stat1Threshold,]
enrich.result <- doDE(g.hihi,transneo.counts)

# Plot MSigDB Hallmarks enrichment
fig3f_2 <- ggplot(enrich.result,aes(x=reorder(process,FDR),y=-log10(FDR),color=Direction))+
  geom_point(size=3)+
  labs(x="",y=expression(-log[10]~FDR~(RD)))+
  scale_colour_manual(name="Expression", values = c("#FB6542","#375E97"), limits=c("Up","Down"))+
  scale_size(name="Number of Genes")+
  coord_flip()+
  theme_manuscript(base_size = 11.5)+
  theme(panel.grid.major.y = element_blank(),
        strip.background =element_blank(),
        plot.margin = unit(c(0.1,0.2,0,0), "lines"))+
  guides(color="none")
fig3f_2


#=========================================================================
# TIDE dysfunction and exclusion scores in immune/proliferation 
# high tumours (Figure 3f)
#=========================================================================

tide <- rnadata
tide$HER2.status <- factor(tide$HER2.status,levels=c("NEG","POS"),labels=c("HER2-","HER2+"))

#T cell dysfunction scores stratified by response and across HER2+ and - subtypes
fig3f_3_a <- 
  ggplot(tide[tide$STAT1.gsva>stat1Threshold & tide$GGI.gsva>ggiThreshold,],
         aes(x=pCR.RD,y=TIDE.Dysfunction,fill=pCR.RD))+
  geom_boxplot(width=0.6, outlier.size = 0.6)+
  labs(x="",y="T cell dysfunction")+
  scale_fill_pCR_RD()+
  facet_grid(~HER2.status)+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,label.y.npc = 0.9,color="#FB6542")+
  theme_manuscript(base_size = 12)+
  guides(fill="none")+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 11),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(-0.1,0,-0.5,1), "lines"))
fig3f_3_a

#T cell dysfunction is associated with response in HER2- tumours, p=0.006
a <- tide[tide$RCB.category=="pCR" & tide$STAT1.gsva> stat1Threshold & tide$GGI.gsva>ggiThreshold 
          &tide$HER2.status=="HER2-","TIDE.Dysfunction"]
b <- tide[tide$RCB.category!="pCR" & tide$STAT1.gsva> stat1Threshold & tide$GGI.gsva>ggiThreshold 
          &tide$HER2.status=="HER2-","TIDE.Dysfunction"]
wilcox.test(a,b)

#T cell exclusion scores stratified by response and across HER2+ and - subtypes
fig3f_3_b <- ggplot(tide[tide$STAT1.gsva>stat1Threshold & tide$GGI.gsva>ggiThreshold,],
                    aes(x=pCR.RD,y=TIDE.Exclusion,fill=pCR.RD))+
  geom_boxplot(width=0.6, outlier.size = 0.6)+
  labs(x="",y="T cell exclusion")+
  scale_fill_pCR_RD()+
  facet_grid(~HER2.status)+
  stat_compare_means(ref.group="pCR",label = "p.signif", hide.ns = T,size=5.5,label.y.npc = 0.9,color="#FB6542")+
  theme_manuscript(base_size = 12)+
  guides(fill="none")+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 11),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(-0.1,0,-0.5,1), "lines"))
fig3f_3_b


o1 <- ggarrange(ggplot()+theme_void(),fig3f_1,ggplot()+theme_void(),nrow=3, heights=c(0.2,0.4,0.2))
o2 <- ggarrange(fig3f_3_a,fig3f_3_b,nrow=1,align = "hv")
o3 <- ggarrange(fig3f_2,ggplot()+theme_void(),o2,nrow=3,heights=c(1.05,0.05,1))
fig3f <- ggarrange(o1,o3,nrow=1,widths=c(0.4,1),labels = "e",
                   font.label = list(size = 12, family="Helvetica"))

pdf(paste0(outputDir,"Fig3f.pdf"),height=11.5/2.54, width=18/2.54, useDingbats = F, onefile = T)
fig3f
dev.off()


# Identify other immunosuppressive cell types in Danaher's gene sets 
# Extended Data Figure 8

m  <- tide[tide$STAT1.gsva>stat1Threshold & tide$GGI.gsva>ggiThreshold & tide$HER2.status=="HER2-",]
md <- melt(m,measure = c("Danaher.T.cells","Danaher.NK.CD56dim.cells","Danaher.Treg"))
md$variable <- factor(md$variable,levels=c("Danaher.T.cells","Danaher.NK.CD56dim.cells","Danaher.Treg"),
                      labels=c("T-cells","NK CD56 dim","T-reg"))

eFig8a<- ggplot(md,aes(x=pCR.RD,y=value,fill=pCR.RD))+
  geom_boxplot(width=0.5, outlier.size = 0.7)+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = F,size=5,color="#FB6542", label.y.npc = 0.96)+
  theme_manuscript()+
  labs(y="Score",x="Response",title="HER2- tumours, high GGI and STAT1 (n=28)")+
  scale_fill_pCR_RD()+
  facet_wrap(~variable,scales="free",nrow=1)+
  guides(fill="none")+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 11),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.1,1,0.5,1), "lines"))
eFig8a

#NK CD56 dim cells associated with response in HER2- tumours, p=0.01
wilcox.test(md[md$variable=="NK CD56 dim" & md$pCR.RD=="pCR","value"],
            md[md$variable=="NK CD56 dim" & md$pCR.RD!="pCR","value"])

#T regs associated with response in HER2- tumours, p=0.02
wilcox.test(md[md$variable=="T-reg" & md$pCR.RD=="pCR","value"],
            md[md$variable=="T-reg" & md$pCR.RD!="pCR","value"])

table(m[m$HER2.status=="HER2-","pCR.RD"])

#=========================================================================
# T cell exclusion metrics
#=========================================================================

colnames(tide) <- gsub("TIDE.","",colnames(tide))
toKeep <- c("Exclusion","CAF","TAM.M2","MDSC")

h <- reshape2::melt(tide,measure.vars=toKeep)
h$variable <- factor(h$variable,levels=toKeep,labels=toKeep)

eFig8b <-
  ggplot(h,aes(x=pCR.RD,y=value))+
  geom_boxplot(outlier.size = 0.6, width=0.5,aes(fill=pCR.RD))+
  geom_hline(yintercept = 0, linetype="dotted")+
  stat_compare_means(ref.group = 1, label = "p.format",color="tomato",hide.ns =F, label.y.npc = 0.96,size=3)+
  labs(x="Response",y="Enrichment",title="All tumours (n=149)")+
  theme_manuscript(base_size = 12)+
  scale_fill_pCR_RD()+
  guides(fill="none")+
  facet_wrap(~variable,scales="free_y",ncol=4)+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 12),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.1,1,0.5,1), "lines"))
eFig8b

#Tumours that did not attain pCR had higher levels of T-cell exclusion 
#p=0.02
wilcox.test(tide[tide$RCB.category=="pCR" ,"Exclusion"],tide[tide$RCB.category!="pCR","Exclusion"])

#CAF: p=0.009
wilcox.test(tide[tide$RCB.category=="pCR","CAF"],tide[tide$RCB.category!="pCR","CAF"])

#TAM M2: p=0.0009
wilcox.test(tide[tide$RCB.category=="pCR","TAM.M2"],tide[tide$RCB.category!="pCR","TAM.M2"])

pdf(paste0(outputDir,"eFig7.pdf"),height=8/2.54, width=30/2.54, useDingbats = F, onefile = T)
ggarrange(eFig8a,eFig8b,widths=c(1,1.5),labels=c("a","b"))
dev.off()



