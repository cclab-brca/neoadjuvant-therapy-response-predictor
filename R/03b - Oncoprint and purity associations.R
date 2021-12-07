# Manuscript:   Multi-omic machine learning predictor of breast cancer therapy response
# Author:       Stephen-John Sammut
# Description:  Oncoprint and purity associations
# Section:      Results - Genomic landscapes associate with response
#=======================================================================================

rm (list=ls())

library (gridExtra)
library (lemon)
library (MASS)
library (readxl)
library (reshape2)


# change the contents of variable baseDir to the root analysis folder 
baseDir <- "~/Neoadjuvant-predictor/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/R/LoadDirectoryStructure.R"))

figure_font_size <- 12
oncoprint_base_size <- 10

# load clinical metadata (Supplementary Table 1)
metadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 1))
# combine ER and HER2 status
metadata$ERHER2.status<-ifelse(metadata$ER.status=="POS","ER+ HER2-","ER- HER2-")
metadata$ERHER2.status<-ifelse(metadata$HER2.status=="POS","HER2+",metadata$ERHER2.status)

# load mutation data (Supplementary Table 2)
mutations  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 2))

# load driver gene list
driverGenes<- scan(paste0(resourcesDir,"breast-cancer-driver-genes.txt"), what=character(),skip = 1)

#=========================================================================
# Oncoprint
#=========================================================================

p <- metadata[metadata$RCB.category!="NA",]
m <- mutations[,c(1,4,8)]
m <- m[m$Donor.ID %in% p$Donor.ID,]
allPatients <- unique(m$Donor.ID)

m$classification<-4
m[m$MAF_Variant %in% c("Missense_Mutation"),"classification"] <- 1
m[m$MAF_Variant %in% c("Nonsense_Mutation", "Splice_Site","Frame_Shift_Del", 
                       "Frame_Shift_Ins"), "classification"] <- 2
m[m$MAF_Variant %in% c("In_Frame_Del", "In_Frame_Ins"), "classification"] <- 3

# subset driver mutations
all.driver.mutations <- m[m$Hugo_Symbol %in% driverGenes,]

# create mutation matrix
mutations.subset <- all.driver.mutations

mutMatrix <- matrix("", nrow = length(unique(mutations.subset$Hugo_Symbol)),ncol=length(allPatients), 
                    dimnames = list(c(unique(mutations.subset$Hugo_Symbol)),allPatients))

for (r in 1:nrow(mutations.subset)){
  mut <- mutMatrix[mutations.subset[r,"Hugo_Symbol"],mutations.subset[r,"Donor.ID"]]
  if (mut=="")  {
    mutMatrix[mutations.subset[r,"Hugo_Symbol"],mutations.subset[r,"Donor.ID"]] <- mutations.subset[r,"classification"]
  } else {
    mutMatrix[mutations.subset[r,"Hugo_Symbol"],mutations.subset[r,"Donor.ID"]] <- paste0(mut,";",mutations.subset[r,"classification"])
  }
}

#remove mutations present in fewer than 2 samples
mutMatrix<-mutMatrix[!apply(mutMatrix,1, function(x) sum(x==""))>(ncol(mutMatrix)-2),]

# sort Matrix in an oncoprint style
sortMatrix <- mutMatrix
sortMatrix[sortMatrix!=""]<-1
sortMatrix[sortMatrix==""]<-0
sortMatrix<-apply(sortMatrix,2,as.numeric)
rownames(sortMatrix)<-rownames(mutMatrix)

#sort rows
rowOrder<- order(rowSums(sortMatrix), decreasing = TRUE)

#sort columns
scoreCol = function(x) {
  score = 0
  for (i in 1:length(x)) {
    if (x[i]) {
      score = score + 2^(length(x) - i * 1/x[i])
    }
  }
  return(score)
}

colOrder<-character()
for (r in c("pCR","RCB-I","RCB-II","RCB-III")){
  s <- sortMatrix[,colnames(sortMatrix) %in% p[p$RCB.category==r,"Donor.ID"]]
  scores <- apply(s[order(rowSums(s), decreasing = TRUE), , drop = FALSE], 2, scoreCol)
  colOrder<-append(colOrder,colnames(s)[order(scores, decreasing = TRUE)])
}

mutMatrix<-mutMatrix[rowOrder,match(colOrder,colnames(mutMatrix))]

mutMatrix2 <- mutMatrix
allM       <- reshape2::melt(mutMatrix2)
allM$value <- as.character(allM$value)
allM       <- merge(allM,p,by.x="Var2",by.y="Donor.ID")
allM$Var1  <- factor(allM$Var1, levels=rev(rownames(mutMatrix2)))
allM$Var2  <- factor(allM$Var2, levels=colnames(mutMatrix2))
f <- allM
f <- f[grep(";",f$value),]

allM[grep("^2;2$|^2;1$|^1;2$|^2;3$",allM$value),"value"]<-2
allM[grep("^1;1$|^1;4$|^4;1$|^1;1;1$|^1;1;4;1$",allM$value),"value"]<-1
allM[grep("^3;4$",allM$value),"value"]<-3
allM[grep("^4;4$",allM$value),"value"]<-4
allM[grep("^4;1;2$|^2;4;1$",allM$value),"value"]<-2

p1<-
  ggplot(allM, aes(Var2, Var1)) +
  geom_tile(aes(fill = value),colour = "white", size=0.6) +
  geom_point(data=f,aes(x=Var2,y=Var1),pch=4,color="white",size=2)+
  scale_fill_manual(name="Mutation classification",values = c('#d7edf5','#217CA3','#D70026','#F5BE41',"#258039","red"),
                    breaks=c(1,2,3,4),
                    labels=c("Missense","Truncating","In-frame","Other"))+
  
  theme_grey(base_size = oncoprint_base_size) + 
  scale_x_discrete(expand = c(0, 0), aes(limits=Donor.ID)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x=element_text(size = oncoprint_base_size*0.7, vjust=0.5,angle = 90, hjust = 1,color="black"),
        axis.text.y=element_text(size = oncoprint_base_size*0.7,face="italic",color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=oncoprint_base_size*1.1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank())+    
  facet_grid(~RCB.category, scales = "free_x",space = "free_x")

m2 <- merge(x=all.driver.mutations,y=p,by="Donor.ID",all.y = T)
m2 <- m2[,c("Donor.ID","classification","ERHER2.status","RCB.category")]
m2$Donor.ID <- factor(m2$Donor.ID, levels=colnames(mutMatrix2))
m2$classification <-factor(m2$classification, levels=c(1,2,3,4,5))

p2<-
  ggplot(m2, aes(Donor.ID, fill=classification)) +
  geom_bar()+
  theme_grey(base_size = oncoprint_base_size) + 
  labs( y = "Mutations in\ndriver genes") + 
  scale_fill_manual(values = c('#217CA3','#D70026','#F5BE41',"#258039","red"))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(name = waiver(), breaks=c(0,2,4,6,8), limits = c(0,9))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = oncoprint_base_size*0.8,face="bold"),
        
        legend.position = "none",
        axis.text.y=element_text(size = oncoprint_base_size*0.7, colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank())+    
  facet_grid(~RCB.category, scales="free",space="free",drop = TRUE)

m3 <- merge(x=m,y=p,by="Donor.ID")
m3 <- m3[,c("Donor.ID","classification","ERHER2.status","RCB.category")]

m2$Donor.ID<- factor(m2$Donor.ID, levels=colnames(mutMatrix2))
m2$classification <- factor(m2$classification, levels=c(1,2,3,4,5))

m3 <- m3[m3$Donor.ID %in% colnames(mutMatrix2),]
m3$Donor.ID <- factor(m3$Donor.ID, levels=colnames(mutMatrix2))
m3$classification <-factor(m3$classification, levels=c(1,2,3,4,5))

p3 <- 
  ggplot(m3, aes(Donor.ID, fill=classification)) +
  geom_bar()+
  theme_grey(base_size = oncoprint_base_size) + 
  labs(y = "Mutations in\nall genes") + 
  scale_fill_manual(values = c('#217CA3','#D70026','#F5BE41',"#258039","red"))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(name = waiver(), breaks=waiver())+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none",
        axis.text.y=element_text(size = oncoprint_base_size*0.7, colour = "black"),
        axis.title.y=element_text(size = oncoprint_base_size*0.8,face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = oncoprint_base_size,face = "bold"),
        strip.text.y = element_text(size = oncoprint_base_size,face = "bold")) + 
  facet_grid(~RCB.category, scales="free",space="free",drop = TRUE)

gA <- ggplot_gtable(ggplot_build(p1))
gB <- ggplot_gtable(ggplot_build(p2))
gC <- ggplot_gtable(ggplot_build(p3))

maxWidth = grid::unit.pmax(gA$widths, gB$widths,gC$widths)
gA$widths <- as.list(maxWidth)
gB$widths <- as.list(maxWidth)
gC$widths <- as.list(maxWidth)

eFig3 <- grid.arrange(
  arrangeGrob(gC,gB,gA,nrow=3,heights=c(.13,.07,.8))
)

pdf(paste0(outputDir,"EFig3.pdf"), onefile=T, useDingbats = F,height=20/2.54, width=36.6/2.54)
ggarrange(eFig3)
dev.off()


#=========================================================================
# Associations with purity
#=========================================================================

#retain cases with RCB data and received >1 cycle of aHER2/chemotherapy
metadata <- metadata[metadata$RCB.category!="NA",]
metadata <- metadata[metadata$Chemo.cycles>1 & metadata$aHER2.cycles>1,]

purity <- read.table(paste0(dataDir,"transneo-diagnosis-ASCAT-purity.tsv.gz"),header = T,stringsAsFactors = F)
purity <- merge(purity,p[,c("Donor.ID","RCB.category","pCR.RD","HER2.status")],by.x="Trial.ID",by.y="Donor.ID")
purityRCB <- purity[purity$Trial.ID %in% metadata$Donor.ID,]

#purity was not associated with ordinal RCB, p = 0.1
rcb    <- factor(purityRCB$RCB.category,labels=c("pCR","RCB-I","RCB-II","RCB-III"), ordered = T)
m      <- polr(rcb ~ purityRCB$purity,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
round(as.numeric(pval[1]),2)

# TMB and %CNA were both independent of tumour purity
totRegion <- 45.54094
tmb <- data.frame(table(mutations$Donor.ID)/totRegion)
purityTMB <- merge(purity,tmb,by.x="Trial.ID",by.y="Var1")

#CIN
cin <- read.table(paste0(dataDir,"transneo-diagnosis-ASCAT-CIN.tsv.gz"), sep="\t",header = T, stringsAsFactors = F, quote = "")
purityTMBCIN <- merge(purityTMB,cin,by="Trial.ID")

t1 <- cor.test(purityTMBCIN$Freq,purityTMBCIN$purity)
t2 <- cor.test(purityTMBCIN$cin,purityTMBCIN$purity)

eFig4_b1 <-
  ggplot(purityTMBCIN,aes(x=purity,y=Freq))+
  geom_point(size=1)+
  annotate("text", y = 12, x = 0.5, label = paste0("italic(R)==",round(t1$estimate,2),"~~italic(p)==",round(t1$p.value,2)),parse=T, size=3)+
  stat_smooth(method = "lm")+
  labs(x="Tumour purity",y="TMB")+
  scale_y_continuous(limits=c(0,12))+
  theme_manuscript(base_size = figure_font_size)+
  theme(plot.margin = unit(c(1.5,1,0.5,1), "lines"))
eFig4_b1


eFig4_b2 <-
  ggplot(purityTMBCIN,aes(x=purity,y=cin))+
  geom_point(size=1)+
  stat_smooth(method = "lm")+
  annotate("text", y = 1, x = 0.5, label = paste0("italic(R)==",round(t2$estimate,2),"~~italic(p)==",round(t2$p.value,2)),parse=T, size=3)+
  labs(x="Tumour purity",y="%CNA")+
  theme_manuscript(base_size = figure_font_size)+
  theme(plot.margin = unit(c(1.5,1,0.5,1), "lines"))
eFig4_b2


