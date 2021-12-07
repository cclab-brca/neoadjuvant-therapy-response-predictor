# Manuscript:   Multi-omic machine learning predictor of breast cancer therapy response
# Author:       Stephen-John Sammut
# Description:  Transcriptomic features associated with response
# Section:      Machine learning integrates multi-omic features
#=======================================================================================

rm (list=ls())

library (data.table)
library (MASS)
library (readxl)
library (reshape2)
library (rstatix)

# change the contents of variable baseDir to the root analysis folder 
baseDir <- "~/Neoadjuvant-predictor/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/R/LoadDirectoryStructure.R"))

# combine metadata used for training and validation
metadata.train  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 1))[,c("Donor.ID","pCR.RD","RCB.category")]
metadata.valid  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 4))[,c("Donor.ID","pCR.RD","RCB.category")]
metadata<-rbind(metadata.valid,metadata.train)

# load ml predictor scores
mlscores <- data.frame(fread(paste0(dataDir,"transneo-diagnosis-MLscores.tsv.gz")),stringsAsFactors = F)

l <- merge(y=metadata,x=mlscores,by.y="Donor.ID",by.x="Trial.ID")
l <- melt(l,measure.vars=c(3:8))
l$variable <- factor(l$variable, levels=c(
  "Clinical","Clinical.DNA", "Clinical.RNA",
  "Clinical.DNA.RNA","Clinical.DNA.RNA.DigPath","Integrated.model"),
  labels=c("Clinical","Clinical+DNA","Clinical+RNA","Clinical+DNA+RNA","Clinical+DNA+RNA+DigPath","Integrated model"))

colnames(l)[5]<-"pipeline"
l=l[l$RCB.category!="NA",]
l$RCB.category<-factor(l$RCB.category,levels=c("pCR","RCB-I","RCB-II","RCB-III"))

stat.test <- l %>%
  group_by(Class, pipeline) %>%
  wilcox_test(value ~ RCB.category)

stat.test<-stat.test[stat.test$group1=="pCR",] %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "RCB.category")

stat.test$p.adj.short<- signif(stat.test$p.adj,2)

eFig9 <- 
  ggboxplot(
    l, x = "RCB.category", y = "value", fill="RCB.category",outlier.size=0.5,
    facet = c("Class", "pipeline")) +
  stat_pvalue_manual(stat.test, hide.ns = F,label = "p.adj.short")+
  labs(x="RCB class",y="Predictor score")+
  scale_fill_RCB()+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  guides(fill="none")+
  theme_bw(base_size = 12)+
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold",size = 10.5),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(colour = "black",size = 12))

pdf(paste0(outputDir,"EFig9.pdf"),height=16/2.54, width=36.6/2.54, useDingbats = F, onefile = T)
eFig9
dev.off()


#Monotonic association with RCB - Integrated model in training dataset
m <- l[l$pipeline=="Integrated model" & l$Class=="Training",]
m <- polr(factor(m$RCB.category,ordered = T) ~ m$value,  Hess=TRUE)
ctable <- coef(summary(m))
pval <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
#p = 3e-10 
pval[1]

#Monotonic association with RCB - Integrated model in validation dataset
m <- l[l$pipeline=="Integrated model" & l$Class!="Training",]
m <- polr(factor(m$RCB.category,ordered = T) ~ m$value,  Hess=TRUE)
ctable <- coef(summary(m))
pval <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
#p = 1e-06 
pval[1]

