# Manuscript:   Multi-omic machine learning predictor of breast cancer therapy response
# Author:       Stephen-John Sammut
# Description:  Cohort description
# Section:      Results - Multi-platform profiling of tumour biopsies
#=======================================================================================

rm(list=ls())

library (reshape2)
library (readxl)

# change the contents of variable baseDir to the root analysis folder 
baseDir <- "~/Neoadjuvant-predictor/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/R/LoadDirectoryStructure.R"))

# load clinical metadata (Supplementary Table 1)
metadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 1))

#=========================================================================
# COSORT diagram and clinical characteristics of cohort
#=========================================================================

#168 cases have molecular or digital pathology data
nrow (metadata)

#168 cases have sWGS and WES data
table(metadata$DNA.sequenced)["YES"]
#162 cases have RNA-seq data
table(metadata$RNA.sequenced)["YES"]
#166 cases have digital pathology data
table(metadata$Digital.pathology)["YES"]

#proportion of er/her cases
#39% ER+HER2-, 30% ER+HER2+, 23% ER-HER2- and 8% ER-HER2+
sort(round(prop.table(table(paste(metadata$ER.status,metadata$HER2.status))),2),decreasing = T)

#Chemotherapy regimens delivered
ch <- metadata$NAT.regimen

#regimens containing block sequential taxane/anthracycline
ch[ch %in% c("EC-T","FEC-T","P-EC","P-FEC","T-EC","T-FEC","EC-T + Trastuzumab","FEC-T + Trastuzumab",
             "T-FEC + Trastuzumab","T-FEC + Trastuzumab + Pertuzumab","FEC-T + Trastuzumab + Pertuzumab")]<-"T+A"
# 145 cases received block sequential chemotherapy with a taxane and anthracycline for a median of 6 cycles
length(ch[ch=="T+A"])
median(metadata[ch=="T+A","Chemo.cycles"])


#regimens that contain a taxane, but no anthracycline
table(ch)
ch[ch %in% c("P-Carboplatin","T","T-Carboplatin","TC","P + Trastuzumab","T + Pertuzumab + Trastuzumab",
             "TC + Trastuzumab","TC + Pertuzumab + Trastuzumab")]<-"TaxaneNoAnthra"

#22 cases received a taxane but not an anthracycline for a median of 4 cycles
length(ch[ch=="TaxaneNoAnthra"])
median(metadata[ch=="TaxaneNoAnthra","Chemo.cycles"])
#given in combination with carboplatin in 3 cases and cyclophosphamide in 13 cases
table(metadata$NAT.regimen[ch=="TaxaneNoAnthra"])

#2 patients received one cycle of chemotherapy due to drug toxicities
table(sort(metadata$Chemo.cycles))

# 65 patients with ERBB2-amplified tumours received a median of 3 cycles of anti-HER2 therapy
median(as.numeric(metadata$aHER2.cycles),na.rm = T)

#4 cases had one cycle of anti-HER2 therapy
nrow(metadata[metadata$aHER2.cycles<2,])



#====== subset to cases with RCB assessment

#161 of 168 cases had an RCB assessment
x <- metadata[which(metadata$RCB.category!="NA"),]
nrow(x)

# Upon completion of NAT, 42 (26%) attained pCR, 25 (16%) had minimal RD (RCB-I), 65 (40%) moderate RD (RCB-II) and 29 (18%) extensive RD (RCB-III). 	
table(x$RCB.category)
round(prop.table(table(x$RCB.category)),2)


#retain cases treated with more than one cycle of aHer2 therapy/chemo
x <- x[x$aHER2.cycles>1 &x$Chemo.cycles>1,]

#155 cases with RCB assessment have sWGS and WES data
table(x$DNA.sequenced)["YES"]
#149 cases with RCB assessment have RNA-seq data
table(x$RNA.sequenced)["YES"]
#153 cases with RCB assessment have digital pathology data
table(x$Digital.pathology)["YES"]

#98 HER2+, 57 HER2-
table(x$ER.status,x$HER2.status)

#DNA: 98 HER2- 57 HER2+ 
table(x$DNA.sequenced,x$HER2.status)
table(x[x$HER2.status=="NEG","DNA.sequenced"],x[x$HER2.status=="NEG","ER.status"])

#RNA: 94 HER2- 55 HER2+ 
table(x$RNA.sequenced,x$HER2.status)["YES",]
table(x[x$HER2.status=="NEG","RNA.sequenced"],x[x$HER2.status=="NEG","ER.status"])

#Dig path: 96 HER2- 57 HER2+ 
table(x$Digital.pathology,x$HER2.status)["YES",]
table(x[x$HER2.status=="NEG","Digital.pathology"],x[x$HER2.status=="NEG","ER.status"])


#=========================================================================
# Distribution of RCB score components across classes
#=========================================================================

# Perform analyses with cases that had an RCB assessment and received more than one cycle of therapy
# as detailed in Methods (Statistical testing)
rcb <- metadata[which(metadata$RCB.category!="NA"),]

rcb[,c(15:22)] <- as.numeric(as.character(unlist(rcb[,c(15:22)])))

#calculate RCB components
rcb$dprim <- sqrt(rcb$Tumour.dimension.surgery.1*rcb$Tumour.dimension.surgery.2)
rcb$finv  <- (1-(rcb$Percent.CIS/100))*(rcb$Percent.cellularity/100)
rcb$var.tumour <- 1.4*((rcb$dprim*rcb$finv)^0.17)
rcb$var.ln <- (rcb$Max.LN.met.size* (1-(0.75^rcb$Number.of.positive.LN)) *4)^0.17

rcbm <- melt(rcb,measure.vars = c(30,31,20,21))
rcbm$variable <- factor(rcbm$variable,
                      levels=c("dprim","finv","Number.of.positive.LN","Max.LN.met.size"),
                      labels=c("Tumour bed\narea","Tumour\ncellularity","Number of\npositive LN","Maximum LN\nmetastasis size"))


#Plot distribution of RCB components across RCB categories
EFig_2b_top <- ggplot(rcbm,aes(y=value,x=RCB.category,fill=RCB.category))+
  geom_boxplot(width=0.5, outlier.size = 0.7)+
  facet_wrap(~variable,scales="free",ncol=4)+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  scale_fill_RCB()+
  theme_manuscript(base_size = 12)+
  labs(y="Score",x="RCB category")+
  guides(fill="none")+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 11),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0,0.5,0.5,0.5), "lines"))
print(EFig_2b_top)

df_poly3 <- data.frame( x=c(0, 3.28, 3.28), y=c(3.28, 3.28, 0) )
df_poly2 <- data.frame( x=c(-1, -1, 4.28, 2.36), y=c(2.36, 4.28, -1, -1) )
df_poly1 <- data.frame( x=c(-Inf, -1, 2.36), y=c(-Inf, 2.37, -1) )

#Plot relationship between primary tumour score, lymph node score and RCB category
EFig_2b_bottom <- 
  ggplot()+
  geom_polygon(data=df_poly1, aes(x, y), fill="#ffe671", alpha=0.3) +
  geom_polygon(data=df_poly2, aes(x, y), fill="#fdb462", alpha=0.3) +
  geom_polygon(data=df_poly3, aes(x, y), fill="#ef3b2c", alpha=0.3) +
  geom_point(data=rcb[rcb$pCR.RD!="pCR",],aes(x=var.tumour,y=var.ln,fill=RCB.category),color="black",shape=21,size=1.5)+
  geom_abline(slope=-1, intercept=1.36, linetype="dashed",color="gray20")+
  geom_abline(slope=-1, intercept=3.28, linetype="dashed",color="gray20")+
  annotate(geom = "text",x = 0.4,y=0.5,label="RCB-I",color="#ffe671", fontface="bold")+
  annotate(geom = "text",x = 0.4,y=2.2,label="RCB-II",color="#fdb462", fontface="bold")+
  annotate(geom = "text",x = 2,y=2.2,label="RCB-III",color="#ef3b2c", fontface="bold")+
  coord_cartesian(xlim=c(0,2.6),ylim=c(0,2.6))+
  scale_x_continuous(expand = c(0,0.05))+
  scale_y_continuous(expand = c(0,0.05))+
  theme_manuscript(base_size = 12)+
  labs(x="Primary tumour score",y="Lymph node score")+
  scale_fill_manual(name="RCB category",values=c("#ffe671","#fdb462","#ef3b2c"))+
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 11),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(1,0.5,0.5,2), "lines"))+
  guides(fill="none")

EFig_2b_bottom <- ggarrange(ggplot()+theme_void(),EFig_2b_bottom,ggplot()+theme_void(),nrow=1,widths=c(0.2,0.8,0.25))
print(EFig_2b_bottom)

#Print extended figure 2b
pdf(paste0(outputDir,"EFig2b.pdf"), height=5,width=7,useDingbats = F,onefile = T)
ggarrange(EFig_2b_top,EFig_2b_bottom, nrow=2, heights=c(0.95,1))
dev.off()

