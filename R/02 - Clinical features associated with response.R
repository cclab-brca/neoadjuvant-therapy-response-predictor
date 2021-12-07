# Manuscript:   Multi-omic machine learning predictor of breast cancer therapy response
# Author:       Stephen-John Sammut
# Description:  Associations between clinical features and response
# Section:      Results - Clinical phenotypes are limited predictors
#=======================================================================================

rm(list=ls())

library (readxl)
library (reshape2)

# change the contents of variable baseDir to the root analysis folder 
baseDir <- "~/Neoadjuvant-predictor/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/R/LoadDirectoryStructure.R"))

# load clinical metadata (Supplementary Table 1)
metadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 1))

#161 of 168 cases had an RCB assessment
metadata <- metadata[which(metadata$RCB.category!="NA"),]
nrow(metadata) 

#155 cases had more than 1 cycle of therapy
metadata <- metadata[which(as.numeric(metadata$Chemo.cycles)>1),]
metadata$aHER2.cycles=as.numeric(metadata$aHER2.cycles)
metadata <- metadata[c(which(is.na(metadata$aHER2.cycles)),which(metadata$aHER2.cycles>1)),]
nrow(metadata)

#=========================================================================
# Association of clinical variables with response
#=========================================================================

y <- metadata
rownames(y) <- y$Donor.ID
clinVariables <- c("T.stage","Histology","ER.status","HER2.status","LN.status.at.diagnosis","Grade.pre.NAT","RCB.category")
y <- y[,colnames(y) %in% clinVariables]

y$RCB.category  <- factor(y$RCB.category,levels=c("RCB-III","RCB-II","RCB-I","pCR"),labels=c(0,0,0,1))
y$T.stage       <- as.numeric(substr(y$T.stage,2,2))
y$Histology     <- factor(y$Histology=="IDC")
y$ER.status     <- factor(y$ER.status,levels=c("POS","NEG"))
y$HER2.status   <- factor(y$HER2.status,levels=c("NEG","POS"))
y$Grade.pre.NAT <- factor(y$Grade.pre.NAT,levels=c(2,3))
y$LN.status.at.diagnosis <- factor(y$LN.status.at.diagnosis,levels=c("POS","NEG"))

#simple logistic regression
simpleLogisticRegression <- lapply( y[,-ncol(y)], function(t) {
  model      <- glm( y$RCB.category ~ t, family="binomial")
  lrm        <- summary(model)
  pval       <- lrm$coefficients[c(2:nrow(lrm$coefficients)),"Pr(>|z|)"]
  odds.ratio <- exp(lrm$coefficients[c(2:nrow(lrm$coefficients)),"Estimate"])
  conf       <- suppressMessages(confint(model))
  ci.low     <- exp(conf[c(2:nrow(conf)),1])
  ci.hi      <- exp(conf[c(2:nrow(conf)),2])
  data.frame(pval,odds.ratio,ci.low,ci.hi)
} )
regression.univariable <- do.call(rbind,simpleLogisticRegression)
regression.univariable$variable <- names(simpleLogisticRegression)
regression.univariable$class <- "pCR"
regression.univariable$type <- "Simple logist."

#multiple logistic regression
model      <- glm(RCB.category~., family="binomial",data=y)
lrm        <- summary(model)
pval       <- lrm$coefficients[c(2:nrow(lrm$coefficients)),"Pr(>|z|)"]
odds.ratio <- exp(lrm$coefficients[c(2:nrow(lrm$coefficients)),"Estimate"])
conf       <- suppressMessages(confint(model))
ci.low     <- exp(conf[c(2:nrow(conf)),1])
ci.hi      <- exp(conf[c(2:nrow(conf)),2])
regression.multivariable <- data.frame(pval,odds.ratio,ci.low,ci.hi)
regression.multivariable$variable <- rownames(regression.multivariable)
regression.multivariable$class <- "pCR"
regression.multivariable$type <- "Multiple logist."
regression.univariable$variable <- regression.multivariable$variable

lr <- rbind(regression.univariable,regression.multivariable)
lr$type <- factor(lr$type, levels=c("Simple logist.","Multiple logist."))
lr <- lr[complete.cases(lr),]

#adjust p value
lr$pval <- p.adjust(lr$pval,method = "BH")
lr$col <- "INSIG"
lr[lr$odds.ratio > 1 & lr$pval < 0.05,"col"] <- "UP"
lr[lr$odds.ratio < 1 & lr$pval < 0.05,"col"] <- "DOWN"

lr$variable<-factor(lr$variable,
                    levels = rev(c("T.stage","LN.status.at.diagnosisNEG","Grade.pre.NAT3",
                                   "HistologyTRUE","ER.statusNEG","HER2.statusPOS")),
                    labels = rev(c("Tumour size","LN- vs +","Grade 3 vs 2",
                                   "IDC histology", "ER- vs +",
                                   "HER2+ vs -")))

EFig2c <-
  ggplot(data=lr,aes(y=variable, x= (odds.ratio))) + 
  geom_errorbarh(aes(xmax = (ci.hi), xmin = (ci.low),color=col), size = .5, height = .2) +
  geom_point(aes(color=col),size=2)+
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dotted") +
  labs (x="Odds ratio", y="") +
  scale_color_manual(values=c("#375E97","gray70","#FB6542"),breaks=c("DOWN","INSIG","UP"))+
  guides(color="none",size="none")+
  facet_grid(~type)+
  theme_manuscript()+
  coord_cartesian(xlim=c(0,15))+
  theme(legend.position = "bottom",
        axis.line = element_blank(),
        strip.text = element_text(size = 11),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0,1,0.5,1), "lines"))
EFig2c


#Print Extended Figure 2c
pdf(paste0(outputDir,"EFig2c.pdf"), height=2,width=4.5,useDingbats = F,onefile = T)
print(EFig2c)
dev.off()

#print out significance results
for (test in c("Simple logist.","Multiple logist.")){
  g <- lr[lr$type==test,]
  g <- g[,c(2,1,5,3,4)]
  g <- g[g$pval<0.05,]
  g <- g[order(g$odds.ratio,decreasing = T),]
  g$label <- paste0("(OR: ",round(g$odds.ratio,1),", CI: ",round(g$ci.low,1),"-",round(g$ci.hi,1),", FDR=",signif(g$pval,1),")")
  g <- g[,c(3,6),drop=F]
  cat(test,",pCR\n")
  for (z in 1:nrow(g)){
    cat(paste0(g[z,1]," ",g[z,2], ", "))
  }
  cat("\n")
}

#=========================================================================
# Generate bar plot summarising clinical characteristics of cohort
#=========================================================================
x <- metadata
x$Histology   <- ifelse(grepl("IDC$",x$Histology),"NST","Oth")
x$ER.status   <- as.character(factor(x$ER.status, levels=c("POS","NEG"),labels=c("ER+","ER-")))
x$HER2.status <- as.character(factor(x$HER2.status, levels=c("POS","NEG"),labels=c("HER2+","HER2-")))
x$LN.status.at.diagnosis <- ifelse(grepl("POS",x$LN.status.at.diagnosis),"N1+","N0")

#classify NAT regimens
x[grep("Trastuzumab",x$NAT.regimen), "NAT.regimen"]<-"aHER2"
x[x$NAT.regimen %in% c("EC-T","FEC-T","P-EC","P-FEC","T-EC","T-FEC"),"NAT.regimen"]<-"T+A"
x[x$NAT.regimen %in% c("P-Carboplatin","T","T-Carboplatin","TC"), "NAT.regimen"]<-"T"

#compute ypT and ypN stage
x$Tumour.dimension.surgery.1 <- as.numeric(x$Tumour.dimension.surgery.1)
x$ypTstage<-"Not assessed"
x[which(x$Tumour.dimension.surgery.1==0),"ypTstage"]<-"ypT0"
x[which(x$Tumour.dimension.surgery.1<=20 & x$ypTstage=="Not assessed"),"ypTstage"]<-"ypT1"
x[which(x$Tumour.dimension.surgery.1<=50 & x$ypTstage=="Not assessed"),"ypTstage"]<-"ypT2"
x[which(x$Tumour.dimension.surgery.1>50 & x$ypTstage=="Not assessed"),"ypTstage"]<-"ypT3"
x$ypTstage<-factor(x$ypTstage,levels=c("ypT0","ypT1","ypT2","ypT3"),
                   labels=c("ypT0","ypT1","ypT2","ypT3"))

x$Number.of.positive.LN<-as.numeric(x$Number.of.positive.LN)
x$ypNstage<-"Not assessed"
x[which(as.numeric(x$Number.of.positive.LN)==0),"ypNstage"]<-"ypN0"
x[x$LN.met.size<0.2,"ypNstage"]<-"ypN0"
x[which(x$Number.of.positive.LN<=3 & x$ypNstage=="Not assessed"),"ypNstage"]<-"ypN1"
x[which(x$Number.of.positive.LN<=9 & x$ypNstage=="Not assessed"),"ypNstage"]<-"ypN2"
x[which(x$Number.of.positive.LN>9 & x$ypNstage=="Not assessed"),"ypNstage"]<-"ypN3"
x$ypNstage<-factor(x$ypNstage,levels=c("ypN0","ypN1","ypN2","ypN3"),
                   labels=c("ypN0","ypN1","ypN2","ypN3"))

x$PAM50<-factor(x$PAM50,levels=c("LumA","LumB","Basal","Her2","Normal","Unk"),
                labels=c("A","B","Ba","H","N","U"))

x$Surgery.type<-factor(x$Surgery.type,levels=c("Mastectomy","WLE"),labels=c("Mx","WLE"))
x$LVI<-factor(x$LVI,levels=c("YES","NO","N/A"),labels=c("LVI+","LVI-","N/A"))

e<-x[,c("Donor.ID","Grade.pre.NAT","Histology", "ER.status","HER2.status","LVI",
        "RCB.category","LN.status.at.diagnosis","T.stage","NAT.regimen","PAM50",
        "Surgery.type","ypTstage","ypNstage")]

df<-character()
for (g in colnames(e)[2:ncol(e)]){
  df<-rbind(df,cbind(melt(prop.table(table(e$RCB.category,e[,g]),margin = 1)),g))
}
df<-df[complete.cases(df),]
colnames(df)<-c("RCB","value","prop","variable")
df<-df[df$prop!=0,]
df<-df[df$variable!="RCB.category",]

df$value<-factor(df$value,levels = rev(c("T1","T2","T3","T4",
                                         "N0","N1+",
                                         "ER-","ER+",
                                         "HER2-","HER2+",
                                         "2","3",
                                         "NST","Oth",
                                         "T+A","T", "aHER2",
                                         "Mx","WLE",
                                         "ypT0","ypT1","ypT2","ypT3",
                                         "ypN0","ypN1","ypN2","ypN3",
                                         "LVI-","LVI+",
                                         "A","B","Ba","H","N","U")))

df$type <- "ALL"
df[df$variable %in% c("T.stage","LN.status.at.diagnosis"),"type"]<-"Clinical"
df[df$variable %in% c("ER.status","HER2.status","Grade.pre.NAT","Histology"),"type"]<-"Histological"
df[df$variable %in% c("NAT.regimen","Surgery.type"),"type"]<-"Treatment"
df[df$variable %in% c("ypTstage","ypNstage","LVI"),"type"]<-"Post-Treatment"
df[df$variable %in% c("PAM50"),"type"]<-"Biomark"
df$type<-factor(df$type,levels=c("Clinical","Histological","Treatment","Post-Treatment","Biomark"))
df<-df[rev(order(df$value)),]

cols.tstage    <- c('#deebf7','#c6dbef','#9ecae1','#6baed6')
cols.lymphNode <- c("#c6dbef","#6baed6")
cols.er        <- c("#bae4b3","#74c476")  
cols.her2      <- c("#bae4b3","#74c476")  
cols.grade     <- c("#bae4b3","#74c476")  
cols.histo     <- c("#bae4b3","#74c476")  
cols.nac       <- c('#fcbba1','#fc9272','#fb6a4a')
cols.surgery   <- c('#fcbba1','#fb6a4a')
cols.pam50     <- c('#fffcd5','#fff9c9','#fff3b3','#ffed9c','#ffe785','#ffffe0')
cols.postRxT   <- c("#efe7f7","#d2bbe8","#ab80d5","#8d54c6")
cols.postRxN   <- c("#efe7f7","#d2bbe8","#ab80d5","#8d54c6")
cols.lvi       <- c("#efe7f7","#8d54c6")  
cols<-rev( c(cols.tstage,cols.lymphNode,
             cols.er,cols.her2,cols.grade,cols.histo,
             cols.nac,cols.surgery,cols.postRxT,cols.postRxN,cols.lvi,
             cols.pam50))


df$variable<-factor(df$variable,levels=rev(c(
  "T.stage","LN.status.at.diagnosis",
  "ER.status","HER2.status","Grade.pre.NAT","Histology",
  "NAT.regimen","Surgery.type","ypTstage","ypNstage","LVI","PAM50")))

names(cols)<-(levels(df$value))

df$label <- as.character(df$value)
df$label <- gsub("^yp|^ic","",df$label)
df$label <- gsub("^T0|^N0","0",df$label)
df$label <- gsub("^T1|^N1","1",df$label)
df$label <- gsub("^T2|^N2","2",df$label)
df$label <- gsub("^T3|^N3","3",df$label)
df$label <- gsub("^T4|^N4","4",df$label)
df$label <- gsub("^ER\\+|^HER2\\+|^LVI\\+","\\+",df$label)
df$label <- gsub("^ER\\-|^HER2\\-|^LVI\\-","\\-",df$label)
df$label <- gsub("^ST","NST",df$label)
df[df$label=="Oth" &df$variable=="Histology" &df$RCB=="pCR","label"] <- ""
df[df$label=="Oth" &df$variable=="Histology" &df$RCB=="RCB-I","label"] <- ""

baseSize=13
EFig2d <- ggplot(df,aes(y=prop,x=variable))+
  geom_bar(stat="identity", position="fill",aes(fill=value),color="gray30",size=0.2)+
  scale_fill_manual(values = cols)+
  geom_text(aes(label = label), fontface="bold",color="black",size = 3, 
            position = position_stack(vjust = 0.5)) +
  facet_grid(type~RCB, scales="free", space = "free")+
  labs(y="Proportion of cases")+
  scale_x_discrete(expand = c(0, 0), breaks= levels(df$variable),
                   labels=rev(c("T stage","N stage",
                                "ER","HER2","Grade","Histology",
                                "NAT","Surgery",
                                "ypT stage","ypN stage","LVI","PAM50")))+
  scale_y_continuous(breaks=c(0,0.5,1),labels = c("0%","50%","100%"))+
  theme_grey(base_size = baseSize) + 
  theme(
    axis.text=element_text(size = baseSize*0.7, face="bold",colour = "black"),
    axis.title=element_blank(),
    axis.ticks=element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    panel.background = element_blank(),
    strip.text.x = element_text(size = baseSize*0.8, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.00001, "lines"),
    panel.spacing.y = unit(0.8, "lines"))+
  coord_flip()

EFig2d

#Print Extended Figure 2d
pdf(paste0(outputDir,"EFig2d.pdf"), height=5, width=11, useDingbats = F, onefile = T)
EFig2d
dev.off()
