# Manuscript:   Multi-omic machine learning predictor of breast cancer therapy response
# Author:       Stephen-John Sammut
# Description:  Genomic features associated with response
# Section:      Results - Genomic landscapes associate with response
#=======================================================================================

rm (list=ls())

library (ggplot2)
library (MASS)
library (readxl)
library (reshape2)
library (UpSetR)
library (vcd)

# change the contents of variable baseDir to the root analysis folder 
baseDir <- "~/Neoadjuvant-predictor/"

# load directory structure downloaded from github site
source (paste0(baseDir,"/R/LoadDirectoryStructure.R"))

figure_font_size <- 12

# load clinical metadata (Supplementary Table 1)
metadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 1))
# combine ER and HER2 status
metadata$ERHER2.status<-ifelse(metadata$ER.status=="POS","ER+ HER2-","ER- HER2-")
metadata$ERHER2.status<-ifelse(metadata$HER2.status=="POS","HER2+",metadata$ERHER2.status)
metadata$HER2.status<-factor(metadata$HER2.status,levels=c("NEG","POS"),labels=c("HER2-","HER2+"))
metadata$ER.status<-factor(metadata$ER.status,levels=c("NEG","POS"),labels=c("ER-","ER+"))

# Perform analyses with cases that had an RCB assessment and received more than one cycle of therapy
# as detailed in Methods (Statistical testing)
metadata <- metadata[metadata$RCB.category!="NA",]
metadata <- metadata[metadata$Chemo.cycles>1 & metadata$aHER2.cycles>1,]
metadata$RCB.score <- as.numeric(metadata$RCB.score)
metadata$RCB.category <- factor(metadata$RCB.category,labels=c("pCR","RCB-I","RCB-II","RCB-III"), ordered = T)

# load mutation data (Supplementary Table 2)
mutations  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 2))

# load breast cancer driver gene list
driverGenes<- scan(paste0(resourcesDir,"breast-cancer-driver-genes.txt"), what=character(),skip = 1)

#=========================================================================
# Mutation data
#=========================================================================

# 16,134 somatic mutations were identified in the pre-therapy tumour WES data across all cases
nrow(mutations)

# Which driver genes have the highest frequency of coding mutations?
codingMutations <- mutations[!mutations$MAF_Variant %in% 
                               c("3'Flank","3'UTR","5'Flank","5'UTR","IGR","Intron","RNA","Silent",
                                 "Targeted_Region","Translation_Start_Site"),]
m <- codingMutations[,c("Donor.ID","Hugo_Symbol")]
m <- m[!duplicated(m),]
m <- sort(table(m$Hugo_Symbol),decreasing = T)
head(m[names(m) %in% driverGenes],4)
round(head(m[names(m) %in% driverGenes],4)/nrow(metadata),2)


#=========================================================================
# Generate mutation interaction plot
#=========================================================================

# subset coding mutations to those within breast cancer driver genes and
# associate with response
codingMutations <- codingMutations[codingMutations$Donor.ID %in% metadata[metadata$RCB.category!="NA","Donor.ID"],]
codingMutations.driver <-codingMutations[codingMutations$Hugo_Symbol %in% driverGenes,]

m <- as.data.frame.matrix(table(codingMutations.driver$Donor.ID,codingMutations.driver$Hugo_Symbol))
m[m>1] <- 1
#retain driver genes mutated in at least 4 samples
m <- m[,colSums(m)>=4]
m <- merge(m,metadata[metadata$RCB.category!="NA",c("Donor.ID","RCB.category")],by.x=0, by.y="Donor.ID",all.y=T)
m[is.na(m)] <- 0
m$Row.names <- factor(m$Row.names)

m$RCB.category<-factor(m$RCB.category)
colnames(m)[1]<-"Identifier"

table(metadata$RCB.category)

#Generate Extended Data Figure 4a
pdf(paste0(outputDir,"EFig4a_main.pdf"),onefile=FALSE, useDingbats = F, width = 9,height=7)
upset(m, nsets = (ncol(m)-2), point.size = 2,
      mb.ratio =  c(0.5,0.5), nintersects = NA,
      sets.bar.color = "#56B4E9", matrix.color="#474954", shade.color="#7284A8",
      text.scale=1.2, sets.x.label="Number of Cases",
      queries = list(list(query = elements, query.name="pCR", 
                          params = list(c("RCB.category","pCR","RCB-I","RCB-II","RCB-III")), color = "#20A39E",active=T),
                     list(query = elements, query.name="RCB-I", 
                          params = list(c("RCB.category","RCB-I","RCB-II","RCB-III")), color = "#ffe671",active=T),
                     list(query = elements, query.name="RCB-II", 
                          params = list(c("RCB.category","RCB-II","RCB-III")), color = "#fdb462",active=T),
                     list(query = elements, query.name="RCB-III", 
                          params = list(c("RCB.category","RCB-III")), color = "#ef3b2c",active=T)),
      query.legend = "top", order.by = "freq", empty.intersections = NULL)


dev.off()


#=========================================================================
# Somatic mutation associations with response
#=========================================================================

# retain driver genes which are mutated (coding) in at least 25 samples
mm <- m[,c(1,as.numeric(which(colSums(m[,c(2:c(ncol(m)-1))])>=25))+1,ncol(m))]
mm <- mm[mm$Identifier%in% metadata$Donor.ID,]

# Mutations associated with pCR
rcb<-factor(mm$RCB.category,levels=c("RCB-II","RCB-III","RCB-I","pCR"),labels=c(0,0,0,1))
df.pcr<-character()  
for (g in colnames(mm)[2:(ncol(mm)-1)]){
  gene       <- factor(mm[,g])
  model      <- glm(formula = rcb ~ gene, family = "binomial")
  lrm        <- summary(model)
  pval       <- lrm$coefficients[c(2:nrow(lrm$coefficients)),"Pr(>|z|)"]
  odds.ratio <- exp(lrm$coefficients[c(2:nrow(lrm$coefficients)),"Estimate"])
  conf       <- suppressMessages(confint(model))
  ci.low     <- exp(conf[c(2:nrow(conf)),1])
  ci.hi      <- exp(conf[c(2:nrow(conf)),2])
  md         <- data.frame(pval,odds.ratio,ci.low,ci.hi)
  md$variable<- g
  df.pcr <- rbind(df.pcr,md)
}
df.pcr$class<-"pCR"

# Mutations associated with increasing RCB score
mm  <- merge(mm,metadata[,c("Donor.ID","RCB.score")],by.x = "Identifier",by.y="Donor.ID")
df.rcb<-character()  
for (g in colnames(mm)[2:(ncol(mm)-2)]){
  gene       <- factor(mm[,g])
  model      <- glm(formula = mm$RCB.score ~ gene)
  lrm        <- summary(model)
  pval       <- lrm$coefficients[c(2:nrow(lrm$coefficients)),"Pr(>|t|)"]
  odds.ratio <- exp(lrm$coefficients[c(2:nrow(lrm$coefficients)),"Estimate"])
  conf       <- suppressMessages(confint(model))
  ci.low     <- exp(conf[c(2:nrow(conf)),1])
  ci.hi      <- exp(conf[c(2:nrow(conf)),2])
  md         <- data.frame(pval,odds.ratio,ci.low,ci.hi)
  md$variable<- g
  df.rcb <- rbind(df.rcb,md)
}
df.rcb$class<-"Increasing RCB"

mutAssoc <- rbind(df.pcr,df.rcb)
mutAssoc$col <- "INSIG"
mutAssoc[mutAssoc$odds.ratio > 1 & mutAssoc$pval <= 0.05,"col"] <- "UP"
mutAssoc[mutAssoc$odds.ratio < 1 & mutAssoc$pval <= 0.05,"col"] <- "DOWN"
mutAssoc$variable<-factor(mutAssoc$variable,levels=unique(sort(mutAssoc$variable,decreasing = T)))
mutAssoc$class<-factor(mutAssoc$class,levels=c("pCR","Increasing RCB"))

eFig4a_inset <-
  ggplot(data=mutAssoc, aes(y=variable, x=log(odds.ratio))) + 
  geom_errorbarh(aes(xmax = log(ci.hi), xmin = log(ci.low),color=col), size = .5, height = .2) +
  geom_point(aes(color=col),shape=16,size=3)+
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dotted") +
  labs (x="Log odds ratio", y="") + 
  scale_color_manual(values=c("#375E97","gray","#FB6542"),breaks = c("DOWN","INSIG","UP"))+
  scale_x_continuous(limits=c(-2.2,2.2))+
  guides(color="none",size="none")+
  facet_grid(~class,scales="free_y",space = "free")+
  theme_manuscript(base_size = figure_font_size)+
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line = element_blank(),
        strip.background =element_blank(),
        strip.text = element_text(size = 12),
        axis.text.y = element_text(face="italic"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0,0,0.3,-1), "lines"))
eFig4a_inset

#Print Extended Data Figure 4a inset
pdf(paste0(outputDir,"EFig4a_inset.pdf"),height=4.5/2.54, width=8/2.54, useDingbats = F, onefile = T)
eFig4a_inset
dev.off()

#TP53 mutations were associated with pCR (OR: 2.9, p=0.01,
g <- mutAssoc[mutAssoc$variable=="TP53",c(1,2,3,4,5,6)]
paste0("OR: ",round(g[g$class=="pCR",2],2), " p= ",signif(g[g$class=="pCR",1],2))

#PIK3CA-mutations were associated with higher residual disease (OR: 2.1, p=0.002). 
g <- mutAssoc[mutAssoc$variable=="PIK3CA",c(1,2,3,4,5,6)]
paste0("OR: ",round(g[g$class=="Increasing RCB",2],2),", p=",round(g[g$class=="Increasing RCB",1],5))


#=========================================================================
# Association with TMB
# Figure 2a, Extended Data Figure 4c
#=========================================================================

p <- metadata
p <- p[,c("Donor.ID","pCR.RD","RCB.category","RCB.score","Grade.pre.NAT","HER2.status","ER.status","ERHER2.status")]

# The exome kit we have used covers a total of 45.5 Mb of sequence
tmb <- data.frame(table(mutations$Donor.ID)/45.54094)
tmb <- merge(p,tmb,by.x="Donor.ID",by.y="Var1")
colnames(tmb)[ncol(tmb)] <- "tmb"

# TMB is monotonically associated with RCB class, p=0.004
m      <- polr(tmb$RCB.category ~ tmb$tmb,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
round(as.numeric(pval[1]),4)

# tumours that attained pCR had higher mutational loads than those with RD following NAT 
# pCR: 2.3 mutations/Mb
round(median(tmb[tmb$RCB.category=="pCR","tmb"]),1)
# RD: 1.4 mutations/Mb
round(median(tmb[tmb$RCB.category!="pCR","tmb"]),1)
# p=0.0005, Wilcoxon rank sum test
wilcox.test(tmb[tmb$RCB.category=="pCR","tmb"],tmb[tmb$RCB.category!="pCR","tmb"])$p.value

# linear association observed between TMB and RCB score (p=0.002, linear regression)
summary(lm( tmb$RCB.score~tmb$tmb))

#The association between TMB and response was solely observed in HER2- tumours
#HER2-, p=9e-06
wilcox.test(tmb[tmb$HER2.status=="HER2-" & tmb$RCB.category=="pCR","tmb"], tmb[tmb$HER2.status=="HER2-" & tmb$RCB.category!="pCR","tmb"])$p.value
#HER2+, p=1
wilcox.test(tmb[tmb$HER2.status=="HER2+" & tmb$RCB.category=="pCR","tmb"], tmb[tmb$HER2.status=="HER2+" & tmb$RCB.category!="pCR","tmb"])$p.value

#Figure 2a
fig2a <- 
  ggplot(tmb,aes(x=RCB.category,y=tmb,fill=RCB.category))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  labs(x="RCB class",y="TMB")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,color="#FB6542",label.y.npc = 0.92)+
  scale_fill_RCB()+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  coord_cartesian(ylim=c(0,11.5))+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.x = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,1), "lines"))+guides(fill="none")
fig2a


#Extended Data Figure 4c - stratify by HER2 status
table(tmb$HER2.status,tmb$pCR.RD)
eFig4c <- 
  ggplot(tmb,aes(y=tmb,x=pCR.RD,fill=pCR.RD))+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = F,size=3,color="#FB6542",label.y.npc = 0.92)+
  geom_boxplot(outlier.size = 0.6, width=0.6)+
  labs(x="Response",y="TMB")+
  scale_fill_pCR_RD()+
  theme_manuscript(base_size = figure_font_size)+
  facet_grid(~HER2.status,scales = "free")+
  guides(fill="none")+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        axis.title = element_text(),
        axis.line = element_blank(),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))+
  theme(plot.margin = unit(c(1,1,0.5,1), "lines"))
eFig4c
saveRDS(eFig4c,paste0(outputDir,"EFig4c.Rdata"))

#=========================================================================
# Clonal and subclonal mutation distribution
# Figure 2b
#=========================================================================

q <- as.data.frame.matrix(table(mutations$Donor.ID,mutations$ClonalStatus))
q <- q[,c(1,3)]
q <- merge(p,q,by.x="Donor.ID",by.y=0)
colnames(q)[c(ncol(q)-1, ncol(q))] <- c("Clonal","Subclonal")

#Tumours with RD post NAT had a higher percentage of subclonal mutations 
q$ratio=(q$`Subclonal`/(q$`Clonal`+q$`Subclonal`))*100
#RD: 22.2% of mutations were subclonal
median(q[q$pCR.RD!="pCR","ratio"])
#pCR: 10.7% of mutations were subclonal
median(q[q$pCR.RD=="pCR","ratio"])
#Tumours that failed to attain pCR had a higher percentage of subclonal mutations
#p=0.002 Wilcoxon rank sum test
wilcox.test(q[q$pCR.RD=="pCR","ratio"],q[q$pCR.RD!="pCR","ratio"])

fig2b <- 
  ggplot(q,aes(x=RCB.category,y=ratio,fill=RCB.category))+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,color="#FB6542",label.y.npc = 0.92)+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  labs(x="RCB class",y="%Subclonal muts")+
  scale_fill_RCB()+
  theme_manuscript(base_size = figure_font_size)+
  guides(fill="none")+
  theme(panel.grid.major.x = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,1), "lines"))
fig2b

# % clonal mutations monotonically associated with RCB class
# p=0.02
m      <- polr(q$RCB.category ~ q$ratio,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
round(as.numeric(pval[1]),5)


#=========================================================================
# Expressed neoantigen burden distribution
# Figure 2c and Extended Data Figure 4d
#=========================================================================

nagFile <- paste0(dataDir,"transneo-diagnosis-neoantigens.tsv.gz")
nag     <- read.table(nagFile,sep="\t",stringsAsFactors = F, header = T)
nag$id  <- paste0(nag$Donor.ID,"_",nag$Chromosome,":",nag$Start,"_",nag$Reference,"/",nag$Variant)
nag     <- nag[,c(1:6,13,16,20,43)]
nag     <- nag[!duplicated(nag),]
nag     <- data.frame(table(nag$Trial.ID))
colnames(nag) <- c("Donor.ID","numNA")
nag <- merge(tmb,nag,by="Donor.ID")

# The total predicted neoantigen burden strongly correlated with TMB 
# R=0.76
signif(cor.test(nag$tmb,nag$numNA)$estimate,2)
# p=<2e-16
summary(lm(nag$tmb~nag$numNA))

# Tumours that attained pCR had higher expressed neoantigen burdens at diagnosis than those with RD 
# pCR: 28NAg
median(nag[nag$pCR.RD=="pCR","numNA"])
# RD: 17 NAg
median(nag[nag$pCR.RD!="pCR","numNA"])
# p=0.009
signif(wilcox.test(nag[nag$pCR.RD=="pCR","numNA"],nag[nag$pCR.RD!="pCR","numNA"])$p.value,1)

fig2c <-
  ggplot(nag[nag$RCB.category %in% c("pCR","RCB-III"),],aes(x=log10(numNA+1),fill=RCB.category))+
  geom_density(alpha=0.8)+
  annotate("text", x=0.7,y=1.2,label="RCB-III", color="#ef3b2c",size=2.8, family="Helvetica")+
  annotate("text", x=1.9,y=1.2,label="pCR", color="#20A39E",size=2.8, family="Helvetica")+
  labs(y="Density",x=expression(Log[10]~neoantigens))+
  scale_fill_manual(name="",values=c("#20A39E","#ef3b2c"))+
  guides(fill="none")+
  theme_manuscript(base_size = figure_font_size)+
  theme(plot.margin = unit(c(1,0.5,0.3,1), "lines"))
fig2c

# NAg load is monotonically associated with RCB class
# p=0.03
m      <- polr(nag$RCB.category ~ log10(nag$numNA+1),  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
signif(as.numeric(pval[1]),1)

wilcox.test(log10(nag[nag$RCB.category=="pCR","numNA"]+1),log10(nag[nag$RCB.category=="RCB-III","numNA"]+1))

# As with TMB, the association between NAg and response was observed 
# in HER2- (p=0.004, Wilcoxon rank sum test), but not HER2+ (p=0.54) tumours 
round(wilcox.test(log10(nag[nag$HER2.status=="HER2-" & nag$pCR.RD=="pCR","numNA"]+1),
                  log10(nag[nag$HER2.status=="HER2-" & nag$pCR.RD!="pCR","numNA"]+1))$p.value,4)

round(wilcox.test(log10(nag[nag$HER2.status=="HER2+" & nag$pCR.RD=="pCR","numNA"]+1),
                  log10(nag[nag$HER2.status=="HER2+" & nag$pCR.RD!="pCR","numNA"]+1))$p.value,4)

# Extended Data Figure 4c - NAg association with response by HER2 status
eFig4d <- 
  ggplot(nag,aes(y=numNA,x=pCR.RD,fill=pCR.RD))+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,color="#FB6542",label.y.npc = 0.92)+
  geom_boxplot(outlier.size = 0.6, width=0.6)+
  labs(x="Response",y="Expressed NAg")+
  scale_fill_pCR_RD()+
  theme_manuscript(base_size = figure_font_size)+
  facet_grid(~HER2.status,scales = "free")+
  guides(fill="none")+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        axis.title = element_text(),
        axis.line = element_blank(),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))+
  theme(plot.margin = unit(c(1,1,0.5,1), "lines"))
eFig4d
saveRDS(eFig4d,paste0(outputDir,"EFig4d.Rdata"))


#=========================================================================
# Mutational signature distribution
# Figure 2d and Extended Data Figure 4e
#=========================================================================

mutSigs <- read.table(paste0(dataDir,"transneo-diagnosis-mutational-signatures.tsv.gz"), header=T, stringsAsFactors = F, sep="\t",row.names = 1)
# Select a significance threshold of 5% to 'call' a signature
mutSigs[mutSigs<0.05 & mutSigs!=0] <- 0

# Non-clock contribution to response - Extended Data Figure 4d
numutSigs <- apply(mutSigs,2,function(x) sum(x!=0))
mutSigs   <- mutSigs[,colnames(mutSigs) %in% names(numutSigs[numutSigs>=1])]
mutSigs$Unknown <- NULL
mutSigs$NumberMutationsAnalysed <- NULL

ms <- merge(metadata[,c("Donor.ID","pCR.RD","RCB.category","HER2.status")],mutSigs, by.x="Donor.ID",by.y=0)
efig4e <-
  ggplot(ms,aes(x=RCB.category,y=1-(Signature.1+Signature.5),fill=RCB.category))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  labs(x="RCB class",y="%Non-clock",title=NULL)+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,
                     label.y.npc = 0.92,size=5.5,color="#FB6542")+
  scale_fill_RCB()+
  #facet_grid(~HER2.status)+
  scale_y_continuous(limits=c(0,1.1))+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  theme_manuscript(base_size = figure_font_size)+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank())+
  guides(fill="none")+
  theme(plot.margin = unit(c(2,0.5,0.5,1), "lines"))
efig4e
saveRDS(efig4e,paste0(outputDir,"EFig4e.Rdata"))

# Non-clock contribution is monotonically associated with RCB class, p=0.005
y      <- 1-(ms$Signature.1+ms$Signature.5)
m      <- polr(ms$RCB.category ~ y,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
round(as.numeric(pval[1]),5)

# Tumours that attained pCR had a lower contribution of non-clock signatures, p=0.002
wilcox.test(y[ms$RCB.category=="pCR"],y[ms$RCB.category!="pCR"])$p.value

# Associations between mutational signature and response

#retain signatures present in at least 10 samples
m <- mutSigs[,colnames(mutSigs) %in% names(numutSigs[numutSigs>=10])]
#normalise to signature 1
e <- (t(apply(m,1, function(x) log2((x/x[1])+0.001))))
e[is.infinite(e) | is.nan(e)] <- NA
e <- as.data.frame(e)
#remove signature one as we normalised to it
e[,1]<-NULL
m<-merge(metadata[,c("Donor.ID","RCB.category","pCR.RD","ER.status","HER2.status")],e, by.x="Donor.ID",by.y=0,all.x=T)

rcb<-factor(m$pCR.RD,levels=c("RD","pCR"),labels=c(0,1))
er<-factor(m$ER.status)
her2<-factor(m$HER2.status)
signatures<-grep("Signature",colnames(m),value = T)
sigMatrix<-character()
for (s in signatures){
  sig   <- m[,s]
  model <- glm(formula = rcb ~ er+her2+sig, family = binomial)
  lrm   <- summary(model)
  if (dim(lrm$coefficients)[1]>1){
    pval       <- lrm$coefficients[4,"Pr(>|z|)"]
    odds.ratio <- exp(lrm$coefficients[4,"Estimate"])
    conf       <- suppressMessages(confint(model))
    ci.low     <- exp(conf[c(4:nrow(conf)),1])
    ci.hi      <- exp(conf[c(4:nrow(conf)),2])
    md <- data.frame(pval,odds.ratio,ci.low,ci.hi)
    md$variable <- s
    sigMatrix <- rbind(sigMatrix,md)
  }
}

sigMatrix <- sigMatrix[complete.cases(sigMatrix),]
sigMatrix$col <- "INSIG"
sigMatrix[sigMatrix$odds.ratio > 1 & sigMatrix$pval <= 0.05,"col"] <- "UP"
sigMatrix[sigMatrix$odds.ratio < 1 & sigMatrix$pval <= 0.05,"col"] <- "DOWN"
sigMatrix$variable <- gsub("Signature\\.","",sigMatrix$variable)
sigMatrix$variable <- factor(sigMatrix$variable, levels=unique(sigMatrix$variable))

fig2d <-
  ggplot(data=sigMatrix[!sigMatrix$variable %in% c(15),],aes(y=reorder(variable,(odds.ratio)), x= (odds.ratio))) + 
  geom_errorbarh(aes(xmax = (ci.hi), xmin = (ci.low),color=col), size = .5, height = .2) +
  geom_point(aes(color=col),shape=16,size=2.5)+
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  labs (x="Odds ratio (pCR)", y="Mutational signature") +
  scale_color_manual(values=c("gray80","#FB6542","#375E97"))+
  scale_x_continuous(breaks=c(0.8,1,1.2,1.4))+
  guides(color="none", size="none")+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major = element_blank(),
        plot.margin = unit(c(1,0.5,0.5,1), "lines"))
fig2d


#=========================================================================
# Associations with HRD
# Figure 2e, Extended Data Figure 2f
#=========================================================================

hrd <- read.table(paste0(dataDir,"transneo-diagnosis-HRD.tsv.gz"),header = T, stringsAsFactors = F)
hrd <- merge(metadata[,c("Donor.ID","HER2.status","RCB.category","pCR.RD","ER.status")],hrd,by.x="Donor.ID",by.y="Trial.ID")

# HRD is monotonically associated with RCB class
# p=0.00001
m      <- polr(hrd$RCB.category ~ hrd$HRD.sum,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
signif(as.numeric(pval[1]),1)

hrd <- melt(hrd, measure.vars = c("HRD", "Telomeric.AI", "LST", "HRD.sum"))
hrd$variable <- factor(hrd$variable,levels=c("HRD","Telomeric.AI","LST","HRD.sum"),
                       labels=c("HRD-LOH", "Telomeric AI", "LST","HRD score"))

eFig4f <-
  ggplot(hrd[hrd$variable=="HRD score",],aes(x=pCR.RD,y=value,fill=pCR.RD))+
  geom_boxplot(outlier.size = 0.6, width=0.6)+
  scale_fill_pCR_RD()+
  facet_grid(~HER2.status, scales="free")+
  labs(x="Response",y="HRD score")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,
                     size=5.5,color="#FB6542",label.y.npc = 0.92)+
  theme_manuscript(base_size = figure_font_size)+
  guides(fill="none")+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 11),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(1,0.5,0.5,1), "lines"))
eFig4f
saveRDS(eFig4e,paste0(outputDir,"EFig4e.Rdata"))

fig2e <-
  ggplot(hrd[hrd$variable=="HRD score",],aes(x=RCB.category,y=value,fill=RCB.category))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  scale_fill_RCB()+
  labs(x="RCB class",y="HRD Score")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,
                     size=5.5,color="#FB6542",label.y.npc = 0.92)+
  theme_manuscript(base_size = figure_font_size)+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  theme_manuscript(base_size = figure_font_size)+
  guides(fill="none")+
  theme(panel.grid.major.x = element_blank(),
        plot.margin = unit(c(1,0.5,0.5,1), "lines"))
fig2e


# Increasing HRD is associated with pCR in the HER2- cohort, p=3e-06                                                
wilcox.test((hrd[hrd$HER2.status!="HER2+" & hrd$pCR.RD=="pCR" & hrd$variable=="HRD score","value"]),
            (hrd[hrd$HER2.status!="HER2+" & hrd$pCR.RD!="pCR" & hrd$variable=="HRD score","value"]))
wilcox.test((hrd[hrd$HER2.status=="HER2+" & hrd$pCR.RD=="pCR" & hrd$variable=="HRD score","value"]),
            (hrd[hrd$HER2.status=="HER2+" & hrd$pCR.RD!="pCR"& hrd$variable=="HRD score","value"]))


#=========================================================================
# Associations with global CNA
# Figure 2f and Extended Data Figure 4g
#=========================================================================

cin <- read.table(paste0(dataDir,"transneo-diagnosis-ASCAT-CIN.tsv.gz"), sep="\t",header = T, stringsAsFactors = F, quote = "")
pathologyCin <- merge(p,cin, by.x="Donor.ID",by.y="Trial.ID",all.x=T)
pathologyCin$cin <- pathologyCin$cin*100

fig2f <-
  ggplot(pathologyCin,aes(x=RCB.category,y=cin,fill=RCB.category))+
  geom_boxplot(outlier.size = 0.6, width=0.7)+
  labs(x="RCB class",y="%CNA")+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,label.y.npc = 0.92,color="#FB6542")+
  scale_fill_RCB()+
  scale_y_continuous(limits = c(0,90))+
  scale_x_discrete(breaks=c("pCR","RCB-I","RCB-II","RCB-III"),labels=c("pCR","I","II","III"))+
  theme_manuscript(base_size = figure_font_size)+
  guides(fill="none")+
  theme(panel.grid.major.x = element_blank(),
        plot.margin = unit(c(1,0.5,0.3,1), "lines"))
fig2f


# %CIN is monotonically associated with RCB class
# p=0.0002
m      <- polr(pathologyCin$RCB.category ~ pathologyCin$cin,  Hess=TRUE)
ctable <- coef(summary(m))
pval   <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
round(as.numeric(pval[1]),5)

eFig4g <- 
  ggplot(pathologyCin,aes(x=pCR.RD,y=cin,fill=pCR.RD))+
  stat_compare_means(ref.group="pCR",label = "p.format", hide.ns = T,size=5.5,color="#FB6542",label.y.npc = 0.92)+
  geom_boxplot(outlier.size = 0.6, width=0.6)+
  labs(x="Response",y="%CNA")+
  scale_fill_pCR_RD()+
  theme_manuscript(base_size = figure_font_size)+
  facet_grid(~HER2.status,scales = "free")+
  guides(fill="none")+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 11),
        axis.title = element_text(),
        axis.line = element_blank(),
        strip.background =element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(1,1,0.5,1), "lines"))
eFig4g
saveRDS(eFig4g,paste0(outputDir,"EFig4f.Rdata"))


#=========================================================================
# Associations with HLA LOH, Global LOH and CNA burden
# Extended Data Figure 4i
#=========================================================================

hlaloh <- read.table(paste0(dataDir,"transneo-diagnosis-lohhla.tsv.gz"),header = T,stringsAsFactors = F)
hlaloh$LOH.HLA.Allele <- paste0(toupper(sapply(strsplit(hlaloh$LOH.HLA.Allele,"_"),"[",1)),"-",
                                toupper(sapply(strsplit(hlaloh$LOH.HLA.Allele,"_"),"[",2)),"*",
                                sapply(strsplit(hlaloh$LOH.HLA.Allele,"_"),"[",3),":",
                                sapply(strsplit(hlaloh$LOH.HLA.Allele,"_"),"[",4))
hlaloh$ID <- paste0(hlaloh$Trial.ID,"_",hlaloh$LOH.HLA.Allele)

# 29 cases have HLA-LOH
length(unique(hlaloh$Trial.ID))

# A median of 30% of NAg were not presented due to LOH HLA
nag   <- read.table(nagFile,sep="\t",stringsAsFactors = F, header = T)
nag$id<- paste0(nag$Trial.ID,"_",nag$Chromosome,":",nag$Start,"_",nag$Reference,"/",nag$Variant)
nag   <- nag[,c(1:6,13,16,20,43)]
nag   <- nag[!duplicated(nag),]
nag   <- nag[nag$Trial.ID %in% hlaloh$Trial.ID,]
nag$NagLost<-0
nag[paste0(nag$Trial.ID,"_",nag$HLA.Allele) %in% hlaloh$ID,"NagLost"]<-1
zz <- table(nag$Trial.ID,nag$NagLost)
median(zz[,2]/rowSums(zz))

# 69% of LOH events occurring within HLA molecules that presented an equal or greater 
# number of neoepitopes than its retained alternative allele
lost<-numeric()
for (a in c(1:nrow(hlaloh))){
  hlatype<-nag[nag$Trial.ID %in% hlaloh[a,"Trial.ID"] & grepl( substr(hlaloh[a,"LOH.HLA.Allele"],1,5), nag$HLA.Allele),]
  lost<-rbind(lost,table(hlatype$NagLost))
}
sum(lost[,2]>=lost[,1])/nrow(lost)

loh.p<-metadata
loh.p$LOH<-"No"
loh.p[loh.p$Donor.ID %in% hlaloh$Trial.ID,"LOH"]<-"Yes"
loh.p$pCR.RD<-factor(loh.p$pCR.RD,levels=c("RD","pCR"))
loh.p$LOH <- factor(loh.p$LOH, levels=c("No","Yes"))
loh.p$ER.status <- factor(loh.p$ER.status,levels=c("ER+","ER-"))
loh.p$HER2.status <- factor(loh.p$HER2.status,levels=c("HER2+","HER2-"))

# is HLA-LOH associated with pCR?
model <- glm(formula = pCR.RD ~ ER.status+HER2.status+LOH, data=loh.p,family = "binomial")
lrm   <- summary(model)
pval       <- lrm$coefficients[c(2:nrow(lrm$coefficients)),"Pr(>|z|)"]
odds.ratio <- exp(lrm$coefficients[c(2:nrow(lrm$coefficients)),"Estimate"])
conf       <- suppressMessages(confint(model))
ci.low     <- exp(conf[c(2:nrow(conf)),1])
ci.hi      <- exp(conf[c(2:nrow(conf)),2])
md.hlaloh  <- data.frame(pval,odds.ratio,ci.low,ci.hi,type="HLA LOH", variables=names(ci.hi))

# is global LOH associated with pCR?
segments <- read.table(paste0(dataDir,"transneo-diagnosis-ASCAT-segments.tsv.gz"),header = T,stringsAsFactors = F)
purity   <- read.table(paste0(dataDir,"transneo-diagnosis-ASCAT-purity.tsv.gz"),header = T,stringsAsFactors = F)
segments <- segments[!segments$Trial.ID %in% as.character(purity[purity$purity==1,1]),]
segments$length <- (segments$endpos - segments$startpos)+1
loh <- segments[segments$nMinor==0,]
loh <- aggregate(loh$length, by=list(Category=loh$Trial.ID), FUN=sum)
allgenome <- aggregate(segments$length, by=list(Category=segments$Trial.ID), FUN=sum)
loh <- data.frame(Donor.ID=loh$Category, loh=((loh$x/allgenome$x)*100))
globalLoh <- merge(loh.p,loh,by="Donor.ID")
model <- glm(formula = pCR.RD ~ ER.status+HER2.status+loh, data=globalLoh,family = "binomial")
lrm <- summary(model)
pval       <- lrm$coefficients[c(2:nrow(lrm$coefficients)),"Pr(>|z|)"]
odds.ratio <- exp(lrm$coefficients[c(2:nrow(lrm$coefficients)),"Estimate"])
conf       <- suppressMessages(confint(model))
ci.low     <- exp(conf[c(2:nrow(conf)),1])
ci.hi      <- exp(conf[c(2:nrow(conf)),2])
md.globalloh <- data.frame(pval,odds.ratio,ci.low,ci.hi,type="Global LOH", variables=names(ci.hi))

# is global CNA associated with pCR?
globalCIN <- merge(cin,loh.p,by.x="Trial.ID", by.y="Donor.ID")
model <- glm(formula = pCR.RD ~ ER.status+HER2.status+cin, data=globalCIN,family = "binomial")
lrm   <- summary(model)
pval       <- lrm$coefficients[c(2:nrow(lrm$coefficients)),"Pr(>|z|)"]
odds.ratio <- exp(lrm$coefficients[c(2:nrow(lrm$coefficients)),"Estimate"])
conf       <- suppressMessages(confint(model))
ci.low     <- exp(conf[c(2:nrow(conf)),1])
ci.hi      <- exp(conf[c(2:nrow(conf)),2])
md.globalcin<-data.frame(pval,odds.ratio,ci.low,ci.hi,type="CNA Burden", variables=names(ci.hi))

# combine all ORs
md<-rbind(md.hlaloh,md.globalloh,md.globalcin)
md$col="INSIG"
md[md$odds.ratio > 1 & md$pval <= 0.05,"col"] <- "UP"
md[md$odds.ratio < 1 & md$pval <= 0.05,"col"] <- "DOWN"
md$variable<-factor(md$variable,levels=c("ER.statusER-","HER2.statusHER2-","LOHYes","loh","cin"),labels=c("ER-","HER2-","HLA LOH","%LOH","%CNA"))


eFig4i <-
  ggplot(data=md,  aes(y=variable, x= log(odds.ratio))) + 
  geom_errorbarh(aes(xmax = log(ci.hi), xmin = log(ci.low),color=col), size = .5, height = .2) +
  geom_point(aes(color=col),shape=16, size=3)+
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dotted") +
  labs (x="Log odds ratio (pCR)", y="") + 
  scale_color_manual(values=c("#375E97","#FB6542"))+
  guides(color="none")+
  facet_wrap(~type,scales="free_y")+
  theme_manuscript(base_size = figure_font_size)+
  guides(size="none")+
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.line = element_blank(),
        strip.background =element_blank(),
        strip.text = element_text(size = figure_font_size-1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.1,0.5,0), "lines"))
eFig4i
saveRDS(eFig4i,paste0(outputDir,"EFig4i.Rdata"))


#=========================================================================
# iC10 subtypes associated with response
# Extended Data Figure 4h
#=========================================================================

iC10 <- metadata
iC10$iC10 <- factor(iC10$iC10,levels=(c(1:10)))

d<-table(`iC10 classification`=iC10$iC10,`RCB class`=iC10$RCB.category)
e=chisq.test(d)

pdf(paste0(outputDir,"EFig4h.pdf"), onefile=T, useDingbats = F,height=5, width=5)
mosaic(d, shade=TRUE, legend=T, 
       gp = shading_hsv,gp_args = list(interpolate = 1:6),
       labeling_args=list(gp_labels=gpar(fontsize=9,fontfamily = "sans"),
                          gp_varnames=gpar(fontsize=10,fontfamily = "sans")), 
       legend_args=list(fontsize=9),
       margins=c(3,3,3,3))
dev.off()

