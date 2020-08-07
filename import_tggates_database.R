#REMOVE UNDERSCORE FROM FILE NAMES

#HEPATOCYTE HYPERTROPHY
setwd("S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose\\HHyp")
## Getting the tsv file names
csv_files1 <- dir(path= 'S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose\\HHyp', pattern='*.tsv$', recursive = T)
## Merging the tsv files
for(i in 1:length(csv_files1)) {
  if(i == 1)
    df1 <- read.csv(csv_files1[i], sep = '\t')
  else
    df1 <- rbind(df1, read.csv(csv_files1[i], sep = '\t'))
}
## Getting subset of df (24 hr, No 'No Chip', High or Control)
df_24hr_nochip_high_control_hep = subset(df1, (df1$SACRI_PERIOD == '24 hr'& df1$BARCODE != 'No ChipData') & (df1$DOSE_LEVEL == 'High' | df1$DOSE_LEVEL == 'Control') )
## BARCODE Retrieval
files1 <- matrix(data=NA, nrow = 5, ncol= 6)
count <- 1
for( i in 1:length(df_24hr_nochip_high_control_hep$BARCODE)){
  if(i/6. > count){
    count = count + 1
  }
  
  files1[count, (i - (count-1)*6)] <- paste( "S:/Madak-Erdogan Lab/Brandi/DAS_CEL_Files/single(early)_dose/HHyp/",
                                             df_24hr_nochip_high_control_hep$COMPOUND_NAME[i], ".Rat.in_vivo.Liver.Single/celfiles/",
                                             df_24hr_nochip_high_control_hep$BARCODE[i], ".", "CEL", sep="")}

###create labels for CEL/sample files
x <- list()
for (i in 1:length(df_24hr_nochip_high_control_hep$SACRI_PERIOD))
{ x <- paste(df_24hr_nochip_high_control_hep$COMPOUND.Abbr., ".", 
             df_24hr_nochip_high_control_hep$DOSE_LEVEL[], ".",
             seq(1:3), sep="")}
df_24hr_nochip_high_control_hep$label <- x


library(affy)
dimension <- dim(files1)
for(i in 1: dimension[1])
  assign(paste0("test", i), cbind(data.frame(exprs(rma(read.affybatch(files1[i,]))))))

dat<-cbind(test1,test2,test3,test4,test5)
names(dat) <-df_24hr_nochip_high_control_hep$label[1:length(files1)]
library(limma)
design <- model.matrix(~0 + df_24hr_nochip_high_control_hep$DOSE_LEVEL[1:length(files1)])
design <- design[,-c(3,4)]
colnames(design) <- c("Control","High")
rownames(design) <- df_24hr_nochip_high_control_hep$label[1:length(files1)]   
fit <- lmFit(dat, design)
cont.matrix <- makeContrasts(CvsH=Control-High, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEG <-topTable(fit2, n=Inf, adjust.method ="BH")
select <- DEG$adj.P.Val<0.01
length(which(select==TRUE))

#source("https://bioconductor.org/biocLite.R")
#biocLite("rat2302.db")
library(rat2302.db)
Annot <- data.frame(ACCNUM=sapply(contents(rat2302ACCNUM),
                                  paste, collapse=", "), SYMBOL=sapply(contents(rat2302SYMBOL), 
                                                                       paste, collapse=", "), DESC=sapply(contents(rat2302GENENAME), paste, collapse=", "))


DEG_annot <- merge(Annot, DEG, by.x=0, by.y=0, all=T)
rownames(DEG_annot)<-DEG_annot$Row.names
all.equal(rownames(DEG_annot), rownames(dat))
DEG.list<-cbind(DEG_annot, dat)
DEG.list.sig <- DEG.list[which(DEG.list$adj.P.Val < 0.01),]
row.names(DEG.list.sig) <- DEG.list.sig$SYMBOL
#install.packages("gplots")
library(gplots)
heatmap.2(as.matrix(t(DEG.list.sig[11:40]), scale="column", cexRow=1.5, labRow=df_24hr_nochip_high_control_hep$label), tracecol = NULL)

#generate a heatmap 
#Get GO Terms

write.csv(DEG.list, "HepHypDEG.csv")
save.image("HepDEG.RData")

#Centrilobular Hypertrophy
setwd("S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose\\CHyp")
## Getting the tsv file names
csv_files2 <- dir(path= 'S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose\\CHyp', pattern='*.tsv$', recursive = T)
## Merging the tsv files
for(i in 1:length(csv_files2)) {
  if(i == 1)
    df2 <- read.csv(csv_files2[i], sep = '\t')
  else
    df2 <- rbind(df2, read.csv(csv_files2[i], sep = '\t'))
}
## Getting subset of df (24 hr, No 'No Chip', High or Control)
df_24hr_nochip_high_control_cen = subset(df2, (df2$SACRI_PERIOD == '24 hr'& df2$BARCODE != 'No ChipData') & (df2$DOSE_LEVEL == 'High' | df2$DOSE_LEVEL == 'Control') )
## BARCODE Retrieval
files2 <- matrix(data=NA, nrow = 24, ncol= 6)
count <- 1
for( i in 1:length(df_24hr_nochip_high_control_cen$BARCODE)){
  if(i/6. > count){
    count = count + 1
  }
  
  files2[count, (i - (count-1)*6)] <- paste( "S:/Madak-Erdogan Lab/Brandi/DAS_CEL_Files/single(early)_dose/CHyp/",
                                             df_24hr_nochip_high_control_cen$COMPOUND_NAME[i], ".Rat.in_vivo.Liver.Single/celfiles/",
                                             df_24hr_nochip_high_control_cen$BARCODE[i], ".", "CEL", sep="")}
###create labels for CEL/sample files
x <- list()
for (i in 1:length(df_24hr_nochip_high_control_cen$SACRI_PERIOD))
{ x <- paste(df_24hr_nochip_high_control_cen$COMPOUND.Abbr., ".", 
             df_24hr_nochip_high_control_cen$DOSE_LEVEL[], ".",
             seq(1:3), sep="")}
df_24hr_nochip_high_control_cen$label <- x


#library(affy)
dimension <- dim(files2)
for(i in 1: dimension[1])
  assign(paste0("test", i), cbind(data.frame(exprs(rma(read.affybatch(files2[i,]))))))
#if true combine
dat1<-cbind( test1,test2,test3,test4,test5, test6,test7,test8,test9,test10,test11, test12,
             test13,test14,test15,test16,test17,test18,test19, test20,test21,test22,test23,test24)
names(dat1) <-df_24hr_nochip_high_control_cen$label[1:length(files2)]

library(limma)
design <- model.matrix(~0 + df_24hr_nochip_high_control_cen$DOSE_LEVEL[1:length(files2)])
design <- design[,-c(3,4)]
colnames(design) <- c("Control","High")
rownames(design) <- df_24hr_nochip_high_control_cen$label[1:length(files2)]   
fit <- lmFit(dat1, design)
cont.matrix <- makeContrasts(CvsH=Control-High, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEG <-topTable(fit2, n=Inf)

DEG_annot <- merge(Annot, DEG, by.x=0, by.y=0, all=T)
rownames(DEG_annot)<-DEG_annot$Row.names
all.equal(rownames(DEG_annot), rownames(dat1))
DEG.list2<-cbind(DEG_annot, dat1)
DEG.list2.sig <- DEG.list2[which(DEG.list2$adj.P.Val < 0.01),]
row.names(DEG.list2.sig) <- DEG.list2.sig$SYMBOL
#install.packages("gplots")
library(gplots)
heatmap.2(as.matrix(t(DEG.list2.sig[11:40]), scale="column", cexRow=1.5, labRow=df_24hr_nochip_high_control_cen$label), tracecol = NULL)

all.equal(rownames(dat), rownames(test1))
dat<-cbind(dat,dat1)
write.csv(DEG.list2, "CenHypDEG.csv")
save.image("CenHypDEG.RData")
all.equal(rownames(DEG.list), rownames(DEG.list2))
#if true combine
#DEG.list <-cbind(DEG.list, DEG.list2)



#NECROSIS
setwd("S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose\\Nec")
## Getting the tsv file names
csv_files3 <- dir(path= 'S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose\\Nec', pattern='*.tsv$', recursive = T)
## Merging the tsv files
for(i in 1:length(csv_files3)) {
  if(i == 1)
    df3 <- read.csv(csv_files3[i], sep = '\t')
  else
    df3 <- rbind(df3, read.csv(csv_files3[i], sep = '\t'))
}
## Getting subset of df (24 hr, No 'No Chip', High or Control)
df_24hr_nochip_high_control_nec = subset(df3, (df3$SACRI_PERIOD == '24 hr'& df3$BARCODE != 'No ChipData') & (df3$DOSE_LEVEL == 'High' | df3$DOSE_LEVEL == 'Control') )
## BARCODE Retrieval
files3 <- matrix(data=NA, nrow = 42, ncol= 6)
count <- 1
for( i in 1:length(df_24hr_nochip_high_control$BARCODE)){
  if(i/6. > count){
    count = count + 1
  }
  
  files3[count, (i - (count-1)*6)] <- paste( "S:/Madak-Erdogan Lab/Brandi/DAS_CEL_Files/single(early)_dose/Nec/",
                                             df_24hr_nochip_high_control_nec$COMPOUND_NAME[i], ".Rat.in_vivo.Liver.Single/celfiles/",
                                             df_24hr_nochip_high_control_nec$BARCODE[i], ".", "CEL", sep="")}

###create labels for CEL/sample files
x <- list()
for (i in 1:length(df_24hr_nochip_high_control_nec$SACRI_PERIOD))
{ x <- paste(df_24hr_nochip_high_control_nec$COMPOUND.Abbr., ".", 
             df_24hr_nochip_high_control_nec$DOSE_LEVEL[], ".",
             seq(1:3), sep="")}
df_24hr_nochip_high_control_nec$label <- x

#library(affy)
dimension <- dim(files3)
for(i in 1: dimension[1])
  assign(paste0("test", i), cbind(data.frame(exprs(rma(read.affybatch(files3[i,]))))))

dat2<-cbind(test1,test2,test3,test4,test5, test6,test7,test8,test9,test10,test11, test12,
            test13,test14,test15,test16,test17,test18,test19, test20,test21,test22,test23,test24,
            test25,test26,test27,test28,test29,test30,test31,test32,test33,test34,test35,test36,
            test37,test38,test39,test40,test41,test42)

names(dat2) <-df_24hr_nochip_high_control$label[1:length(files3)]

library(limma)
design <- model.matrix(~0 + df_24hr_nochip_high_control_nec$DOSE_LEVEL[1:length(files3)])
design <- design[,-c(3,4)]
colnames(design) <- c("Control","High")
rownames(design) <- df_24hr_nochip_high_control_nec$label[1:length(files3)]   
fit <- lmFit(dat2, design)
cont.matrix <- makeContrasts(CvsH=Control-High, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEG <-topTable(fit2, n=Inf)

DEG_annot <- merge(Annot, DEG, by.x=0, by.y=0, all=T)
rownames(DEG_annot)<-DEG_annot$Row.names
all.equal(rownames(DEG_annot), rownames(dat2))
DEG.list3<-cbind(DEG_annot, dat2)
DEG.list3<-cbind(DEG_annot, dat)
DEG.list3.sig <- DEG.list3[which(DEG.list3$adj.P.Val < 0.01),]
#SOME SYMBOLS UNKNOWN ; CANT USE; row.names(DEG.list2.sig) <- DEG.list2.sig$SYMBOL
#install.packages("gplots")
library(gplots)
heatmap.2(as.matrix(t(log2(DEG.list2.sig[11:40])), scale="row", cexRow=1.5, labRow=df_24hr_nochip_high_control$label), tracecol = NULL)

all.equal(rownames(dat2), rownames(test1))
dat<-cbind(dat,dat2)
write.csv(DEG.list2, "NecDEG.csv")
save.image("NecDEG.RData")
all.equal(rownames(DEG.list), rownames(DEG.list2))
#if true combine
DEG.list <-cbind(DEG.list, DEG.list2)

#CenNec
setwd("S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose\\CenNec")
## Getting the tsv file names
csv_files4 <- dir(path= 'S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose\\CenNec', pattern='*.tsv$', recursive = T)
## Merging the tsv files
for(i in 1:length(csv_files4)) {
  if(i == 1)
    df4 <- read.csv(csv_files4[i], sep = '\t')
  else
    df4 <- rbind(df4, read.csv(csv_files4[i], sep = '\t'))
}
## Getting subset of df4 (24 hr, No 'No Chip', High or Control)
df_24hr_nochip_high_control_cen_nec = subset(df4, (df4$SACRI_PERIOD == '24 hr'& df4$BARCODE != 'No ChipData') & (df4$DOSE_LEVEL == 'High' | df4$DOSE_LEVEL == 'Control') )
## BARCODE Retrieval
files4 <- matrix(data=NA, nrow = 17, ncol= 6)
count <- 1
for( i in 1:length(df_24hr_nochip_high_control_cen_nec$BARCODE)){
  if(i/6. > count){
    count = count + 1
  }
  
  files4[count, (i - (count-1)*6)] <- paste( "S:/Madak-Erdogan Lab/Brandi/DAS_CEL_Files/single(early)_dose/CenNec/",
                                             df_24hr_nochip_high_control_cen_nec$COMPOUND_NAME[i], ".Rat.in_vivo.Liver.Single/celfiles/",
                                             df_24hr_nochip_high_control_cen_nec$BARCODE[i], ".", "CEL", sep="")}

###create labels for CEL/sample files
x <- list()
for (i in 1:length(df_24hr_nochip_high_control_cen_nec$SACRI_PERIOD))
{ x <- paste(df_24hr_nochip_high_control_cen_nec$COMPOUND.Abbr., ".", 
             df_24hr_nochip_high_control_cen_nec$DOSE_LEVEL[], ".",
             seq(1:3), sep="")}
df_24hr_nochip_high_control_cen_nec$label <- x

#library(affy)
dimension <- dim(files4)
for(i in 1: dimension[1])
  assign(paste0("test", i), cbind(data.frame(exprs(rma(read.affybatch(files4[i,]))))))

dat3<-cbind(test1,test2,test3,test4,test5, test6,test7,test8,test9,test10,test11, test12,
            test13,test14,test15,test16,test17)

names(dat3) <-df_24hr_nochip_high_control_cen_nec$label[1:length(files4)]
all.equal(rownames(dat3), rownames(test1))
dat<-cbind(dat,dat3)
library(limma)
design <- model.matrix(~0 + df_24hr_nochip_high_control_cen_nec$DOSE_LEVEL[1:length(files4)])
design <- design[,-c(3,4)]
colnames(design) <- c("Control","High")
rownames(design) <- df_24hr_nochip_high_control_cen_nec$label[1:length(files4)]   
fit <- lmFit(dat, design)
cont.matrix <- makeContrasts(CvsH=Control-High, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEG <-topTable(fit2, n=Inf)

DEG_annot <- merge(Annot, DEG, by.x=0, by.y=0, all=T)
rownames(DEG_annot)<-DEG_annot$Row.names
all.equal(rownames(DEG_annot), rownames(dat))
DEG.list2<-cbind(DEG_annot, dat)
write.csv(DEG.list2, "CenNecDEG.csv")
save.image("CenNecDEG.RData")
all.equal(rownames(DEG.list), rownames(DEG.list2))
#if true combine
DEG.list <-cbind(DEG.list, DEG.list2)



#HepNec
setwd("S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose\\HepNec")
## Getting the tsv file names
csv_files5 <- dir(path= 'S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose\\HepNec', pattern='*.tsv$', recursive = T)
## Merging the tsv files
for(i in 1:length(csv_files5)) {
  if(i == 1)
    df5 <- read.csv(csv_files5[i], sep = '\t')
  else
    df5 <- rbind(df5, read.csv(csv_files5[i], sep = '\t'))
}
## Getting subset of df (24 hr, No 'No Chip', High or Control)
df_24hr_nochip_high_control_hec_nec = subset(df5, (df5$SACRI_PERIOD == '24 hr'& df5$BARCODE != 'No ChipData') & (df5$DOSE_LEVEL == 'High' | df5$DOSE_LEVEL == 'Control') )
## BARCODE Retrieval
files <- matrix(data=NA, nrow = 6, ncol= 6)
count <- 1
for( i in 1:length(df_24hr_nochip_high_control_hec_nec$BARCODE)){
  if(i/6. > count){
    count = count + 1
  }
  
  files[count, (i - (count-1)*6)] <- paste( "S:/Madak-Erdogan Lab/Brandi/DAS_CEL_Files/single(early)_dose/HepNec/",
                                            df_24hr_nochip_high_control_hec_nec$COMPOUND_NAME[i], ".Rat.in_vivo.Liver.Single/celfiles/",
                                            df_24hr_nochip_high_control_hec_nec$BARCODE[i], ".", "CEL", sep="")}

###create labels for CEL/sample files
x <- list()
for (i in 1:length(df_24hr_nochip_high_control_hec_nec$SACRI_PERIOD))
{ x <- paste(df_24hr_nochip_high_control_hec_nec$COMPOUND.Abbr., ".", 
             df_24hr_nochip_high_control_hec_nec$DOSE_LEVEL[], ".",
             seq(1:3), sep="")}
df_24hr_nochip_high_control_hec_nec$label <- x

#library(affy)
dimension <- dim(files)
for(i in 1: dimension[1])
  assign(paste0("test", i), cbind(data.frame(exprs(rma(read.affybatch(files[i,]))))))

dat4<-cbind(test1,test2,test3,test4,test5, test6)

names(dat4) <-df_24hr_nochip_high_control_hec_nec$label[1:length(files)]
all.equal(rownames(dat), rownames(test1))
dat<-cbind(dat,dat4)

library(limma)
design <- model.matrix(~0 + df_24hr_nochip_high_control_hec_nec$DOSE_LEVEL[1:length(files)])
design <- design[,-c(3,4)]
colnames(design) <- c("Control","High")
rownames(design) <- df_24hr_nochip_high_control_hec_nec$label[1:length(files)]   
fit <- lmFit(dat, design)
cont.matrix <- makeContrasts(CvsH=Control-High, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEG <-topTable(fit2, n=Inf)

DEG_annot <- merge(Annot, DEG, by.x=0, by.y=0, all=T)
rownames(DEG_annot)<-DEG_annot$Row.names
all.equal(rownames(DEG_annot), rownames(dat))
DEG.list2<-cbind(DEG_annot, dat)
write.csv(DEG.list2, "HepNecDEG.csv")
save.image("HepNecDEG.RData")
setwd("S:\\Madak-Erdogan Lab\\Brandi\\DAS_CEL_Files\\single(early)_dose")
write.csv(dat, "allexpr.csv")
save.image("allDEG.RData")
