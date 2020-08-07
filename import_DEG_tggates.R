#retrieving CEL file attributes 
setwd("/Users/brandismith/Box/Dow/Nec")#pathToSingleDoseFolders
files <- list.files(path= ".", pattern='*.tsv$', recursive = T)
for(i in 1:length(files)) {
  if(i == 1)
    df <- read.csv(files[i], sep = '\t')
  else
    df <- rbind(df, read.csv(files[i], sep = '\t'))
}
## Retrieving chip BARCODES for 24 hr dose times
df_24hr_nochip_high_control= subset(df, (df$SACRI_PERIOD == '24 hr'& df$BARCODE != 'No ChipData') & (df$DOSE_LEVEL == 'High' | df$DOSE_LEVEL == 'Control') )
bar_codes <- matrix(data=NA, nrow = 42, ncol= 6)
count <- 1
for( i in 1:length(df_24hr_nochip_high_control$BARCODE)){
  if(i/6. > count){
    count = count + 1
  }
  
  bar_codes[count, (i - (count-1)*6)] <- paste( "/Users/brandismith/Box/Dow/Nec",
                                             df_24hr_nochip_high_control$COMPOUND_NAME[i], ".Rat.in_vivo.Liver.Single/celfiles/",
                                             df_24hr_nochip_high_control$BARCODE[i], ".", "CEL", sep="")}

###create labels for CEL/sample files
x <- list()
for (i in 1:length(df_24hr_nochip_high_control$SACRI_PERIOD))
{ x <- paste(df_24hr_nochip_high_control$COMPOUND.Abbr., ".", 
             df_24hr_nochip_high_control$DOSE_LEVEL[], ".",
             seq(1:3), sep="")}
df_24hr_nochip_high_control$label <- x

#install affy and retrieve gene expressions from CEL files
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")
library(affy)
dimension <- dim(bar_codes)
for(i in 1: dimension[1])
  assign(paste0("test", i), cbind(data.frame(exprs(rma(read.affybatch(bar_codes[i,]))))))

dat2<-cbind(test1,test2,test3,test4,test5, test6,test7,test8,test9,test10,test11, test12,
            test13,test14,test15,test16,test17,test18,test19, test20,test21,test22,test23,test24,
            test25,test26,test27,test28,test29,test30,test31,test32,test33,test34,test35,test36,
            test37,test38,test39,test40,test41,test42)

names(dat2) <-df_24hr_nochip_high_control$label[1:length(bar_codes)]

library(limma)
design <- model.matrix(~0 + df_24hr_nochip_high_control$DOSE_LEVEL[1:length(bar_codes)])
design <- design[,-c(3,4)]
colnames(design) <- c("Control","High")
rownames(design) <- df_24hr_nochip_high_control$label[1:length(bar_codes)]   
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
