#retrieve and combine attributes in single dataframe
setwd("/Users/brandismith/Box/Dow/Nec")#pathToSingleDoseFolders
files <- list.files(path= ".", pattern='*.tsv$', recursive = T)

## Retrieving chip BARCODES for 24 hr dose times
df_24hr_nochip_high_control= subset(df, (df$SACRI_PERIOD == '24 hr'& df$BARCODE != 'No ChipData') & (df$DOSE_LEVEL == 'High' | df$DOSE_LEVEL == 'Control') )
bar_codes <- matrix(data=NA, nrow = 42, ncol= 6)
count <- 1
for( i in 1:length(df_24hr_nochip_high_control$BARCODE)){
  if(i/6. > count){
    count = count + 1
    }
  bar_codes[count, (i - (count-1)*6)] <- paste( "/Users/brandismith/Box/Dow/Nec/", 
                                                df_24hr_nochip_high_control$COMPOUND_NAME[i], ".Rat.in_vivo.Liver.Single/celfiles/",
                                                df_24hr_nochip_high_control$BARCODE[i], ".", "CEL", sep="")
  }

###create labels column in attribute file to differentiate control and high dose barcodes
x <- list()
for (i in 1:length(df_24hr_nochip_high_control$SACRI_PERIOD))
{ x <- paste(df_24hr_nochip_high_control$COMPOUND.Abbr., ".", 
             df_24hr_nochip_high_control$DOSE_LEVEL[], ".",
             seq(1:3), sep="")}
df_24hr_nochip_high_control$label <- x

#install affy, limma, and annotations and retrieve gene expressions from CEL files
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("rat2302.db")
library(rat2302.db)
library(affy)
dimension <- dim(bar_codes)
#this may take some time...
for(i in 1: dimension[1])
  assign(paste0("test", i), cbind(data.frame(exprs(rma(read.affybatch(bar_codes[i,]))))))

df_exprs<-cbind(test1,test2,test3,test4,test5, test6,test7,test8,test9,test10,test11, test12,
            test13,test14,test15,test16,test17,test18,test19, test20,test21,test22,test23,test24,
            test25,test26,test27,test28,test29,test30,test31,test32,test33,test34,test35,test36,
            test37,test38,test39,test40,test41,test42)

names(df_exprs) <-df_24hr_nochip_high_control$label[1:length(bar_codes)]

BiocManager::install("limma")
library(limma)
design <- model.matrix(~0 + df_24hr_nochip_high_control$DOSE_LEVEL[1:length(bar_codes)])
design <- design[,-c(3,4)]
colnames(design) <- c("Control","High")
rownames(design) <- df_24hr_nochip_high_control$label[1:length(bar_codes)]   
fit <- lmFit(df_exprs, design)
cont.matrix <- makeContrasts(CvsH=Control-High, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
DEG <-topTable(fit2, n=Inf)

#matching chip names with annotated numbers from rat database to get list of genes
Annot <- data.frame(ACCNUM=sapply(contents(rat2302ACCNUM),
                                  paste, collapse=", "), SYMBOL=sapply(contents(rat2302SYMBOL), 
                                                                       paste, collapse=", "), DESC=sapply(contents(rat2302GENENAME), paste, collapse=", "))


DEG_annot <- merge(Annot, DEG, by.x=0, by.y=0, all=T)
rownames(DEG_annot)<-DEG_annot$Row.names
all.equal(rownames(DEG_annot), rownames(df_exprs))
DEG.list.sig <- DEG_annot[which(DEG_annot$adj.P.Val < 0.01),]
write.csv(DEG.list.sig, "NecDEG.csv")