library(tidyverse)
library(fgsea)
library(Matrix)
library(png)
library(cowplot)
library(calibrate)

#load pathways. sclung100 is top 100 markers of AT1, AT2, Club, and Ciliated cells, taken from Han et al. 2018. Hallmarks is hallmarks from msigdb plus sclung100.
hallmarks <- gmtPathways("~/Dropbox/JEV/GSEA/h.all.v6.1.mcalung.symbols.gmt.txt")
sclung100 <- gmtPathways("~/Dropbox/JEV/GSEA/mca.lungsymbols.v1.top100.gmt.txt")


wk2pi3k <- read.table("~/Dropbox/JEV/GSEA/Galaxy155-[DESeq2_2_week_pi3k].tabular", header = FALSE, colClasses = c("character", "NULL", "numeric", "NULL", "NULL", "NULL", "NULL"), col.names = c("ID", ".", "t", ".", ".", ".", "."))
wk2pi3k <- wk2pi3k[complete.cases(wk2pi3k),]
wk2pi3k[,1] = toupper(wk2pi3k[,1])
wk2pi3k <- setNames(wk2pi3k$t, wk2pi3k$ID)
wk2pi3k.hallmarks <- fgsea(hallmarks, wk2pi3k, minSize = 15, maxSize = 200, nperm = 1000)

wk6pi3k <- read.table("~/Dropbox/JEV/GSEA/Galaxy158-[DESeq2_6_week_pi3k].tabular", header = FALSE, colClasses = c("character", "NULL", "numeric", "NULL", "NULL", "NULL", "NULL"), col.names = c("ID", ".", "t", ".", ".", ".", "."))
wk6pi3k <- wk6pi3k[complete.cases(wk6pi3k),]
wk6pi3k[,1] = toupper(wk6pi3k[,1])
wk6pi3k <- setNames(wk6pi3k$t, wk6pi3k$ID)
wk6pi3k.hallmarks <- fgsea(hallmarks, wk6pi3k, minSize = 15, maxSize = 500, nperm = 10000)

wk12pi3k <- read.table("~/Dropbox/JEV/GSEA/Galaxy161-[DESeq2_12_week_pi3k].tabular", header = FALSE, colClasses = c("character", "NULL", "numeric", "NULL", "NULL", "NULL", "NULL"), col.names = c("ID", ".", "t", ".", ".", ".", "."))
wk12pi3k <- wk12pi3k[complete.cases(wk12pi3k),]
wk12pi3k[,1] = toupper(wk12pi3k[,1])
wk12pi3k <- setNames(wk12pi3k$t, wk12pi3k$ID)
wk12pi3k.hallmarks <- fgsea(hallmarks, wk12pi3k, minSize = 15, maxSize = 500, nperm = 10000)

plotEnrichment(sclung100[["AT2"]], wk2pi3k)
plotEnrichment(sclung100[["AT2"]], wk6pi3k)
plotEnrichment(sclung100[["AT2"]], wk12pi3k)
