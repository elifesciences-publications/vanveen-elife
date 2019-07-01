library(tidyverse)
library(fgsea)
library(Matrix)
library(png)
library(cowplot)
library(calibrate)

#load pathways. sclung100 is top 100 markers of AT1, AT2, Club, and Ciliated cells, taken from Han et al. 2018. Hallmarks is hallmarks from msigdb plus sclung100.
hallmarks <- gmtPathways("~/Dropbox/JEV/GSEA/h.all.v6.1.mcalung.symbols.gmt.txt")
sclung100 <- gmtPathways("~/Dropbox/JEV/GSEA/mca.lungsymbols.v1.top100.gmt.txt")


allpi3krin7 <- read.table("~/Dropbox/JEV/GSEA/Galaxy173-[DESeq2_pi3k_allwk_rin7].tabular", header = FALSE, colClasses = c("character", "NULL", "numeric", "NULL", "NULL", "NULL", "NULL"), col.names = c("ID", ".", "t", ".", ".", ".", "."))
allpi3krin7 <- allpi3krin7[complete.cases(allpi3krin7),]
allpi3krin7[,1] = toupper(allpi3krin7[,1])
allpi3krin7 <- setNames(allpi3krin7$t, allpi3krin7$ID)

sc100allpi3krin7 <- fgsea(sclung100, allpi3krin7, minSize = 15, maxSize = 500, nperm = 1000)
allpi3krin7.hallmarks <- fgsea(hallmarks, allpi3krin7, minSize = 15, maxSize = 500, nperm = 1000)

#bar graph showing all hallmark + identity markers

t <- allpi3krin7.hallmarks[, 1:7]
ggplot(data = t, aes(x = reorder(t$pathway, t$NES), y = t$NES, fill = cut(t$padj, c(0, .05, 1)))) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(labels = c("<.05", "≥.05"), values = c("black", "grey")) +
  coord_flip() + 
  labs(x = "Pathway", y = "NES", fill = "Adj. P Value") + 
  theme(axis.text.y = element_text(size = 10))

#individual mountain plots

plotEnrichment(hallmarks[["HALLMARK_PI3K_AKT_MTOR_SIGNALING"]], allpi3krin7)
plotEnrichment(hallmarks[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]], allpi3krin7)
plotEnrichment(sclung100[["AT2"]], allpi3krin7)
plotEnrichment(sclung100[["AT1"]], allpi3krin7)

#separate time points

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



t1 <- wk2pi3k.hallmarks[, 1:7]
t2 <- wk6pi3k.hallmarks[,1:7]
t0 <- wk12pi3k.hallmarks[, 1:7]
t3 <- merge(t1, t2, by.x = "pathway", by.y = "pathway", all = T)
t3 <- merge(t3, t0, by.x = "pathway", by.y = "pathway", all = T)
t3 <- column_to_rownames(t3, var = "pathway")
#t3 <- subset(t3, padj.x <= .05 | padj.y <= .05 | padj <= .05)



#this one captures the most information - exploded view with time points. Shown in supplement 1 to figure 2
ggplot(data = t3) + 
  geom_point(aes(x = reorder(row.names(t3), -(padj + padj.x + padj.y)), y = NES.x, size = -log(padj.x), color = "a  2 Week"), alpha = .8) +
  geom_point(aes(x = reorder(row.names(t3), -(padj + padj.x + padj.y)), y = NES.y, size = -log(padj.y), color = "b  6 Week"), alpha = .8) +
  geom_point(aes(x = reorder(row.names(t3), -(padj + padj.x + padj.y)), y = NES, size = -log(padj), color = "c  12 Week"), alpha = .8) +
  scale_size_area(name = "P Adj.", breaks = c(1, 1.3, 3, 4), labels = c(".1", ".05", ".001", ".0001")) +
  theme(axis.text.y = element_text(size = 10)) +
  labs(x = "Pathway", y = "NES", color = "Weeks P.I.") +
  coord_flip()



