library(tidyverse)
library(fgsea)
library(Matrix)
library(png)
library(cowplot)
library(calibrate)

#load pathways
hallmarks <- gmtPathways("~/Dropbox/JEV/GSEA/h.all.v6.1.mcalung.symbols.gmt.txt")


pgckohet_101220 <- read.table("~/Dropbox/JEV/GSEA/Galaxy753-[DESeq2_result_file_on_data_694,_data_707,_and_others].tabular", header = FALSE, colClasses = c("character", "NULL", "numeric", "NULL", "NULL", "NULL", "NULL"), col.names = c("ID", ".", "t", ".", ".", ".", "."))
pgckohet_101220 <- pgckohet_101220[complete.cases(pgckohet_101220),]
pgckohet_101220[,1] = toupper(pgckohet_101220[,1])
pgckohet_101220 <- setNames(pgckohet_101220$t, pgckohet_101220$ID)

pgckohet_101220.hallmarks <- fgsea(hallmarks, pgckohet_101220, minSize = 15, maxSize = 500, nperm = 1000)


t <- pgckohet_101220.hallmarks[,1:7]

#bar graph in main figure 7
ggplot(data = t, aes(x = reorder(t$pathway, t$NES), y = t$NES, fill = cut(t$padj, c(0, .05, 1)))) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(labels = c("<.05", "≥.05"), values = c("black", "grey")) +
  coord_flip() + 
  labs(x = "Pathway", y = "NES", fill = "Adj. P Value") + 
  theme(axis.text.y = element_text(size = 10))

#mountain plot in figure supplement 1
plotEnrichment(sclung100[["AT2"]], pgckohet_101220) 
