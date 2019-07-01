library(tidyverse)
library(ggpubr)
library(cowplot)
library(Hmisc)
library(scales)
library(ggrepel)


###read cellprofiler output files and convert metadata to factor ####

#lysozyme quant in pgc tumors

pgclystumors <- read.csv("~/Dropbox/JEV/CellProfiler/pgckolystumorsv1Cytoplasm.csv")
pgclystumors$Metadata_Genotype <- as.factor(pgclystumors$Metadata_Genotype)
pgclystumors$Metadata_Mouse <- as.factor(pgclystumors$Metadata_Mouse)
pgclystumors$Metadata_Tumor <- as.factor(pgclystumors$Metadata_Tumor)
pgclystumors$Metadata_Genotype <- ordered(pgclystumors$Metadata_Genotype, c("PGCWT", "PGCKO"))


###statistics - treat tumors as independent replicates. aggregate individual cell data by mouse, tumor, and week. ####

z <- pgclystumors


#aggregate most common quants for stats
aggregate.z <- z %>% group_by(Metadata_Mouse, Metadata_Tumor, Metadata_Genotype) %>% summarise(mean_intensity = mean(Intensity_MeanIntensity_OrigRed), median_intensity = mean(Intensity_MedianIntensity_OrigRed), area = mean(AreaShape_Area))
wilcox.test(median_intensity ~ Metadata_Genotype, data = aggregate.z)

#aggregate nuclear stains quants for stats
aggregate.nuc <- z %>% group_by(Metadata_Mouse, Metadata_Tumor, Metadata_Genotype) %>% summarise(mean_intensity = mean(Intensity_MeanIntensity_OrigGreen1), median_intensity = mean(Intensity_MedianIntensity_OrigGreen1), area = mean(AreaShape_Area))
wilcox.test(median_intensity ~ Metadata_Genotype, data = aggregate.nuc)


#violin graphs

color.cat <- "#00BFC4"
color.cathr <- "#F8766D"

ggplot(z, aes(x = Metadata_Genotype, y = Intensity_MedianIntensity_OrigRed, fill = Metadata_Genotype)) + 
  geom_violin(trim = T, bw = .015) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = "black", fun.args = list(mult = 1), show.legend = F) +
  scale_fill_manual(values = c(color.cat, color.cathr)) +
  labs(fill = NULL, y = "Median Fluorescent Intensity", title = "pgc tum lys p = .02881") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())



