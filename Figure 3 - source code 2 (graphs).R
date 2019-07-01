library(tidyverse)
library(ggpubr)
library(cowplot)
library(Hmisc)
library(scales)
library(ggrepel)


###read cellprofiler output files and convert metadata to factor ####

#spa quantitation, primary object DAPI

spatumors <- read.csv("~/Dropbox/JEV/CellProfiler/spatumorsv1Cytoplasm.csv")
spatumors$Metadata_Genotype <- as.factor(spatumors$Metadata_Genotype)
spatumors$Metadata_Mouse <- as.factor(spatumors$Metadata_Mouse)

#spc quantitation, primary object DAPI

spctumors <- read.csv("~/Dropbox/JEV/CellProfiler/spctumorsv1Cytoplasm.csv")
spctumors$Metadata_Genotype <- as.factor(spctumors$Metadata_Genotype)
spctumors$Metadata_Mouse <- as.factor(spctumors$Metadata_Mouse)

#lys quant - DAPI

lystumors <- read.csv("~/Dropbox/JEV/CellProfiler/lystumorsv3Cytoplasm.csv")
lystumors$Metadata_Genotype <- as.factor(lystumors$Metadata_Genotype)
lystumors$Metadata_Mouse <- as.factor(lystumors$Metadata_Mouse)


###statistics - treat tumors as independent replicates. aggregate individual cell data by mouse, tumor, and week. ####
z <- spatumors
z <- spctumors
z <- lystumors


aggregate.z <- z %>% group_by(Metadata_Mouse, Metadata_Tumor, Metadata_Genotype) %>% summarise(mean_intensity = mean(Intensity_MeanIntensity_OrigRed), median_intensity = mean(Intensity_MedianIntensity_OrigRed), area = mean(AreaShape_Area))
wilcox.test(median_intensity ~ Metadata_Genotype, data = aggregate.z)

color.cat <- "#00BFC4"
color.cathr <- "#F8766D"

ggplot(z, aes(x = Metadata_Genotype, y = Intensity_MedianIntensity_OrigRed, fill = Metadata_Genotype)) + 
  geom_violin(trim = T, bw = .015) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = "black", fun.args = list(mult = 1), show.legend = F) +
  scale_fill_manual(values = c(color.cat, color.cathr)) +
  labs(fill = NULL, y = "Median Fluorescent Intensity", title = "pgc tum lys p = .02881") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())