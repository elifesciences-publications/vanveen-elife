library(tidyverse)
library(ggpubr)
library(cowplot)
library(Hmisc)
library(scales)
library(ggrepel)


###read cellprofiler output files and convert metadata to factor ####
#spa quantitation with primary object being NKX in cat vs cat;hr. Also measure nkx in nuclei and cyto and sizes

cathrspatumors <- read.csv("~/Dropbox/JEV/CellProfiler/catnkxtumorsv6Cytoplasm.csv")
cathrspatumors$Metadata_Genotype <- as.factor(cathrspatumors$Metadata_Genotype)
cathrspatumors$Metadata_Mouse <- as.factor(cathrspatumors$Metadata_Mouse)
cathrspatumors.nuc <- read.csv("~/Dropbox/JEV/CellProfiler/catnkxtumorsv6Nuclei.csv")
cathrspatumors.nuc$Metadata_Genotype <- as.factor(cathrspatumors.nuc$Metadata_Genotype)
cathrspatumors.nuc$Metadata_Mouse <- as.factor(cathrspatumors.nuc$Metadata_Mouse)

cathrspatumors.nuc.cyto <- bind_cols(cathrspatumors, cathrspatumors.nuc)

#2wk spa quantitation

SPA2wk <- read.csv("~/Dropbox/JEV/CellProfiler/2wktumorsv5Cytoplasm.csv")
SPA2wk$Metadata_Genotype <- as.factor(SPA2wk$Metadata_Genotype)
SPA2wk$Metadata_Mouse <- as.factor(SPA2wk$Metadata_Mouse)
SPA2wk$Metadata_Tumor <- as.factor(SPA2wk$Metadata_Tumor)

SPA2wk.nuc <- read.csv("~/Dropbox/JEV/CellProfiler/2wktumorsv5Nuclei.csv")
SPA2wk.nuc$Metadata_Genotype <- as.factor(SPA2wk.nuc$Metadata_Genotype)
SPA2wk.nuc$Metadata_Mouse <- as.factor(SPA2wk.nuc$Metadata_Mouse)
SPA2wk.nuc$Metadata_Tumor <- as.factor(SPA2wk.nuc$Metadata_Tumor)

SPA2wk.nuc.cyto <- bind_cols(SPA2wk, SPA2wk.nuc)


#spa quantitation in lsl vs lsl;hr - nkx primary object. 


lslspatumors.cyto <- read.csv("~/Dropbox/JEV/CellProfiler/lslspankxtumorsv3Cytoplasm.csv")
lslspatumors.cyto$Metadata_Genotype <- as.factor(lslspatumors.cyto$Metadata_Genotype)
lslspatumors.cyto$Metadata_Mouse <- as.factor(lslspatumors.cyto$Metadata_Mouse)
lslspatumors.nuc <- read.csv("~/Dropbox/JEV/CellProfiler/lslspankxtumorsv3Nuclei.csv")
lslspatumors.nuc$Metadata_Genotype <- as.factor(lslspatumors.nuc$Metadata_Genotype)
lslspatumors.nuc$Metadata_Mouse <- as.factor(lslspatumors.nuc$Metadata_Mouse)

lslspatumors.nuc.cyto <- bind_cols(lslspatumors.cyto, lslspatumors.nuc)



###statistics - treat tumors as independent replicates. aggregate individual cell data by mouse, tumor, and week. ####

z <- SPA2wk
z <- SPA2wk.nuc
z <- cathrspatumors
z <- cathrspatumors.nuc
z <- lslspatumors.nuc.cyto

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



### two factor analyses ####

color.cat <- "blue"
color.cathr <- "red"

#test nkx vs spa in combined analyses

z <- SPA2wk.nuc.cyto
z <- cathrspatumors.nuc.cyto
z <- lslspatumors.nuc.cyto


#2wk - nuclear NKX is Green. 12wk - nuclear NKX is Green1, cytoplasmic NKX is Green. pgc similar to 12 wk - nuclear pgc is grayRed1. lsl similar too, nuclear nkx is green1



p3 <- ggplot(data = filter(z, Metadata_Genotype == "CAT"), aes(x = Intensity_MedianIntensity_OrigRed, y = Intensity_MedianIntensity_OrigGreen)) +
  stat_binhex(bins = 32) +
  xlim(0, .3) +
  ylim(0, .3) +
  labs(x = "Median SFTPA Intensity", y = "Median NKX Intensity", title = "SPA vs NKX staining in CAT tumors 2wk", subtitle = "Spearman's rank correlation Rho = .233")

p4 <- ggplot(data = filter(z, Metadata_Genotype == "CATHR"), aes(x = Intensity_MedianIntensity_OrigRed, y = Intensity_MedianIntensity_OrigGreen)) +
  stat_binhex(bins = 32) +
  xlim(0, .3) +
  ylim(0, .3) +
  labs(x = "Median SFTPA Intensity", y = "Median NKX Intensity", title = "SPA vs NKX staining in CATHR tumors 2wk", subtitle = "Spearman's rank correlation Rho = .078")


p5 <- ggplot(data = filter(z, Metadata_Genotype == "CAT"), aes(x = Intensity_MedianIntensity_OrigRed, y = Intensity_MedianIntensity_OrigGreen1)) +
  stat_binhex(bins = 40) +
  xlim(0, .3) +
  ylim(0, .3) +
  labs(x = "Median SFTPA Intensity", y = "Median NKX Intensity", title = "SPA vs NKX staining in CAT tumors 12wk", subtitle = "Spearman's rank correlation Rho = .279")

p6 <- ggplot(data = filter(z, Metadata_Genotype == "CATHR"), aes(x = Intensity_MedianIntensity_OrigRed, y = Intensity_MedianIntensity_OrigGreen1)) +
  stat_binhex(bins = 40) +
  xlim(0, .3) +
  ylim(0, .3) +
  labs(x = "Median SFTPA Intensity", y = "Median NKX Intensity", title = "SPA vs NKX staining in CATHR tumors 12wk", subtitle = "Spearman's rank correlation Rho = .407")


p7 <- ggplot(data = filter(z, Metadata_Genotype == "LSL"), aes(x = Intensity_MedianIntensity_OrigRed, y = Intensity_MedianIntensity_OrigGreen1)) +
  stat_binhex(bins = 40) +
  xlim(0, .5) +
  labs(x = "Median SFTPA Intensity", y = "Median NKX Intensity", title = "SPA vs NKX staining in LSL tumors 14wk", subtitle = "Spearman's rank correlation Rho = .256")

p8 <- ggplot(data = filter(z, Metadata_Genotype == "LSLHR"), aes(x = Intensity_MedianIntensity_OrigRed, y = Intensity_MedianIntensity_OrigGreen1)) +
  stat_binhex(bins = 40) +
  xlim(0,.5) +
  labs(x = "Median SFTPA Intensity", y = "Median NKX Intensity", title = "SPA vs NKX staining in LSLHR tumors 14wk", subtitle = "Spearman's rank correlation Rho = .203")


cor.test(x = filter(z, Metadata_Genotype == "CATHR")$Intensity_MedianIntensity_OrigRed, y = filter(z, Metadata_Genotype == "CATHR")$Intensity_MedianIntensity_OrigGreen1, method = "spearman")
cor.test(x = filter(z, Metadata_Genotype == "CAT")$Intensity_MedianIntensity_OrigRed, y = filter(z, Metadata_Genotype == "CAT")$Intensity_MedianIntensity_OrigGreen1, method = "spearman")

cor.test(x = filter(z, Metadata_Genotype == "LSLHR")$Intensity_MedianIntensity_OrigRed, y = filter(z, Metadata_Genotype == "LSLHR")$Intensity_MedianIntensity_OrigGreen1, method = "spearman")
cor.test(x = filter(z, Metadata_Genotype == "LSL")$Intensity_MedianIntensity_OrigRed, y = filter(z, Metadata_Genotype == "LSL")$Intensity_MedianIntensity_OrigGreen1, method = "spearman")

ggsave("~/Dropbox/JEV/CAT mouse paper/figures for resubmission/spankx2wkcat.tiff", plot = p3, width = 4, height = 4, dpi = 320, units = "in")
ggsave("~/Dropbox/JEV/CAT mouse paper/figures for resubmission/spankx2wkcathr.tiff", plot = p4, width = 4, height = 4, dpi = 320, units = "in")

ggsave("~/Dropbox/JEV/CAT mouse paper/figures for resubmission/spankx12wkcat.tiff", plot = p5, width = 4, height = 4, dpi = 320, units = "in")
ggsave("~/Dropbox/JEV/CAT mouse paper/figures for resubmission/spankx12wkcathr.tiff", plot = p6, width = 4, height = 4, dpi = 320, units = "in")

ggsave("~/Dropbox/JEV/CAT mouse paper/figures for resubmission/spankx14wklsl.tiff", plot = p7, width = 4, height = 4, dpi = 320, units = "in")
ggsave("~/Dropbox/JEV/CAT mouse paper/figures for resubmission/spankx14wklslhr.tiff", plot = p8, width = 4, height = 4, dpi = 320, units = "in")



#define quadrants by mean of CAT tumors minus one SD

z.cat <- filter(z, Metadata_Genotype == "CAT")
z.cathr <- filter(z, Metadata_Genotype == "CATHR")
cat.nkx.line <- (mean(z.cat$Intensity_MedianIntensity_OrigGreen) - sd(z.cat$Intensity_MedianIntensity_OrigGreen))
cat.sftpa.line <- (mean(z.cat$Intensity_MedianIntensity_OrigRed) - sd(z.cat$Intensity_MedianIntensity_OrigRed))


#quantitate fraction of cells in each quadrant
q1 <- count(filter(z, Metadata_Genotype == "CAT"), Intensity_MedianIntensity_OrigRed < cat.sftpa.line & Intensity_MedianIntensity_OrigGreen >= cat.nkx.line)
q2 <- count(filter(z, Metadata_Genotype == "CAT"), Intensity_MedianIntensity_OrigRed >= cat.sftpa.line & Intensity_MedianIntensity_OrigGreen >= cat.nkx.line)
q3 <- count(filter(z, Metadata_Genotype == "CAT"), Intensity_MedianIntensity_OrigRed >= cat.sftpa.line & Intensity_MedianIntensity_OrigGreen < cat.nkx.line)
q4 <- count(filter(z, Metadata_Genotype == "CAT"), Intensity_MedianIntensity_OrigRed < cat.sftpa.line & Intensity_MedianIntensity_OrigGreen < cat.nkx.line)
q5 <- count(filter(z, Metadata_Genotype == "CATHR"), Intensity_MedianIntensity_OrigRed < cat.sftpa.line & Intensity_MedianIntensity_OrigGreen >= cat.nkx.line)
q6 <- count(filter(z, Metadata_Genotype == "CATHR"), Intensity_MedianIntensity_OrigRed >= cat.sftpa.line & Intensity_MedianIntensity_OrigGreen >= cat.nkx.line)
q7 <- count(filter(z, Metadata_Genotype == "CATHR"), Intensity_MedianIntensity_OrigRed >= cat.sftpa.line & Intensity_MedianIntensity_OrigGreen < cat.nkx.line)
q8 <- count(filter(z, Metadata_Genotype == "CATHR"), Intensity_MedianIntensity_OrigRed < cat.sftpa.line & Intensity_MedianIntensity_OrigGreen < cat.nkx.line)



Xsq <- data.frame("CAT" = numeric(), "CATHR" = numeric())
Xsq[1,"CAT"] <- q1[2,2]
Xsq[2,"CAT"] <- q2[2,2]
Xsq[3,"CAT"] <- q3[2,2]
Xsq[4,"CAT"] <- q4[2,2]
Xsq[1,"CATHR"] <- q5[2,2]
Xsq[2,"CATHR"] <- q6[2,2]
Xsq[3,"CATHR"] <- q7[2,2]
Xsq[4,"CATHR"] <- q8[2,2]

chisq <- chisq.test(Xsq)
round(chisq$residuals, 3)

ggplot(data = z, aes(x = Intensity_MedianIntensity_OrigRed, y = Intensity_MedianIntensity_OrigGreen)) +
  stat_binhex(bins = 40, data = filter(z, Metadata_Genotype == "CAT"), fill = color.cat, aes(alpha = ..count..)) +
  stat_binhex(bins = 40, data = filter(z, Metadata_Genotype == "CATHR"), fill = color.cathr, aes(alpha = ..count..)) +
  labs(x = "Median SFTPA Intensity", y = "Median NKX Intensity", title = "X-squared = 421.24, df = 3, p-value < 2.2e-16") +
  xlim(0.05,.25) + 
  ylim(0.05,.25) +
  geom_hline(yintercept = cat.nkx.line, linetype="dashed", color = "black") +
  geom_vline(xintercept = cat.sftpa.line, linetype="dashed", color = "black")  +
  geom_text(size = 4, data = data.frame(x = .055, y = .25), aes(x = x, y = y, label = percent(as.numeric(q1[2,2]/(q1[1,2]+q1[2,2])))), inherit.aes = F, color = color.cat) +
  geom_text(size = 4, data = data.frame(x = .2, y = .25), aes(x = x, y = y, label = percent(as.numeric(q2[2,2]/(q1[1,2]+q1[2,2])))), inherit.aes = F, color = color.cat) +
  geom_text(size = 4, data = data.frame(x = .2, y = .055), aes(x = x, y = y, label = percent(as.numeric(q3[2,2]/(q1[1,2]+q1[2,2])))), inherit.aes = F, color = color.cat) +
  geom_text(size = 4, data = data.frame(x = .055, y = .055), aes(x = x, y = y, label = percent(as.numeric(q4[2,2]/(q1[1,2]+q1[2,2])))), inherit.aes = F, color = color.cat) +
  geom_text(size = 4, data = data.frame(x = .085, y = .25), aes(x = x, y = y, label = percent(as.numeric(q5[2,2]/(q5[1,2]+q5[2,2])))), inherit.aes = F, color = color.cathr) +
  geom_text(size = 4, data = data.frame(x = .23, y = .25), aes(x = x, y = y, label = percent(as.numeric(q6[2,2]/(q5[1,2]+q5[2,2])))), inherit.aes = F, color = color.cathr) +
  geom_text(size = 4, data = data.frame(x = .23, y = .055), aes(x = x, y = y, label = percent(as.numeric(q7[2,2]/(q5[1,2]+q5[2,2])))), inherit.aes = F, color = color.cathr) +
  geom_text(size = 4, data = data.frame(x = 0.085, y = 0.055), aes(x = x, y = y, label = percent(as.numeric(q8[2,2]/(q5[1,2]+q5[2,2])))), inherit.aes = F, color = color.cathr) 
