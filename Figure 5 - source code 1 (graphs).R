library(tidyverse)
library(ggpubr)
library(cowplot)
library(Hmisc)
library(scales)
library(ggrepel)


###read cellprofiler output files and convert metadata to factor ####

#start aqp5 quant - pipeline is dapi centered and same as spa otherwise. 2 wk and 12 wk were done at same time, so can include Week in analysis and do ANOVA

aqp5tumors <- as_tibble(read.csv("~/Dropbox/JEV/CellProfiler/aqp5tumorsCytoplasm.csv"))
aqp5tumors$Metadata_Geno <- as.factor(aqp5tumors$Metadata_Geno)
aqp5tumors$Metadata_Mouse <- as.factor(aqp5tumors$Metadata_Mouse)
aqp5tumors$Metadata_Week <- as.factor(aqp5tumors$Metadata_Week)

aqp5tumors2wk <- aqp5tumors %>% filter(Metadata_Week == 2)
aqp5tumors12wk <- aqp5tumors %>% filter(Metadata_Week == 12)

aqp5tumors <- unite(aqp5tumors, col = "WeekXGeno", c("Metadata_Week", "Metadata_Geno"), sep = ".")
aqp5tumors$WeekXGeno <- as.factor(aqp5tumors$WeekXGeno)
aqp5tumors$WeekXGeno <- ordered(aqp5tumors$WeekXGeno, c("2.CAT", "2.CATHR", "12.CAT", "12.CATHR"))

#more aqp5 quant - now aqp5 and lys at same time to look for any correlation

lysaqp5tumors <- as_tibble(read.csv("~/Dropbox/JEV/CellProfiler/aqp5lystumorsCytoplasm.csv"))
lysaqp5tumors$Metadata_Geno <- as.factor(lysaqp5tumors$Metadata_Geno)
lysaqp5tumors$Metadata_Mouse <- as.factor(lysaqp5tumors$Metadata_Mouse)
lysaqp5tumors$Metadata_Week <- as.factor(lysaqp5tumors$Metadata_Week)


###statistics - treat tumors as independent replicates. aggregate individual cell data by mouse, tumor, and week. ####

z <- aqp5tumors


#aqp5 is a little different because quant all done simultaneuoulsy so 4 graphs together
aggaqp5tumors <- aqp5tumors %>% group_by(Metadata_Mouse, Metadata_Tumor, WeekXGeno) %>% summarise(mean_intensity = mean(Intensity_MeanIntensity_AQP5))
aggaqp5tumors.aov <- aov(mean_intensity ~ WeekXGeno, data = aggaqp5tumors)
summary.aov(aggaqp5tumors.aov)
TukeyHSD(aggaqp5tumors.aov)


#violin graphs

color.cat <- "#00BFC4"
color.cathr <- "#F8766D"

p2 <- ggplot(z, aes(x = WeekXGeno, y = Intensity_MedianIntensity_AQP5, fill = WeekXGeno)) + 
  geom_violin(trim = T, bw = .01) +
  stat_summary(size = 1, fun.data = mean_sdl, geom = "pointrange", color = "black", fun.args = list(mult = 1), show.legend = F) +
  coord_cartesian(ylim = c(0,.8)) + 
  scale_fill_manual(values = c(color.cat, color.cathr, color.cat, color.cathr)) +
  labs(fill = NULL, y = "Median Fluorescent Intensity", title = "aqp5 tumors, anova ***, 12cat - 2cat p = .0013576. worse ns, better *** < 1e-5")

ggsave("~/Dropbox/JEV/CAT mouse paper/figures for resubmission/aqp5tumors.tiff", plot = p2, width = 9, height = 8, dpi = 320, units = "in")




### two factor analyses ####

color.cat <- "blue"
color.cathr <- "red"

#test aqp5 vs lyz in combined analyses


p3 <- ggplot(data = filter(lysaqp5tumors, Metadata_Geno == "CAT"), aes(x = Intensity_MedianIntensity_LYS, y = Intensity_MedianIntensity_AQP5)) + 
  stat_binhex(bins = 50) + 
  xlim(0, .5) +
  ylim(0, .5) +
  labs(title = "LYS vs AQP5 staining in CAT tumors", subtitle = "Spearman's rank correlation Rho = .54")


p4 <- ggplot(data = filter(lysaqp5tumors, Metadata_Geno == "CATHR"), aes(x = Intensity_MedianIntensity_LYS, y = Intensity_MedianIntensity_AQP5)) + 
  stat_binhex(bins = 50) + 
  xlim(0, .5) +
  ylim(0, .5) +
  labs(title = "LYS vs AQP5 staining in CAT;PI3K tumors", subtitle = "Spearman's rank correlation Rho = .13")


ggsave("~/Dropbox/JEV/CAT mouse paper/figures for resubmission/lysaqp5cat.tiff", plot = p3, width = 8, height = 8, dpi = 320, units = "in")
ggsave("~/Dropbox/JEV/CAT mouse paper/figures for resubmission/lysaqp5cathr.tiff", plot = p4, width = 8, height = 8, dpi = 320, units = "in")



#define quadrants by mean of CAT tumors minus one SD

lysaqp5tumors.cat <- filter(lysaqp5tumors, Metadata_Geno == "CAT")

cat.lys.line <- (mean(lysaqp5tumors.cat$Intensity_MedianIntensity_LYS) - sd(lysaqp5tumors.cat$Intensity_MedianIntensity_LYS))
cat.aqp5.line <- (mean(lysaqp5tumors.cat$Intensity_MedianIntensity_AQP5) - sd(lysaqp5tumors.cat$Intensity_MedianIntensity_AQP5))


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


p1 <- ggplot(data = lysaqp5tumors, aes(x = Intensity_MedianIntensity_LYS, y = Intensity_MedianIntensity_AQP5)) + 
  stat_binhex(bins = 50, data = filter(lysaqp5tumors, Metadata_Geno == "CAT"), fill = color.cat, aes(alpha = ..count..)) +
  stat_binhex(bins = 50, data = filter(lysaqp5tumors, Metadata_Geno == "CATHR"), fill = color.cathr, aes(alpha = ..count..)) +
  geom_hline(yintercept = cat.lys.line, linetype="dashed", color = "black") +
  geom_vline(xintercept = cat.aqp5.line, linetype="dashed", color = "black") + 
  xlim(0, .5) + 
  ylim(0, .5) +
  geom_text(size = 4, data = data.frame(x = .01, y = .5), aes(x = x, y = y, label = percent(as.numeric(q1[2,2]/(q1[1,2]+q1[2,2])))), inherit.aes = F, color = color.cat) +
  geom_text(size = 4, data = data.frame(x = .45, y = .5), aes(x = x, y = y, label = percent(as.numeric(q2[2,2]/(q1[1,2]+q1[2,2])))), inherit.aes = F, color = color.cat) +
  geom_text(size = 4, data = data.frame(x = .45, y = .035), aes(x = x, y = y, label = percent(as.numeric(q3[2,2]/(q1[1,2]+q1[2,2])))), inherit.aes = F, color = color.cat) +
  geom_text(size = 4, data = data.frame(x = 0.01, y = .025), aes(x = x, y = y, label = percent(as.numeric(q4[2,2]/(q1[1,2]+q1[2,2])))), inherit.aes = F, color = color.cat) +
  geom_text(size = 4, data = data.frame(x = .01, y = .475), aes(x = x, y = y, label = percent(as.numeric(q5[2,2]/(q5[1,2]+q5[2,2])))), inherit.aes = F, color = color.cathr) +
  geom_text(size = 4, data = data.frame(x = .45, y = .475), aes(x = x, y = y, label = percent(as.numeric(q6[2,2]/(q5[1,2]+q5[2,2])))), inherit.aes = F, color = color.cathr) +
  geom_text(size = 4, data = data.frame(x = .45, y = .01), aes(x = x, y = y, label = percent(as.numeric(q7[2,2]/(q5[1,2]+q5[2,2])))), inherit.aes = F, color = color.cathr) +
  geom_text(size = 4, data = data.frame(x = 0.01, y = 0), aes(x = x, y = y, label = percent(as.numeric(q8[2,2]/(q5[1,2]+q5[2,2])))), inherit.aes = F, color = color.cathr) +
  labs(title = "chisq = 7314.3 p val ****")
ggsave("~/Dropbox/JEV/CAT mouse paper/figures for resubmission/lysaqp5overlay.tiff", plot = p1, width = 6, height = 5, dpi = 320, units = "in")






