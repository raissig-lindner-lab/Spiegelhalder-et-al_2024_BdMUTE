##### R Script for Spiegelhalder et al. 2024

#### Required data can be found at https://github.com/lea-berg/spiegelhalder_2024_BdMUTE
### Data for plots is "MUTE2024_Data.xlsx"
### Li-6800 files corrected by leaf area are in "Li-6800-files_Spiegelhalder-et-al_2024.zip"

### load packages required to run the script
library(tidyverse)
library(readxl)
library(agricolae)
library(MetBrewer)
library(licornetics)
library(ggpubr)
library(ggtext)


### Colourblind-friendly palettes -----------------------------------------------------------------------
## to distinguish between genotypes (Fig. 6, Figs S6-9)
as.vector(met.brewer("Morgenstern", n=8))
# WT = #98768e
# MYM = #b08ba5
# M3GM++ = #c7a2b6
# M3GM+ = #dfbbc8
# M3GM- = #ffc680
# M3GM-- = #ffb178
# MMY = #db8872
# sid = #a56457
c("#98768e", "#b08ba5", "#c7a2b6", "#dfbbc8", "#ffc680", "#ffb178", "#db8872", "#a56457")

## to distinguish between phenotypes (Figs 1, 2, 3, 5, Figs S2, S5, S8)
as.vector(met.brewer("Hokusai3", n=4))
# 2 SC = #95c36e
# 1 SC = #d8d97a
# 0 SC = #5a97c1
# aborted = #295384
c("#95c36e", "#d8d97a", "#5a97c1", "#295384")





#### Fig. 1E, F -------------------------------------------------------------------------------------
### load data
distances <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_Fig1E_F")
distances <- na.omit(distances)
colnames(distances) <- c("genotype", "individual", "ACD_GCD_distance", "ACD_GCD_counts")


### statistics (Two-sided t-test)
t.test(distances$ACD_GCD_distance ~distances$genotype)
# p-value = 0.6112 --> no significant difference between WT and sid

t.test(distances$ACD_GCD_counts ~distances$genotype)
# p-value = 0.237 --> no significant difference between WT and sid

sigacdgcd <- data.frame(genotype=c("wt", "sid"),
                        significances=c("a", "a"))


### plot distance
plot1_E <- ggplot(distances, mapping=aes(x=factor(genotype, levels = c("wt", "sid")), y=ACD_GCD_distance, colour = genotype))+
  geom_boxplot(outlier.shape = NA, colour = "black")+
  geom_text(sigacdgcd, mapping=aes(y=1450, label=paste(significances)), show.legend = F, colour = "black", size = 8/.pt)+
  geom_jitter(width=0.1, height = 0, show.legend = F, size = 2.5, alpha = 0.8)+
  scale_y_continuous(limits = c(0, 1450))+
  scale_x_discrete(labels = c("WT", "*sid*"))+
  scale_colour_manual(values = c("grey", "grey40"))+
  theme_classic()+
  theme(axis.text.x = element_markdown())+
  labs(x=NULL, y="ACD to GCD [µm]")

plot1_E_final <- plot1_E + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8))


### plot counts
plot1_F <- ggplot(distances, mapping=aes(x=factor(genotype, levels = c("wt", "sid")), y=ACD_GCD_counts, colour = genotype))+
  geom_boxplot(outlier.shape = NA, colour = "black")+
  geom_text(sigacdgcd, mapping=aes(y=180, label=paste(significances)), show.legend = F, colour = "black", size = 8/.pt)+
  geom_jitter(width=0.1, height = 0, show.legend = F, size = 2.5, alpha = 0.8)+
  scale_x_discrete(labels = c("WT", "*sid*"))+
  scale_y_continuous(limits = c(0, 180))+
  scale_colour_manual(values = c("grey", "grey40"))+
  theme_classic()+
  theme(axis.text.x = element_markdown())+
  labs(x=NULL, y="Number of cells ACD to GCD")

plot1_F_final <- plot1_F + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8))


### plot both together
ggarrange(plot1_E_final, plot1_F_final, labels = c("E", "F"), font.label = list(size = 12))





#### Fig. 2B, C -------------------------------------------------------------------------------------
### load data
intensity <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_Fig2B_C")


### label the groups for different panels in the plot
int_m3gm <- subset(intensity, genotype != "MYM" & genotype != "MMY")
int_m3gm$group <- "3xGFP"

int_yfp <- subset(intensity, genotype == "MYM" | genotype == "MMY")
## filter out individuals with different laser settings
int_yfp <- subset(int_yfp, individual != 6 & individual !=7 & individual !=8)
int_yfp$group <- "YFP"

intensity <- rbind(int_m3gm, int_yfp)


### average ctcf and sd per line and individual
mean_int<- intensity %>% group_by(line, genotype, group, individual) %>% summarise(means=mean(ctcf, na.rm=T),means_sd=sd(ctcf, na.rm=T))
mean_int_m3gm <- int_m3gm %>% group_by(line, genotype, individual) %>% summarise(means=mean(ctcf, na.rm=T),means_sd=sd(ctcf, na.rm=T))
mean_int_yfp <- int_yfp %>% group_by(line, genotype, individual) %>% summarise(means=mean(ctcf, na.rm=T),means_sd=sd(ctcf, na.rm=T))


### statistical analysis with Anova and Tukey's test
## compare all lines
# ANOVA
anovaint<-aov(mean_int$means ~mean_int$genotype)
summary(anovaint)
# Tukey's test
tukeyint <- HSD.test(anovaint, trt="mean_int$genotype")
tukeyint$groups
sigint<- data.frame(genotype=c("MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY"), 
                    group=c("YFP", "3xGFP", "3xGFP", "3xGFP", "3xGFP", "YFP"),
                    significance=c("a", "a", "a", "b", "b", "a"))

## compare only 3xGFP (M3GM) lines
# ANOVA
anovaintm3gm<-aov(mean_int_m3gm$means ~mean_int_m3gm$genotype)
summary(anovaintm3gm)
# Tukey's test
tukeyintm3gm <- HSD.test(anovaintm3gm, trt="mean_int_m3gm$genotype")
tukeyintm3gm$groups

## compare only YFP (MYM and MMY) lines
# ANOVA
anovaintyfp <-aov(mean_int_yfp$means ~mean_int_yfp$genotype)
summary(anovaintyfp)
# Tukey's test
tukeyintyfp <- HSD.test(anovaintyfp, trt="mean_int_yfp$genotype")
tukeyintyfp$groups


### order by complementation phenotype
mean_int$genotype <- ordered(mean_int$genotype , levels=c("MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY"))
mean_int_m3gm$genotype <- ordered(mean_int_m3gm$genotype , levels=c("M3GM++", "M3GM+", "M3GM-", "M3GM--"))
mean_int_yfp$genotype <- ordered(mean_int_yfp$genotype , levels=c("MYM", "MMY"))


#plot (dots are individuals) - all lines split by fluorophore
plot2_B_C <- ggplot(mean_int, aes(x=factor(genotype), y=means, colour=genotype))+
  geom_boxplot(outlier.shape = NA, show.legend = F, colour = "black")+
  geom_jitter(shape=16, fill=NA, width=0.1, alpha=0.8, show.legend = F)+
  geom_text(sigint, mapping=aes(y=4000, label=paste(significance)), show.legend = F, colour = "black", size = 8/.pt)+
  scale_colour_manual(values=c("M3GM++" = "#282c29", "M3GM+" = "#6d756f", "M3GM-" = "#b2beb5","M3GM--" = "#cdd6d0",
                               "MYM" = "#282c29", "MMY" = "#cdd6d0"))+
  scale_y_continuous(limits = c(0, 4000))+
  scale_x_discrete(labels=c("M3GM++" = "M3GM++", "M3GM+" = "M3GM+", "M3GM-" = "M3GM-","M3GM--" = "M3GM- -",
                            "MYM" = "MYM", "MMY" = "MMY"))+
  facet_grid(~group, scales = "free", space = "free")+
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 4), 
                             title = "Line"))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-10)),
        legend.background = element_rect(fill=NA))+
  theme(axis.text.x= element_markdown())+
  labs(x=NULL, y="Mean cell fluorescence (CTCF)")

plot2_B_C <- plot2_B_C + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))





#### Fig. 2E, Fig. S2 -------------------------------------------------------------------------------------
### load data
pheno3rd <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_Fig2E_FigS2")
pheno3rd <- na.omit(pheno3rd)


### set colour palette
colourblindpalette <- c("#295384", "#5a97c1", "#d8d97a", "#95c36e")


### Fig. 2E
## calculate means for all types for each line
means1<- pheno3rd %>% group_by(type) %>% 
  summarize(wt_mean=mean(`wt`),
            wt_sd=sd(`wt`),
            one_SC=mean(`1sub`),
            one_SC_sd=sd(`1sub`),
            no_SC=mean(`0sub`),
            no_SC_sd=sd(`0sub`),
            aborted=mean(abort),
            aborted_sd=sd(abort),
            pictures=n())

## reorder columns to have one column 'type' (contains name of phenotype) and another with the percentages ('percent') of the phenotypes
means2<- means1 %>% select(-pictures) %>%
  pivot_longer(., cols = c(wt_mean, `one_SC`, `no_SC`, aborted), names_to = "phenotype", values_to = "percent")

## plot phenotypes displayed as stacked plot (types in one bar)
plot2_E <- means2 %>% mutate(phenotype = fct_relevel(phenotype, "aborted", "no_SC", "one_SC", "wt_mean")) %>%
  ggplot(aes(x=factor(type, levels = c("wt", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=percent, fill=phenotype))+
  geom_col(position="fill", show.legend = T)+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 25, 50, 75, 100))+
  scale_fill_manual(labels=c("wt_mean"="WT", "one_SC"="1 SC", "no_SC"="0 SC", "aborted"="aborted"),
                    values=colourblindpalette)+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM- -", "MMY", "*sid*"))+
  guides(fill = guide_legend(override.aes = list(shape = 22), 
                             title = "Phenotype",
                             reverse = TRUE))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_markdown())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-16)),
        legend.background = element_rect(fill=NA))+
  labs(x=NULL,y="Percentage [%]",title=NULL)

plot2_E_final <- plot2_E + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))


### Fig. S2
## reorder data frame
phenopercentages <- data.frame(line=pheno3rd$type, identifier=pheno3rd$line, individual=pheno3rd$individual, phenotype="WT", percent=pheno3rd$wt)
phenopercentages <- rbind(phenopercentages, data.frame(line=pheno3rd$type, identifier=pheno3rd$line, individual=pheno3rd$individual, phenotype="one_SC", percent=pheno3rd$`1sub`))
phenopercentages <- rbind(phenopercentages, data.frame(line=pheno3rd$type, identifier=pheno3rd$line, individual=pheno3rd$individual, phenotype="no_SC", percent=pheno3rd$`0sub`))
phenopercentages <- rbind(phenopercentages, data.frame(line=pheno3rd$type, identifier=pheno3rd$line, individual=pheno3rd$individual, phenotype="aborted", percent=pheno3rd$abort))

plotS2 <- ggplot(phenopercentages, aes(x=factor(line, levels = c("wt", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=percent, colour=factor(phenotype, levels = c("WT", "one_SC", "no_SC", "aborted"))))+
  geom_boxplot(position="dodge", show.legend = F, outlier.shape = NA)+
  geom_point(position = position_dodge(width = 0.75), show.legend = T, alpha = 0.8)+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 25, 50, 75, 100))+
  scale_colour_manual(labels=c("wt"="WT", "one_SC"="1 SC", "no_SC"="0 SC", 
                               "aborted"="aborted"),
                      values=rev(colourblindpalette))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM- -", "MMY", "*sid*"))+
  guides(colour = guide_legend(title = "Phenotype"))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_markdown())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-16)),
        legend.background = element_rect(fill=NA))+
  labs(x=NULL, y="Percentage [%]",title=NULL)

plotS2final <- plotS2 + theme(text = element_text(size = 8),
                                      axis.text = element_text(size = 8),
                                      strip.text = element_text(size = 8),
                                      legend.text = element_text(size = 8))





#### Fig. 3B -------------------------------------------------------------------------------------
### load data
mobility <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_Fig3B")


### calculate mobility ratio
mobility$SMC_means <- rowMeans(mobility[ , c(20,28)], na.rm=TRUE)
mobility$SMC_means <- gsub("NaN", NA, mobility$SMC_means)
mobility$SMC_means <- as.numeric(mobility$SMC_means)

mobility$ratio <-mobility$SMC_means/mobility$`GMC intden`


### statistical test
anovamob <- aov(mobility$ratio ~mobility$geno)
summary(anovamob)

tukeymob <- HSD.test(anovamob, trt = "mobility$geno")
tukeymob$groups

sigmob <- data.frame(geno=c("MYM", "M3GM++"),
                     significance = c("a", "b"))


#plot 
plot3_B <- ggplot(mobility,aes(x=factor(geno, levels = c("MYM", "M3GM++")), y=ratio, colour=geno))+
  geom_boxplot(outlier.shape = NA, show.legend = F, colour="black")+
  geom_jitter(width=0.1, height=0, alpha=0.7, show.legend = F)+
  geom_text(sigmob, mapping=aes(y=1, label=paste(significance)), show.legend = F, colour="black", size = 8/.pt)+
  scale_y_continuous(limits = c(0,1))+
  scale_colour_manual(values = c("grey", "grey40"))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  labs(x=NULL, y=" SMC nuclear  / GMC cytoplasmic ")

plot3_B_final <- plot3_B + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))





#### Fig. 3D -------------------------------------------------------------------------------------
### load data
screc <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_Fig3D")
screc$number_of_SC <- as.character(screc$number_of_SC)
screc2 <- subset(screc, !is.na(screc$new_ident))

### modify facet labels
facetlabels <- c('MYM' = "MYM", '3xGFP-MUTE' = "M3GM++")

### plot
plot3_D <- ggplot(screc2, aes(x=factor(new_ident), y=distance_from_first_GMC_division, colour=factor(number_of_SC), shape=factor(divided)))+
  geom_point(size=1.5)+
  facet_grid(~factor(genotype, levels = c("MYM", "3xGFP-MUTE")),
             labeller = as_labeller(facetlabels))+
  #scale_shape_manual(values = c("\u25A1", "\u25E7"))+
  scale_shape_manual(values = c(1, 16))+
  scale_colour_manual(values = c("2" = "#95d36e", "1" = "#d8d97a", "0" = "#5a97c1"))+
  theme_classic()+
  theme(legend.position = "bottom",
        legend.box = "vertical")+
  guides(shape=guide_legend(title="Divided GMC?"),
         colour=guide_legend(title="Number of SCs"))+
  labs(x="Individual", y="Distance from first GMC division")
# save as 400x400

plot3_D_final <- plot3_D + theme(text = element_text(size = 8),
                                   axis.text = element_text(size = 8),
                                   strip.text = element_text(size = 8),
                                   legend.text = element_text(size = 8))





#### Fig. 3E, Fig. S3 -------------------------------------------------------------------------------------
### load data
gc1 <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_Fig3E_FigS3")
gc1$SC <- as.character(gc1$SC)
gc1$LWratio <- as.numeric(gc1$LWratio)

### Fig. 3E
## statistical analysis
## ANOVA
intergc <- interaction(gc1$Genotype, gc1$SC)
anovagc <- aov(gc1$LWratio ~intergc)
summary(anovagc)
## Tukey's test
tukeygc <- HSD.test(anovagc, trt="intergc")
tukeygc$groups
siggc<- data.frame(Genotype=c("MYM", "MYM", "MYM", "M3GM++", "M3GM++", "M3GM++"), 
                   SC=c("0", "1", "2", "0", "1", "2"),
                   significance=c("a", "ab", "d", "a", "bc", "c"))

## plot
plot3_E <- ggplot(gc1,aes(x=SC,y=LWratio, fill = factor(Genotype, levels = c("MYM", "M3GM++")), colour = factor(Genotype, levels = c("MYM", "M3GM++"))))+
  geom_boxplot(outlier.shape = NA, show.legend = F, colour = "black")+
  geom_point(alpha=0.7, position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), show.legend = T)+
  geom_text(siggc, mapping = aes(y=4, label=paste(significance)), show.legend = F, colour = "black", size = 8/.pt, position = position_dodge(width=0.7))+
  scale_colour_manual(values = c("grey40", "grey"))+
  scale_fill_manual(values = c("white", "white"), guide = "none")+
  scale_y_continuous(limits = c(0, 4))+
  theme_bw()+
  theme(axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-10)),
        legend.background = element_rect(fill=NA))+
  guides(colour = guide_legend(title = "Genotype"))+
  labs(x="Number of SCs", y="Length / width ratio of GMC",title=NULL)
# save as 350x350
plot3_E_final <- plot3_E + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))


### Fig. S3
## statistical analysis
# ANOVA
interlwgmc <- interaction(gc1$Genotype, gc1$Divided)
anovalwgmc <- aov(gc1$LWratio ~interlwgmc)
summary(anovalwgmc)
# Tukey's test
tukeylwgmc <- HSD.test(anovalwgmc, trt = "interlwgmc")
tukeylwgmc$groups

siglwgmc <- data.frame(Genotype=c("MYM", "M3GM++", "MYM", "M3GM++", "MYM", "M3GM++"),
                       Divided=c("yes", "yes", "dividing", "dividing", "no", "no"),
                       significance = c("a", "b", "ab", "bc", "c", "c"))

## plot
plotS3 <- ggplot(gc1,aes(x=factor(Divided, levels = c("no", "dividing", "yes")), y=LWratio, fill = factor(Genotype, levels = c("MYM", "M3GM++")), colour = factor(Genotype, levels = c("MYM", "M3GM++"))))+
  geom_boxplot(outlier.shape = NA, show.legend = F, colour = "black")+
  geom_point(alpha=0.7, position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), show.legend = T)+
  geom_text(siglwgmc, mapping=aes(y=4, label=paste(significance)), show.legend = F, colour="black", size = 8/.pt, position = position_dodge(width=0.7))+
  scale_colour_manual(values = c("grey40", "grey"))+
  scale_fill_manual(values = c("white", "white"), guide = "none")+
  scale_y_continuous(limits = c(0, 4))+
  theme_bw()+
  theme(axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-10)),
        legend.background = element_rect(fill=NA))+
  guides(colour = guide_legend(title = "Genotype"))+
  labs(x="Divided?", y="Length / width ratio of GMC",title=NULL)
# save as 350x350
plotS3_final <- plotS3 + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))





#### Fig. 4C, Fig. S4A, C -------------------------------------------------------------------------------------
### load data
division_data <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_Fig4C_FigS4A_C")


### exclude complexes that are not divided
division <- division_data %>% filter(!transverse=="undivided")


### calculate ratios
division$distance_left_top <- as.numeric(division$distance_left_top)
division$distance_right_top <- as.numeric(division$distance_right_top)
division$distance_left_bottom <- as.numeric(division$distance_left_bottom)
division$distance_right_bottom <- as.numeric(division$distance_right_bottom)

division$top_ratio <- NA
division$bottom_ratio <- NA

for(row in 1:length(division$distance_left_top)) {
  if(!is.na(division$distance_left_top[row])){
    if(division$distance_left_top[row] > division$distance_right_top[row]) {
      division$top_ratio[row] <- division$distance_right_top[row]/division$distance_left_top[row]
    }
    else {
      division$top_ratio[row] <- division$distance_left_top[row]/division$distance_right_top[row]
    }
  }
  if(!is.na(division$distance_left_bottom[row])){
    if(division$distance_left_bottom[row] > division$distance_right_bottom[row]) {
      division$bottom_ratio[row] <- division$distance_right_bottom[row]/division$distance_left_bottom[row]
    }
    else {
      division$bottom_ratio[row] <- division$distance_left_bottom[row]/division$distance_right_bottom[row]
    }
  }
} 

for(row in 1:length(division$top_ratio)) {
  if(is.nan(division$top_ratio[row])) {
    division$top_ratio[row] <- 0
  }
  
  if(is.nan(division$bottom_ratio[row])) {
    division$bottom_ratio[row] <- 0
  }
}

## statistical analysis
# top ratio
anovadiv_top <- aov(division$top_ratio ~division$line)
tukeydiv_top <- HSD.test(anovadiv_top, trt="division$line")
tukeydiv_top$groups

# bottom ratio
anovadiv_bot <- aov(division$bottom_ratio ~division$line)
tukeydiv_bot <- HSD.test(anovadiv_bot, trt="division$line")
tukeydiv_bot$groups

sig_div <- data.frame(line = c("WT", "M3GM--", "MMY", "sid"),
                      sig_top = c("a", "a", "b", "c"),
                      sig_bottom = c("a", "a", "b", "c"))

## set categories for the plots
division$classes_top <- NA

for(row in 1:length(division$top_ratio)) {
  if(division$top_ratio[row] > 0.75) {
    division$classes_top[row] <- "Above 0.75"
  }
  
  else {
    division$classes_top[row] <- "Below 0.75"
  }
}

division$classes_bottom <- NA

for(row in 1:length(division$bottom_ratio)) {
  if(division$bottom_ratio[row] > 0.75) {
    division$classes_bottom[row] <- "Above 0.75"
  }
  
  else {
    division$classes_bottom[row] <- "Below 0.75"
  }
}

## plot the ratio of the distances between division plane and stomatal borders at the top apex (Fig. S4A)
plotS4_A <- ggplot(division, aes(x=factor(line, levels = c("WT", "M3GM--", "MMY", "sid")), colour=classes_top, y=top_ratio))+
  geom_violin(colour="black")+
  geom_jitter(width=0.1, height=0, alpha=0.5, show.legend = F)+
  geom_text(sig_div, mapping = aes(y=1.05, label=paste(sig_top)), show.legend = F, colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0,1.05))+
  scale_colour_manual(values=c("Above 0.75"="grey40", "Below 0.75"="red"))+
  scale_x_discrete(labels = c("WT", "M3GM- -", "MMY", "*sid*"))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_markdown())+
  #labs(x=NULL,y="Ratio",title="Ratio division plane to stomatal borders (top)")
  labs(x=NULL,y="Ratio division plane to stomatal borders (top)",title=NULL)

plotS4_A_final <- plotS4_A + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))

## plot the ratio of the distances between division plane and stomatal borders at the bottom apex (Fig. S4C)
plotS4_C <- ggplot(division, aes(x=factor(line, levels = c("WT", "M3GM--", "MMY", "sid")), colour=classes_bottom, y=bottom_ratio))+
  geom_violin(colour="black")+
  geom_jitter(width=0.1, height=0, alpha=0.5, show.legend = F)+
  geom_text(sig_div, mapping = aes(y=1.05, label=paste(sig_bottom)), show.legend = F, colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0,1.05))+
  scale_colour_manual(values=c("Above 0.75"="grey40", "Below 0.75"="red"))+
  scale_x_discrete(labels = c("WT", "M3GM- -", "MMY", "*sid*"))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_markdown())+
  #labs(x=NULL,y="Ratio",title="Ratio division plane to stomatal borders (bottom)")
  labs(x=NULL,y="Ratio division plane to stomatal borders (bottom)",title=NULL)


plotS4_C_final <- plotS4_C + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))

## plot together for Fig. S4
ggarrange(plotS4_A_final, plotS4_C_final, labels = c("A", "C"), font.label = list(size = 12))




### Calculate combined ratio for "total skewness" - Fig. 4C
division$combined_ratio <- division$top_ratio + division$bottom_ratio

## statistical analysis
anovadiv_com <- aov(division$combined_ratio ~division$line)
tukeydiv_com <- HSD.test(anovadiv_com, trt="division$line")
tukeydiv_com$groups

sig_div2 <- data.frame(line = c("WT", "M3GM--", "MMY", "sid"),
                      sig_comb = c("a", "a", "b", "c"))

## add classification cutoff ratio at 1.6
division$classes <- NA

for(row in 1:length(division$combined_ratio)) {
  if(division$combined_ratio[row] > 1.6) {
    division$classes[row] <- "Above 1.6"
  }
  
  else {
    division$classes[row] <- "Below 1.6"
  }
}

## plot combined ratio (Fig. 4, panel C)
plot4_C <- ggplot(division, aes(x=factor(line, levels = c("WT", "M3GM--", "MMY", "sid")), colour=classes, y=combined_ratio))+
  geom_violin(colour="black")+
  geom_jitter(width=0.1, height=0, alpha=0.5, show.legend = F)+
  geom_hline(yintercept=1.6, colour = "black", linetype = "dashed")+
  geom_text(sig_div2, mapping = aes(y=2.10, label=paste(sig_comb)), show.legend = F, colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0,2.10))+
  scale_colour_manual(values=c("Above 1.6"="grey40", "Below 1.6"="red"))+
  scale_x_discrete(labels = c("WT", "M3GM- -", "MMY", "*sid*"))+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_markdown(),
        text = element_text(size = 8),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8))+
  #labs(x=NULL,y="Ratio",title="Ratio division plane to stomatal borders (bottom)")
  labs(x=NULL,y="Combined ratio of division plane to stomatal borders",title=NULL)





#### Fig. 5B, Fig. S5 -------------------------------------------------------------------------------------
### load data
morpho <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_Fig5B_FigS5")
morpho <- morpho %>% filter(complex_nr != "NA")
morpho$gc_with_sc <- as.character(morpho$gc_with_sc)
morpho$gc_width_top <- as.numeric(morpho$gc_width_top)
morpho$gc_width_bottom <- as.numeric(morpho$gc_width_bottom)
morpho$gc_length <- as.numeric(morpho$gc_length)
morpho$gc_width_middle_1 <- as.numeric(morpho$gc_width_middle_1)


### Fig. 5B
## calculate average at the apex
morpho$apex_average <- rowMeans(morpho[, c("gc_width_top", "gc_width_bottom")], na.rm=TRUE)
morpho$apex_average <- gsub("NaN", NA, morpho$apex_average)
morpho$apex_average <- as.numeric(morpho$apex_average)


## calculate ratio of GC width at the apices to width at the middle
morpho$apex_middle_ratio <- morpho$apex_average/morpho$gc_width_middle_1

## statistical analysis
# ANOVA
gcinter <- interaction(morpho$sc_count,morpho$gc_with_sc)
anovagcw <- aov(morpho$apex_middle_ratio ~gcinter)
summary(anovagcw)
# Tukey's test
tukeygcw <- HSD.test(anovagcw, trt="gcinter")
tukeygcw

siggcw <- data.frame(sc_count=c(0,1,1,2), 
                    gc_with_sc=c("no", "no", "yes", "yes"),
                    significance=c("a", "a", "b", "c"))

## plot - Fig. 5B
facetlabels <- c('0' = "0 SC",
                 '1' = "1 SC",
                 '2' = "2 SC")

plot5_B <- ggplot(morpho,aes(x=gc_with_sc, y=apex_middle_ratio, colour=factor(sc_count)))+
  geom_boxplot(outlier.shape=NA, show.legend = F, colour = "black")+
  geom_jitter(width=0.1, alpha=0.7, show.legend = F)+
  geom_text(siggcw, mapping=aes(y=3.5, label=paste(significance)), show.legend = F, colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0, 3.5))+
  scale_colour_manual(values = c("#5a97c1", "#d8d97a", "#95d36e"))+
  facet_grid(~sc_count, labeller=as_labeller(facetlabels), scales = "free", space = "free")+
  guides(colour = guide_legend(title = "Guard cell with subsidiary cell?"))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-16)),
        legend.background = element_rect(fill=NA))+
  labs(x="Guard cell with subsidiary cell?",y="Guard cell width at apices / center")

plot5_B_final <- plot5_B + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))


### Fig. S5
## guard cell length in the complexes is the same for both guard cells so i will only plot the length for one GC
morphoshort1 <- morpho %>% filter(gc_with_sc == "no")
morphoshort2 <- morpho %>% filter(gc_with_sc == "yes" & sc_count == 2)
morphoshort <- rbind(morphoshort1, morphoshort2)

## statistical analysis
#ANOVA
anovagcl <- aov(morphoshort$gc_length ~morphoshort$sc_count)
summary (anovagcl)
#Tukey's test
tukeygcl <- HSD.test(anovagcl, trt="morphoshort$sc_count")
tukeygcl$groups
siggcl <- data.frame(sc_count=c(0,1,2), 
                     significance=c("a", "b", "c"))

## plot - Fig. S5
plotS5 <- ggplot(morpho,aes(x=factor(sc_count), y=gc_length, colour=factor(sc_count)))+
  geom_boxplot(outlier.shape=NA, show.legend = F, colour = "black")+
  geom_jitter(width=0.1, alpha=0.7, show.legend = F)+
  geom_text(siggcl, mapping=aes(y=33, label=paste(significance)), show.legend = F, colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0, 33))+
  scale_colour_manual(values = c("#5a97c1", "#d8d97a", "#95d36e"))+
  guides(colour = guide_legend(title = "Guard cell with subsidiary cell?"))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-16)),
        legend.background = element_rect(fill=NA))+
  labs(x="Number of SCs",y="Guard cell length [µm]")

plotS5_final <- plotS5 + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))





#### Fig. 5D -------------------------------------------------------------------------------------
### load data
mor3d <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_Fig5D")
mor3d <- na.omit(mor3d)
mor3d$area_apex1 <- as.numeric(mor3d$area_apex1)
mor3d$area_apex2 <- as.numeric(mor3d$area_apex2)
mor3d$area_middle <- as.numeric(mor3d$area_middle)


### calculate average of area at the apices and the ratio of area at the apex/middle
mor3d$apex_average <- rowMeans(mor3d[, c("area_apex1", "area_apex2")], na.rm=TRUE)
mor3d$apex_average <- as.numeric(mor3d$apex_average)

mor3d$apex_middle_ratio <-mor3d$apex_average/mor3d$area_middle


### statistical analysis
## ANOVA
gcinter<-interaction(mor3d$sc_count,mor3d$gc_with_sc)
anovamor3d<-aov(mor3d$apex_middle_ratio ~gcinter)
summary(anovamor3d)
## Tukey's test
tukeymor3d <- HSD.test(anovamor3d, trt="gcinter")
tukeymor3d$groups
sigmor3d <- data.frame(sc_count=c(0,1,1,2), 
                       gc_with_sc=c("no", "no", "yes", "yes"),
                       significance=c("a", "a", "ab", "b"))


### plot 
facetlabels <- c('0' = "0 SC",
                 '1' = "1 SC",
                 '2' = "2 SC")

plot5_D <- ggplot(mor3d,aes(x=gc_with_sc, y=apex_middle_ratio, colour=factor(sc_count)))+
  geom_boxplot(outlier.shape=NA, show.legend = F, colour = "black")+
  geom_jitter(width=0.1, alpha=0.7, show.legend = F)+
  geom_text(sigmor3d, mapping=aes(y=4.5, label=paste(significance)), show.legend = F, colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0, 4.5))+
  scale_colour_manual(values = c("#5a97c1", "#d8d97a", "#95d36e"))+
  facet_grid(~sc_count, labeller=as_labeller(facetlabels), scales = "free", space = "free")+
  guides(colour = guide_legend(title = "Guard cell with subsidiary cell?"))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-16)),
        legend.background = element_rect(fill=NA))+
  labs(x="Guard cell with subsidiary cell?",y="Guard cell area at apices / center")

plot5_D_final <- plot5_D + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))





#### Fig. 6, Figs S6, S7, S9 -------------------------------------------------------------------------------------
### use licornetics package on the Li-6800 excel files (manually corrected by leaf area)

### Fig. 6 (physiology not corrected by stomatal density)
## Plot 6A: Stomatal conductance (gsw) of WT, MYM, MMY and sid
plot6_A <- licorplots(c("wt", "mym", "mmy", "sid"), timestamps = c(20,40,60), timeframe = c(11:75), 
                    legend_labels = c("WT", "MYM", "MMY", "*sid*"), legend_title = NULL,
                    remove_outliers = "yes",
                    type = "gsw", colours = c("#98768e", "#b08ba5", "#db8872", "#a56457"))

## Plot 6B: Stomatal conductance (gsw) of WT, M3GM (++, +, -, --) and sid
plot6_B <- licorplots(c("wt", "LB37", "LB14", "m3gm-LB3", "LB2", "sid"), timestamps = c(20,40,60), timeframe = c(11:75), 
                    legend_labels = c("WT", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "*sid*"), legend_title = NULL,
                    remove_outliers = "yes",
                    type = "gsw", colours = c("#98768e", "#c7a2b6", "#dfbbc8", "#ffc680", "#ffb178", "#a56457"))

## Plot 6C: Carbon assimilation (A) of WT, MYM, MMY and sid
plot6_C <- licorplots(c("wt", "mym", "mmy", "sid"), timestamps = c(20,40,60), timeframe = c(11:59), 
                      legend_labels = c("WT", "MYM", "MMY", "*sid*"), legend_title = NULL,
                      remove_outliers = "yes",
                      type = "A", colours = c("#98768e", "#b08ba5", "#db8872", "#a56457"))

## Plot 6D: Carbon assimilation (A) of WT, M3GM (++, +, -, --) and sid
plot6_D <- licorplots(c("wt", "LB37", "LB14", "m3gm-LB3", "LB2", "sid"), timestamps = c(20,40, 60), timeframe = c(11:59), 
                    legend_labels = c("WT", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "*sid*"), legend_title = NULL,
                    remove_outliers = "yes",
                    type = "A", colours = c("#98768e", "#c7a2b6", "#dfbbc8", "#ffc680", "#ffb178", "#a56457"))

## Plot 6E: Intrinsic water-use efficiency (iWUE) of WT, MYM, MMY and sid
plot6_E <- licorplots(c("wt", "mym", "mmy", "sid"), timestamps = c(20,40, 60), timeframe = c(11:59), 
                    legend_labels = c("WT", "MYM", "MMY", "*sid*"), legend_title = NULL,
                    remove_outliers = "yes",
                    type = "WUE", colours = c("#98768e", "#b08ba5", "#db8872", "#a56457"))

## Plot 6F: Intrinsic water-use efficiency (iWUE) of WT, M3GM (++, +, -, --) and sid
plot6_F <- licorplots(c("wt", "LB37", "LB14", "m3gm-LB3", "LB2", "sid"), timestamps = c(20,40,60), timeframe = c(11:59), 
                    legend_labels = c("WT", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "*sid*"), legend_title = NULL,
                    remove_outliers = "yes",
                    type = "WUE", colours = c("#98768e", "#c7a2b6", "#dfbbc8", "#ffc680", "#ffb178", "#a56457"))

## correct text sizes and assemble Fig. 6
plot6_A_final <- plot6_A + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8))
plot6_B_final <- plot6_B + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8))
plot6_C_final <- plot6_C + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8))
plot6_D_final <- plot6_D + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8))
plot6_E_final <- plot6_E + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8))
plot6_F_final <- plot6_F + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8))

ggarrange(ggarrange(plot6_A_final, plot6_C_final, plot6_E_final, nrow = 3, 
                    common.legend = T, legend = "none", labels = c("A", "C", "E"), font.label = list(size = 12), align = "v"), 
          ggarrange(plot6_B_final, plot6_D_final, plot6_F_final, nrow = 3, 
                    common.legend = T, legend = "none", labels = c("B", "D", "F"), font.label = list(size = 12), align = "v"), 
          ncol=2)



### Fig. S6 and Fig. S7
## Calculate means and sd of steady state gsw, A and iWUE per line per light intensity (for Table 1)
# Low light = 21:38
# High light = 41:58
# Darkness = 61:75
gasexchange <- licorvalues(c("wt", "mym", "LB37", "LB14", "m3gm-LB3", "LB2", "mmy", "sid"),
            transition = list(c(21:38), c(41:58), c(61:75)),
            label = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM- -", "MMY", "*sid*"), 
            colours = c("#98768e", "#b08ba5", "#c7a2b6", "#dfbbc8", "#ffc680", "#ffb178", "#db8872", "#a56457"), 
            remove_outliers = "yes")

## add a "genotype" column
unique(gasexchange$ID)
gasexchange$genotype <- c("M3GM+", "M3GM+", "M3GM+", "M3GM--", "M3GM--", "M3GM--",
                          "M3GM++", "M3GM++", "M3GM++", "M3GM-", "M3GM-", "M3GM-",
                          "MMY", "MMY", "MMY", "MYM", "MYM", "MYM",
                          "sid", "sid", "sid", "WT", "WT", "WT")

## For statistics and plots (Fig. S6 and Fig. S7)
individual_data <- licorvalues(c("LB1.1", "LB1.2", "LB1.4", "LB1.5",
                                 "LB8.1", "LB8.2", "LB8.3", "LB8.4", "LB8.5",
                                 "LB11.3", "LB11.7", "LB11.2", "LB11.4", "LB11.5", "LB11.8",
                                 "LB13.1", "LB13.2", "LB13.3", "LB13.4", "LB13.5", "LB13.8",
                                 "LB36.1", "LB36.2", "LB36.3", "LB36.4", "LB36.5",
                                 "LB9.1", "LB9.2", "LB9.3", "LB9.4", "LB9.5",
                                 "LB37.1", "LB37.2", "LB37.3", "LB37.4", "LB37.5",
                                 "LB14.1", "LB14.2", "LB14.3", "LB14.4", "LB14.5", "LB14.8",
                                 "LB3.1", "LB3.2", "LB3.3", "LB3.4", "LB3.5",
                                 "LB2.1", "LB2.2", "LB2.3", "LB2.4", "LB2.5",
                                 "LB10.1", "LB10.2", "LB10.3", "LB10.4", "LB10.5",
                                 "LB12.1", "LB12.2", "LB12.3", "LB12.4", "LB12.5"),
                               transition = list(c(21:38), c(41:58), c(61:75)),
                               remove_outliers = "yes")

## add columns for genotype and line
for(row in 1:length(individual_data$ID)) {
  if(grepl("LB1.", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "WT"
    individual_data$line[row] <- "LB1"
  }
  
  if(grepl("LB8.", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "WT"
    individual_data$line[row] <- "LB8"
  }
  
  if(grepl("LB11", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "WT"
    individual_data$line[row] <- "LB11"
  }
  
  if(grepl("LB13", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "WT"
    individual_data$line[row] <- "LB13"
  }
  
  if(grepl("LB3", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "M3GM-"
    individual_data$line[row] <- "LB3"
  }   
  
  if(grepl("LB36", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "WT"
    individual_data$line[row] <- "LB36"
  }
  
  if(grepl("LB37", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "M3GM++"
    individual_data$line[row] <- "LB37"
  }
  
  if(grepl("LB2", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "M3GM--"
    individual_data$line[row] <- "LB2"
  }
  
  if(grepl("LB9", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "MYM"
    individual_data$line[row] <- "LB9"
  }
  
  if(grepl("LB10", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "MMY"
    individual_data$line[row] <- "LB10"
  }
  
  if(grepl("LB12", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "sid"
    individual_data$line[row] <- "LB12"
  }
  
  if(grepl("LB14", individual_data$ID[row])==T) {
    individual_data$genotype[row] <- "M3GM+"
    individual_data$line[row] <- "LB14"
  }
}

## statistical analysis
sigphys <- data.frame(genotype= c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid"))

## low light (LL)
shortdata_LL <- individual_data %>% filter(transition_zone == "21:38")
# gsw
anova_gsw <- aov(shortdata_LL$gsw ~shortdata_LL$genotype)
summary(anova_gsw)
tukey_gsw <- HSD.test(anova_gsw, trt = "shortdata_LL$genotype")
tukey_gsw$groups
# no significant difference
sigphys$gsw_LL <- c("a") 

# A - LL
anova_A <- aov(shortdata_LL$A ~shortdata_LL$genotype)
summary(anova_A)
tukey_A <- HSD.test(anova_A, trt = "shortdata_LL$genotype")
tukey_A$groups
# no significant difference
sigphys$A_LL <- c("a")

# iWUE - LL
anova_WUE <- aov(shortdata_LL$iWUE ~shortdata_LL$genotype)
summary(anova_WUE)
tukey_WUE <- HSD.test(anova_WUE, trt = "shortdata_LL$genotype")
tukey_WUE$groups
# no significant difference
sigphys$iWUE_LL <- c("a")


## high light
shortdata_HL <- individual_data %>% filter(transition_zone == "41:58")
# gsw - HL
anova_gsw <- aov(shortdata_HL$gsw ~shortdata_HL$genotype)
summary(anova_gsw)
tukey_gsw <- HSD.test(anova_gsw, trt = "shortdata_HL$genotype")
tukey_gsw$groups
#shortdata_HL$gsw groups
#MYM       0.27835001      a
#M3GM++    0.26825945      a
#WT        0.25337923      a
#M3GM+     0.19859154     ab
#M3GM-     0.13988734     bc
#M3GM--    0.10678567     bc
#MMY       0.08386875      c
#sid       0.06949059      c
sigphys$gsw_HL <- c("a", "a", "a", "ab", "bc", "bc", "c", "c")

# A - HL
anova_A <- aov(shortdata_HL$A ~shortdata_HL$genotype)
summary(anova_A)
tukey_A <- HSD.test(anova_A, trt = "shortdata_HL$genotype")
tukey_A$groups
#shortdata_HL$A groups
#M3GM++    20.84342      a
#WT        17.67942     ab
#MYM       17.54528    abc
#M3GM+     16.39369    abc
#M3GM-     13.22839   abcd
#M3GM--    12.66383    bcd
#MMY        9.99713     cd
#sid        8.10210      d
sigphys$A_HL <- c("ab", "abc", "a", "abc", "abcd", "bcd", "cd", "d")

# iWUE - HL
anova_WUE <- aov(shortdata_HL$iWUE ~shortdata_HL$genotype)
summary(anova_WUE)
tukey_WUE <- HSD.test(anova_WUE, trt = "shortdata_HL$genotype")
tukey_WUE$groups
#shortdata_HL$iWUE groups
#MMY         120.29855      a
#sid         119.95831      a
#M3GM--      118.91023      a
#M3GM-        96.29050     ab
#M3GM+        86.42987     bc
#M3GM++       77.90774     bc
#WT           70.92771      c
#MYM          64.16201      c
sigphys$iWUE_HL <- c("a", "a", "ab", "ab", "bc", "c", "c", "c")

## Darkness
shortdata_D <- individual_data %>% filter(transition_zone == "61:75")
# gsw - D
anova_gsw <- aov(shortdata_D$gsw ~shortdata_D$genotype)
summary(anova_gsw)
tukey_gsw <- HSD.test(anova_gsw, trt = "shortdata_D$genotype")
tukey_gsw$groups
#shortdata_D$gsw groups
#M3GM-     0.06408050      a
#M3GM--    0.04767559     ab
#MMY       0.04179032     ab
#M3GM++    0.03990685     ab
#sid       0.02828226     ab
#WT        0.02478372      b
#MYM       0.02467923      b
#M3GM+     0.02220491      b
sigphys$gsw_D <- c("a", "a", "ab", "a", "b", "ab", "ab", "ab")

# A - D
anova_A <- aov(shortdata_D$A ~shortdata_D$genotype)
summary(anova_A)
tukey_A <- HSD.test(anova_A, trt = "shortdata_D$genotype")
tukey_A$groups
# no significant difference
sigphys$A_D <- c("a")


## Plots for Fig. S6
plotS6_A <- ggplot(shortdata_LL, aes(x=factor(genotype, levels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=gsw, colour=genotype))+
  geom_boxplot(outlier.shape = NA, colour = "black", show.legend = F)+
  geom_jitter(height = 0, width = 0.1, alpha = 0.7, show.legend = F)+
  geom_text(sigphys, mapping=aes(x=genotype, y=0.45, label=paste(gsw_LL)), colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0, 0.45))+
  scale_colour_manual(values = c("WT"="#98768e", "MYM"="#b08ba5", "M3GM++"="#c7a2b6", "M3GM+"="#dfbbc8", 
                                 "M3GM-"="#ffc680", "M3GM--"="#ffb178", "MMY"="#db8872", "sid"="#a56457"))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "*sid*"))+
  theme_grey()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(text = element_text(size = 8), 
        axis.title.y = element_markdown(),
        axis.text.x = element_markdown())+
  theme(panel.background = element_rect(colour = 'black'))+
  labs(x=NULL,y="Absolute *g*<sub>SW</sub> [mol m<sup>-2</sup> s<sup>-1</sup>]",title=NULL)

plotS6_B <- ggplot(shortdata_LL, aes(x=factor(genotype, levels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=A, colour=genotype))+
  geom_boxplot(outlier.shape = NA, colour = "black", show.legend = F)+
  geom_jitter(height = 0, width = 0.1, alpha = 0.7, show.legend = F)+
  geom_text(sigphys, mapping=aes(x=genotype, y=30, label=paste(A_LL)), colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0, 30))+
  scale_colour_manual(values = c("WT"="#98768e", "MYM"="#b08ba5", "M3GM++"="#c7a2b6", "M3GM+"="#dfbbc8", 
                                 "M3GM-"="#ffc680", "M3GM--"="#ffb178", "MMY"="#db8872", "sid"="#a56457"))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "*sid*"))+
  theme_grey()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(text = element_text(size = 8), 
        axis.title.y = element_markdown(),
        axis.text.x = element_markdown())+
  theme(panel.background = element_rect(colour = 'black'))+
  labs(x=NULL,y="*A* [µmol m<sup>-2</sup> s<sup>-1</sup>]",title=NULL)

plotS6_C <- ggplot(shortdata_HL, aes(x=factor(genotype, levels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=gsw, colour=genotype))+
  geom_boxplot(outlier.shape = NA, colour = "black", show.legend = F)+
  geom_jitter(height = 0, width = 0.1, alpha = 0.7, show.legend = F)+
  geom_text(sigphys, mapping=aes(x=genotype, y=0.45, label=paste(gsw_HL)), colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0, 0.45))+
  scale_colour_manual(values = c("WT"="#98768e", "MYM"="#b08ba5", "M3GM++"="#c7a2b6", "M3GM+"="#dfbbc8", 
                                 "M3GM-"="#ffc680", "M3GM--"="#ffb178", "MMY"="#db8872", "sid"="#a56457"))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "*sid*"))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(text = element_text(size = 8), 
        axis.title.y = element_markdown(),
        axis.text.x = element_markdown())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-16)),
        legend.background = element_rect(fill=NA))+
  labs(x=NULL,y="Absolute *g*<sub>SW</sub> [mol m<sup>-2</sup> s<sup>-1</sup>]",title=NULL)

plotS6_D <- ggplot(shortdata_HL, aes(x=factor(genotype, levels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=A, colour=genotype))+
  geom_boxplot(outlier.shape = NA, colour = "black", show.legend = F)+
  geom_jitter(height = 0, width = 0.1, alpha = 0.7, show.legend = F)+
  geom_text(sigphys, mapping=aes(x=genotype, y=30, label=paste(A_HL)), colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0, 30))+
  scale_colour_manual(values = c("WT"="#98768e", "MYM"="#b08ba5", "M3GM++"="#c7a2b6", "M3GM+"="#dfbbc8", 
                                 "M3GM-"="#ffc680", "M3GM--"="#ffb178", "MMY"="#db8872", "sid"="#a56457"))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "*sid*"))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(text = element_text(size = 8), 
        axis.title.y = element_markdown(),
        axis.text.x = element_markdown())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-16)),
        legend.background = element_rect(fill=NA))+
  labs(x=NULL,y="*A* [µmol m<sup>-2</sup> s<sup>-1</sup>]",title=NULL)

plotS6_E <- ggplot(shortdata_D, aes(x=factor(genotype, levels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=gsw, colour=genotype))+
  geom_boxplot(outlier.shape = NA, colour = "black", show.legend = F)+
  geom_jitter(height = 0, width = 0.1, alpha = 0.7, show.legend = F)+
  geom_text(sigphys, mapping=aes(x=genotype, y=0.45, label=paste(gsw_D)), colour = "white", size = 8/.pt)+
  scale_y_continuous(limits = c(0, 0.45))+
  scale_colour_manual(values = c("WT"="#98768e", "MYM"="#b08ba5", "M3GM++"="#c7a2b6", "M3GM+"="#dfbbc8", 
                                 "M3GM-"="#ffc680", "M3GM--"="#ffb178", "MMY"="#db8872", "sid"="#a56457"))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "*sid*"))+
  theme_dark()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(text = element_text(size = 8), 
        axis.title.y = element_markdown(),
        axis.text.x = element_markdown())+
  theme(panel.background = element_rect(colour = 'black'))+
  labs(x=NULL,y="Absolute *g*<sub>SW</sub> [mol m<sup>-2</sup> s<sup>-1</sup>]",title=NULL)

## Correct text size and assemble Fig. S6
plotS6_A <- plotS6_A + theme(axis.text = element_text(size = 8),
                           legend.text = element_text(size = 8))
plotS6_B <- plotS6_B + theme(axis.text = element_text(size = 8),
                           legend.text = element_text(size = 8))
plotS6_C <- plotS6_C + theme(axis.text = element_text(size = 8),
                           legend.text = element_text(size = 8))
plotS6_D <- plotS6_D + theme(axis.text = element_text(size = 8),
                           legend.text = element_text(size = 8))
plotS6_E <- plotS6_E + theme(axis.text = element_text(size = 8),
                           legend.text = element_text(size = 8))

ggarrange(ggarrange(plotS6_A, plotS6_C, plotS6_E, nrow = 3, labels = c("A", "C", "E"), font.label = list(size = 12)),
          ggarrange(plotS6_B, plotS6_D, nrow = 3, labels = c("B", "D"), font.label = list(size = 12)), nrow = 1)


## Plots for Fig. S7
plotS7_A <- ggplot(shortdata_LL, aes(x=factor(genotype, levels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=iWUE, colour=genotype))+
  geom_boxplot(outlier.shape = NA, colour = "black", show.legend = F)+
  geom_jitter(height = 0, width = 0.1, alpha = 0.8, show.legend = F)+
  geom_text(sigphys, mapping=aes(x=genotype, y=150, label=paste(iWUE_LL)), colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0, 150))+
  scale_colour_manual(values = c("WT"="#98768e", "MYM"="#b08ba5", "M3GM++"="#c7a2b6", "M3GM+"="#dfbbc8", 
                                 "M3GM-"="#ffc680", "M3GM--"="#ffb178", "MMY"="#db8872", "sid"="#a56457"))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "*sid*"))+
  theme_grey()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(text = element_text(size = 8), 
        axis.title.y = element_markdown(),
        axis.text.x = element_markdown())+
  theme(panel.background = element_rect(colour = 'black'))+
  labs(x=NULL,y="iWUE [µmol(C) mol(H<sub>2</sub>O)<sup>-1</sup>]",title=NULL)

plotS7_B <- ggplot(shortdata_HL, aes(x=factor(genotype, levels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=iWUE, colour=genotype))+
  geom_boxplot(outlier.shape = NA, colour = "black", show.legend = F)+
  geom_jitter(height = 0, width = 0.1, alpha = 0.8, show.legend = F)+
  geom_text(sigphys, mapping=aes(x=genotype, y=150, label=paste(iWUE_HL)), colour = "black", size = 8/.pt)+
  scale_y_continuous(limits = c(0, 150))+
  scale_colour_manual(values = c("WT"="#98768e", "MYM"="#b08ba5", "M3GM++"="#c7a2b6", "M3GM+"="#dfbbc8", 
                                 "M3GM-"="#ffc680", "M3GM--"="#ffb178", "MMY"="#db8872", "sid"="#a56457"))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "*sid*"))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank())+
  theme(text = element_text(size = 8), 
        axis.title.y = element_markdown(),
        axis.text.x = element_markdown())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-16)),
        legend.background = element_rect(fill=NA))+
  labs(x=NULL,y="iWUE [µmol(C) mol(H<sub>2</sub>O)<sup>-1</sup>]",title=NULL)

## Correct text size and assemble Fig. S7
plotS7_A <- plotS7_A + theme(axis.text = element_text(size = 8),
                             legend.text = element_text(size = 8))
plotS7_B <- plotS7_B + theme(axis.text = element_text(size = 8),
                             legend.text = element_text(size = 8))

ggarrange(plotS7_A, plotS7_B, nrow = 2, labels = c("A", "B"), font.label = list(size = 12))


### Fig. S9 (physiology corrected by stomatal density, i.e. physiology per stoma)
## Densities (based on the mean density of Fig. S8C)
# line   density   
# M3GM+  109 
# M3GM++ 109 
# M3GM-   94
# M3GM--  96 
# MMY     93 
# MYM    110
# WT     105
# sid     89

## Plot S9A: Stomatal conductance (gsw) of WT, MYM, MMY and sid
plotS9_A <- licorplots(c("wt", "mym", "mmy", "sid"), timestamps = c(20,40,60), timeframe = c(11:75), 
                       legend_labels = c("WT", "MYM", "MMY", "*sid*"), legend_title = NULL,
                       remove_outliers = "yes", stomden = c(105, 110, 93, 89),
                       type = "gsw", colours = c("#98768e", "#b08ba5", "#db8872", "#a56457"))

## Plot S9B: Stomatal conductance (gsw) of WT, M3GM (++, +, -, --) and sid
plotS9_B <- licorplots(c("wt", "LB37", "LB14", "m3gm-LB3", "LB2", "sid"), timestamps = c(20,40,60), timeframe = c(11:75), 
                       legend_labels = c("WT", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "*sid*"), legend_title = NULL,
                       remove_outliers = "yes", stomden = c(105, 109, 109, 94, 96, 89),
                       type = "gsw", colours = c("#98768e", "#c7a2b6", "#dfbbc8", "#ffc680", "#ffb178", "#a56457"))

## Plot S9C: Carbon assimilation (A) of WT, MYM, MMY and sid
plotS9_C <- licorplots(c("wt", "mym", "mmy", "sid"), timestamps = c(20,40,60), timeframe = c(11:59), 
                       legend_labels = c("WT", "MYM", "MMY", "*sid*"), legend_title = NULL,
                       remove_outliers = "yes", stomden = c(105, 110, 93, 89),
                       type = "A", colours = c("#98768e", "#b08ba5", "#db8872", "#a56457"))

## Plot S9D: Carbon assimilation (A) of WT, M3GM (++, +, -, --) and sid
plotS9_D <- licorplots(c("wt", "LB37", "LB14", "m3gm-LB3", "LB2", "sid"), timestamps = c(20,40, 60), timeframe = c(11:59), 
                       legend_labels = c("WT", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "*sid*"), legend_title = NULL,
                       remove_outliers = "yes", stomden = c(105, 109, 109, 94, 96, 89),
                       type = "A", colours = c("#98768e", "#c7a2b6", "#dfbbc8", "#ffc680", "#ffb178", "#a56457"))

## Plot S9E: Intrinsic water-use efficiency (iWUE) of WT, MYM, MMY and sid
plotS9_E <- licorplots(c("wt", "mym", "mmy", "sid"), timestamps = c(20,40, 60), timeframe = c(11:59), 
                       legend_labels = c("WT", "MYM", "MMY", "*sid*"), legend_title = NULL,
                       remove_outliers = "yes", stomden = c(105, 110, 93, 89),
                       type = "WUE", colours = c("#98768e", "#b08ba5", "#db8872", "#a56457"))

## Plot S9F: Intrinsic water-use efficiency (iWUE) of WT, M3GM (++, +, -, --) and sid
plotS9_F <- licorplots(c("wt", "LB37", "LB14", "m3gm-LB3", "LB2", "sid"), timestamps = c(20,40,60), timeframe = c(11:59), 
                       legend_labels = c("WT", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "*sid*"), legend_title = NULL,
                       remove_outliers = "yes", stomden = c(105, 109, 109, 94, 96, 89),
                       type = "WUE", colours = c("#98768e", "#c7a2b6", "#dfbbc8", "#ffc680", "#ffb178", "#a56457"))

## correct text sizes and assemble Fig. S9
plotS9_A_final <- plotS9_A + theme(text = element_text(size = 8),
                                 axis.text = element_text(size = 8))
plotS9_B_final <- plotS9_B + theme(text = element_text(size = 8),
                                 axis.text = element_text(size = 8))
plotS9_C_final <- plotS9_C + theme(text = element_text(size = 8),
                                 axis.text = element_text(size = 8))
plotS9_D_final <- plotS9_D + theme(text = element_text(size = 8),
                                 axis.text = element_text(size = 8))
plotS9_E_final <- plotS9_E + theme(text = element_text(size = 8),
                                 axis.text = element_text(size = 8))
plotS9_F_final <- plotS9_F + theme(text = element_text(size = 8),
                                 axis.text = element_text(size = 8))

ggarrange(ggarrange(plotS9_A_final, plotS9_C_final, plotS9_E_final, nrow = 3, 
                    common.legend = T, legend = "none", labels = c("A", "C", "E"), font.label = list(size = 12), align = "v"), 
          ggarrange(plotS9_B_final, plotS9_D_final, plotS9_F_final, nrow = 3, 
                    common.legend = T, legend = "none", labels = c("B", "D", "F"), font.label = list(size = 12), align = "v"), 
          ncol=2)





#### Fig. S8 -------------------------------------------------------------------------------------
### load data
pheno <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_FigS8")


### set colour palette
colourblindpalette <- c("#295384", "#5a97c1", "#d8d97a", "#95c36e")


### Fig. S8A
## calculate means for all phenotypes for each line
means1<- pheno %>% group_by(line) %>% 
  summarize(wt=mean(`wt-like`),
            one_SC=mean(`one_SC`),
            no_SC=mean(`no_SC`),
            aborted=mean(aborted),
            pictures=n())


## reorder columns to have one column 'type' (contains name of phenotype) and another with the percentages ('percent')
means2<- means1 %>% select(-pictures) %>%
  pivot_longer(., cols = c(wt, `one_SC`, `no_SC`, aborted), names_to = "phenotype", values_to = "percent")


## plot phenotypes displayed as stacked plot (types in one bar)
plotS8_A <- means2 %>% mutate(phenotype = fct_relevel(phenotype, "aborted", "no_SC", "one_SC", "wt")) %>%
  ggplot(aes(x=factor(line, levels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=percent, fill=phenotype))+
  geom_col(position="fill", show.legend = T)+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 25, 50, 75, 100))+
  scale_fill_manual(labels=c("wt"="WT", "one_SC"="1 SC", "no_SC"="0 SC", "aborted"="aborted"),
                    values=colourblindpalette)+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--"="M3GM- -", "MMY", "sid"="*sid*"))+
  guides(fill = guide_legend(override.aes = list(shape = 22), 
                             title = "Phenotype",
                             reverse = TRUE))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_markdown())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-16)),
        legend.background = element_rect(fill=NA))+
  labs(x=NULL,y="Percentage [%]",title=NULL)


### Fig. S8B
## calculate percentages for each stomatal phenotype
pheno$percentwt <- pheno$`wt-like`/pheno$stomata
pheno$percent1sc <- pheno$one_SC/pheno$stomata
pheno$percent0sc <- pheno$no_SC/pheno$stomata
pheno$percentaborted <- pheno$aborted/pheno$stomata

## reorder data
phenopercentages <- data.frame(line=pheno$line, identifier=pheno$identifier, individual=pheno$individual, type="WT", percent=pheno$percentwt)
phenopercentages <- rbind(phenopercentages, data.frame(line=pheno$line, identifier=pheno$identifier, individual=pheno$individual, type="one_SC", percent=pheno$percent1sc))
phenopercentages <- rbind(phenopercentages, data.frame(line=pheno$line, identifier=pheno$identifier, individual=pheno$individual, type="no_SC", percent=pheno$percent0sc))
phenopercentages <- rbind(phenopercentages, data.frame(line=pheno$line, identifier=pheno$identifier, individual=pheno$individual, type="aborted", percent=pheno$percentaborted))

percentmeans <- phenopercentages %>% group_by(line, identifier, individual, type) %>% summarise(mean_percent=mean(percent))

## plot phenotypes as box- and jitterplots
plotS8_B <- ggplot(percentmeans, aes(x=factor(line, levels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid")), y=mean_percent, colour=factor(type, levels = c("WT", "one_SC", "no_SC", "aborted"))))+
  geom_boxplot(position="dodge", show.legend = F, outlier.shape = NA)+
  geom_point(position = position_dodge(width = 0.75), show.legend = T, alpha = 0.8)+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 25, 50, 75, 100))+
  scale_colour_manual(labels=c("wt"="WT", "one_SC"="1 SC", "no_SC"="0 SC", 
                               "aborted"="aborted"),
                      values=rev(colourblindpalette))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--"="M3GM- -", "MMY", "sid"="*sid*"))+
  guides(colour = guide_legend(title = "Phenotype"))+
  theme_bw()+
  theme(strip.background = element_rect(colour=NA,fill=NA))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_markdown())+
  theme(legend.position = "bottom",
        legend.justification="left",
        legend.box.margin = margin(c(-16)),
        legend.background = element_rect(fill=NA))+
  labs(x=NULL, y="Percentage [%]",title=NULL)


### Fig. S8C
### test for significant differences in comparison to the wild type that was used in the same round
stomden1 <- pheno %>% group_by(line, identifier, individual) %>% summarise(mean_dens= mean(stomatal_density_functional))

## Statistical analysis (significant differences between a genotype and the respective WT individuals in the same experiment)
# LB1 (WT x4) with LB2 (M3GM-- x5) and LB3 (M3GM- x4) - no significant difference
LB1_2_3 <- stomden1 %>% filter(identifier %in% c("LB1", "LB2", "LB3"))
a1 <- aov(LB1_2_3$mean_dens ~LB1_2_3$line)
summary(a1)
t1 <- HSD.test(a1, trt="LB1_2_3$line")
t1$groups

# LB8 (WT x5) with LB9 (MYM x5) and LB10 (MMY x5) - no significant difference
LB8_9_10 <- stomden1 %>% filter(identifier %in% c("LB8", "LB9", "LB10"))
a2 <- aov(LB8_9_10$mean_dens ~LB8_9_10$line)
summary(a2)
t2 <- HSD.test(a2, trt="LB8_9_10$line")
t2$groups

# LB11 (WT x6) with LB12 (sid x5) - significant difference
LB11_12 <- stomden1 %>% filter(identifier %in% c("LB11", "LB12"))
a3 <- aov(LB11_12$mean_dens ~LB11_12$line)
summary(a3)
t3 <- HSD.test(a3, trt="LB11_12$line")
t3$groups
#LB11_12$mean_dens groups
#WT          114.61187      a
#sid          89.31507      b

# LB13 (WT) with LB14 (M3GM+) - significant difference
LB13_14 <- stomden1 %>% filter(identifier %in% c("LB13", "LB14"))
a4 <- aov(LB13_14$mean_dens ~LB13_14$line)
summary(a4)
t4 <- HSD.test(a4, trt="LB13_14$line")
t4$groups
#LB13_14$mean_dens groups
#M3GM+         108.90411      a
#WT             90.63927      b

# LB36 (WT x6) with LB37 (M3GM++ x6) - no significant difference
LB36_37 <- stomden1 %>% filter(identifier %in% c("LB36", "LB37"))
a5 <- aov(LB36_37$mean_dens ~LB36_37$line)
summary(a5)
t5 <- HSD.test(a5, trt="LB36_37$line")
t5$groups

## create data frame with significant differences marked (in comparison to WT)
sigden <- data.frame(line=c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid"),
                     significance=c("", "", "", "*", "", "", "", "*"))

## reorder data
stomden1$line <- ordered(stomden1$line , levels=c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM--", "MMY", "sid"))

## plot
plotS8_C <- ggplot(stomden1, mapping=aes(x=line, y=mean_dens, colour = line))+
  geom_boxplot(colour = "black")+
  geom_text(sigden, mapping=aes(y=150, label=paste(significance)), show.legend = F, colour = "black")+
  geom_jitter(width=0.1, height = 0, alpha = 0.5, show.legend = F)+
  scale_y_continuous(limits = c(0, 150))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "M3GM+", "M3GM-", "M3GM- -", "MMY", "*sid*"))+
  scale_colour_manual(values = c("#98768e", "#b08ba5", "#c7a2b6", "#dfbbc8", "#ffc680", "#ffb178", "#db8872", "#a56457"))+
  theme_classic()+
  guides(colour = guide_legend(title = NULL))+
  theme(axis.text.x = element_markdown())+
  labs(x=NULL, y=expression(paste("Stomatal density [stomata mm"^-2, "]")))



### correct text size and assemble Fig. S8
plotS8_A_final <- plotS8_A + theme(text = element_text(size = 8),
                            axis.text = element_text(size = 8),
                            legend.text = element_text(size = 8))
plotS8_B_final <- plotS8_B + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              legend.text = element_text(size = 8))
plotS8_C_final <- plotS8_C + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              legend.text = element_text(size = 8))

ggarrange(ggarrange(plotS8_A_final, plotS8_B_final, nrow=2, common.legend = T, legend = "right", labels = c("A", "B"), font.label = list(size = 12)),
          plotS8_C_final, nrow = 2, labels = c("", "C"), heights = c(1.5,1), font.label = list(size = 12))





#### Fig. S10 -------------------------------------------------------------------------------------
### load data
stomin <- read_excel("MUTE2024_Data.xlsx", sheet = "Data_FigS10")
stomin <- subset(stomin, !is.na(stomin$stomatal_index))


### calculate means per individual and line and days after germination
stominmean <- stomin %>% group_by(line, individual, dap) %>% summarise(si = mean(stomatal_index), 
                                                                       sli = mean(stomatal_lineage_index),
                                                                       sd = mean(stomatal_density),
                                                                       sld = mean(stomatal_lineage_density),
                                                                       sum_stomata=sum(stomatal_lineage))


### statistical analysis
## Stomatal index
# ANOVA
indexinter <- interaction(stominmean$line, stominmean$dap)
anovasi <- aov(stominmean$si ~indexinter)
summary(anovasi)
# Tukey's test
tukeysi <- HSD.test(anovasi, trt="indexinter")
tukeysi$groups

# create data frame with significance levels
sigsi <- data.frame(line=c("WT", "MYM", "M3GM++", "sid", "WT", "MYM", "M3GM++", "sid"),
                    dap=c(12, 12, 12, 12, 28, 28, 28, 28),
                    significance=c("a", "ab", "a", "c", "a", "a", "ab", "bc"))

## Stomatal lineage index
# ANOVA
indexinter <- interaction(stominmean$line, stominmean$dap)
anovasli <- aov(stominmean$sli ~indexinter)
summary(anovasli)
# Tukey's test
tukeysli <- HSD.test(anovasli, trt="indexinter")
tukeysli$groups

# create data frame with significance levels
sigsli <- data.frame(line=c("WT", "MYM", "M3GM++", "sid", "WT", "MYM", "M3GM++", "sid"),
                     dap=c(12, 12, 12, 12, 28, 28, 28, 28),
                     significance=c("a"))

## Stomatal density
# ANOVA
densityinter <- interaction(stominmean$line, stominmean$dap)
anovasd <- aov(stominmean$sd ~densityinter)
summary(anovasd)
# Tukey's test
tukeysd <- HSD.test(anovasd, trt="densityinter")
tukeysd$groups

# create data frame with significance levels
sigsd <- data.frame(line=c("WT", "MYM", "M3GM++", "sid", "WT", "MYM", "M3GM++", "sid"),
                    dap=c(12, 12, 12, 12, 28, 28, 28, 28),
                    significance=c("a", "a", "a", "c", "b", "b", "b", "a"))

## Stomatal lineage density
# ANOVA
densityinter <- interaction(stominmean$line, stominmean$dap)
anovasld <- aov(stominmean$sld ~densityinter)
summary(anovasld)
# Tukey's test
tukeysld <- HSD.test(anovasld, trt="densityinter")
tukeysld$groups

# create data frame with significance levels
sigsld <- data.frame(line=c("WT", "MYM", "M3GM++", "sid", "WT", "MYM", "M3GM++", "sid"),
                     dap=c(12, 12, 12, 12, 28, 28, 28, 28),
                     significance=c("a", "a", "a", "a", "b", "b", "b", "b"))


### reorder data frame
stominmean$line <- ordered(stominmean$line , levels=c("WT", "MYM", "M3GM++", "sid"))


### plots
facets <- c("12" = "12 days after germination", "28" = "28 days after germination")

## stomatal index
plotS10_A <- ggplot(stominmean, mapping=aes(x=line, y=si))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(sigsi, mapping=aes(y=12.5, label=paste(significance)), show.legend = F, size = 8/.pt)+
  geom_jitter(width=0.1, height = 0, alpha = 0.5)+
  facet_grid(~dap, label = as_labeller(facets))+
  scale_y_continuous(limits = c(0, 12.5))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "*sid*"))+
  theme_classic()+
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
        text = element_text(size = 8))+
  labs(x=NULL, y="Stomatal index [%]")

## stomatal lineage index
plotS10_B <- ggplot(stominmean, mapping=aes(x=line, y=sli))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(sigsli, mapping=aes(y=12.5, label=paste(significance)), show.legend = F, size = 8/.pt)+
  geom_jitter(width=0.1, height = 0, alpha = 0.5)+
  facet_wrap(~dap, label = as_labeller(facets))+
  scale_y_continuous(limits = c(0, 12.5))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "*sid*"))+
  theme_classic()+
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
        text = element_text(size = 8))+
  labs(x=NULL, y="Stomatal lineage index [%]")

## stomatal density
plotS10_C <- ggplot(stominmean, mapping=aes(x=line, y=sd))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(sigsd, mapping=aes(y=125, label=paste(significance)), show.legend = F, size = 8/.pt)+
  geom_jitter(width=0.1, height = 0, alpha = 0.5)+
  facet_wrap(~dap, label = as_labeller(facets))+
  scale_y_continuous(limits = c(0, 125))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "*sid*"))+
  theme_classic()+
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.title.y = element_markdown(),
        text = element_text(size = 8))+
  labs(x=NULL, y="Stomatal density [mm<sup>-2</sup>]")

## stomatal lineage density
plotS10_D <- ggplot(stominmean, mapping=aes(x=line, y=sld))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(sigsld, mapping=aes(y=125, label=paste(significance)), show.legend = F, size = 8/.pt)+
  geom_jitter(width=0.1, height = 0, alpha = 0.5)+
  facet_wrap(~dap, label = as_labeller(facets))+
  scale_y_continuous(limits = c(0, 125))+
  scale_x_discrete(labels = c("WT", "MYM", "M3GM++", "*sid*"))+
  theme_classic()+
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1), 
        axis.title.y = element_markdown(),
        text = element_text(size = 8))+
  labs(x=NULL, y="Stomatal lineage density [mm<sup>-2</sup>]")


### correct text size and assemble Fig. S10
plotS10_A_final <- plotS10_A + theme(text = element_text(size = 8),
                              axis.text = element_text(size = 8),
                              strip.text = element_text(size = 8),
                              legend.text = element_text(size = 8))

plotS10_B_final <- plotS10_B + theme(text = element_text(size = 8),
                                     axis.text = element_text(size = 8),
                                     strip.text = element_text(size = 8),
                                     legend.text = element_text(size = 8))

plotS10_C_final <- plotS10_C + theme(text = element_text(size = 8),
                                     axis.text = element_text(size = 8),
                                     strip.text = element_text(size = 8),
                                     legend.text = element_text(size = 8))

plotS10_D_final <- plotS10_D + theme(text = element_text(size = 8),
                                     axis.text = element_text(size = 8),
                                     strip.text = element_text(size = 8),
                                     legend.text = element_text(size = 8))

ggarrange(plotS10_A_final, plotS10_B_final, plotS10_C_final, plotS10_D_final, labels = c("A", "B", "C", "D"), font.label = list(size = 12))
