#load packages

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(ggridges)
library(stringr)
library(lme4)
library(emmeans)
library(nlme)
library(pbkrtest)
library(ggsignif)
library(flextable)
library(readxl)
library(rmarkdown)
library(webshot)

### Do communities remain stable over evolutionary time?

#load dataset

inv.dat <- read.csv("Community_stability_060420.csv", header = T)

summary(inv.dat)
#dataset key
#Res.com.ID = for polyculture-based communities (Res.type = "Com"), this refers to the evolution line from which the clones originated. For monoculture-based communities (Res.type = "Mono"), this refers to which block clones were assigned to form a community
#Inv.sp = invading species. A = Achromobacter sp. O = Ochrobactrum sp. P = Pseudomonas sp. S = Stenotrophomonas sp. V = Variovorax sp.
#coev = evolutionary history of invading species. C = polyculture-evolved. M = monoculture-evolved.
#Res.type = evolutionary history of resident species which matches the evolutionary history of the invading species
#rel.fitness = relative invader fitness

## Are relative invader fitness of all species greater than 1?

#grab data columns needed
inv.dat2 <- inv.dat[c(1:120), c(2, 3, 5)]

#run multiple one-sampled t-tests testing each species' relative invader fitness against one within each community
d_mods <- inv.dat2 %>%
  nest(-c(Inv.sp, coev)) %>%
  mutate(., mod = map(data, ~ t.test(.x$rel.fitness, mu = 1)),
         tidied = map(mod, broom::tidy))

#adjust p-values for multiple testing using false-discovery rate (fdr)
d_mods_summary <- d_mods %>%
  select(-c(mod, data)) %>%
  unnest(tidied) %>%
  mutate(padjust = p.adjust(p.value, "fdr"))

#Table of one-sample t-test statistics

d_mods_summary <- d_mods_summary[order(d_mods_summary$Inv.sp),]

d_mods_summary$estimate <- round(d_mods_summary$estimate, digits = 3)
d_mods_summary$statistic <- round(d_mods_summary$statistic, digits = 3)
d_mods_summary$p.value <- round(d_mods_summary$p.value, digits = 3)
d_mods_summary$conf.low <- round(d_mods_summary$conf.low, digits = 3)
d_mods_summary$conf.high <- round(d_mods_summary$conf.high, digits = 3)
d_mods_summary$padjust <- round(d_mods_summary$padjust, digits = 3)

d_mods_summary <- d_mods_summary[c(1:10), c(1:8, 11)]
d_mods_summary$Inv.sp <- as.character(d_mods_summary$Inv.sp)
d_mods_summary$Inv.sp <- substr(d_mods_summary$Inv.sp, 1, 1)

d_mods_summary$coev <- as.character(d_mods_summary$coev)
d_mods_summary$coev[d_mods_summary$coev == "C"] <- "Polyculture-based"
d_mods_summary$coev[d_mods_summary$coev == "M"] <- "Monoculture-based"

d_mods_summary$estimate <- as.character(d_mods_summary$estimate)
d_mods_summary$statistic <- as.character(d_mods_summary$statistic)
d_mods_summary$p.value <- as.character(d_mods_summary$p.value)
d_mods_summary$conf.low <- as.character(d_mods_summary$conf.low)
d_mods_summary$conf.high <- as.character(d_mods_summary$conf.high)
d_mods_summary$padjust <- as.character(d_mods_summary$padjust)

d_mods_summary$p.value[d_mods_summary$p.value == "0"] <- "<0.001"
d_mods_summary$padjust[d_mods_summary$padjust == "0"] <- "<0.001"

inv_flex <- flextable(d_mods_summary) %>%
  set_header_labels(Inv.sp = "Sp.", coev = "Community type", estimate = 'Estimate', statistic = 't-value', p.value = 'p-value', parameter = 'df', conf.low = 'Lower CI', conf.high = 'Upper CI', padjust = "Adj. p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit() %>%
  hline(i = c(2, 4, 6, 8, 10), border = officer::fp_border(color="black")) %>%
  align(j = c(1:9), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header') 

inv_flex

## does relative invader fitness differ between community assemblies?

#linear model testing relative invader fitness against interactions of species and community assembly. 

invm1 <- lm(rel.fitness ~ Inv.sp * coev, data = inv.dat)
invm2 <- lm(rel.fitness ~ Inv.sp + coev, data = inv.dat)
anova(invm1, invm2, test = "F") #non-significant interaction

invm3 <- lm(rel.fitness ~ Inv.sp, data = inv.dat2)
anova(invm2, invm3, test = "F") #non-significant effect of community assembly

invm4 <- lm(rel.fitness ~ coev, data = inv.dat2)
anova(invm2, invm4, test = "F") #significant effect of species

invn <- lm(rel.fitness ~ 1, data = inv.dat2)
anova(invm3, invn, test = "F") #significant effect of species

#Figure 4

inv.dat$coev <- as.character(inv.dat$coev)
inv.dat$coev[inv.dat$coev == "C"] <- "Polyculture-based"
inv.dat$coev[inv.dat$coev == "M"] <- "Monoculture-based"
cols_2 <- c("#4682b4", "#cd8500")

inv.dat2$Inv.sp <- as.character(inv.dat2$Inv.sp)
inv.dat2$Inv.sp <- substr(inv.dat2$Inv.sp, 1, 1)

sp.labs <- c(expression(italic("Achromobacter sp.")), expression(italic("Ochrobactrum sp.")), expression(italic("Pseudomonas sp.")), expression(italic("Stenotrophomonas sp.")), expression(italic("Variovorax sp.")))

ggplot(inv.dat, aes(x = Inv.sp, y = rel.fitness, col = coev), na.rm = T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1), alpha = 0.5) +
  xlab("Species") +
  theme_bw() +
  ylab("Relative invader fitness") +
  theme(axis.title = element_text(size = 12, colour = "black"), axis.text = element_text(size = 12, colour = "black"), axis.text.x =element_text(size = 12, colour = "black"), strip.background = element_rect(fill="white"), legend.position = 'bottom', legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  labs(colour = "Community") +
  scale_alpha(guide = 'none') +
  scale_colour_manual(values = cols_2) +
  scale_x_discrete(labels = sp.labs)
