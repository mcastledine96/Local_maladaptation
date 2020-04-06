#load packages

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

### Do species become adapted to the presence of other community members?

#read in dataset
biotic_dat <- read.csv("Biotic_adaptation_060420.csv", header = T)

summary(biotic_dat)
##dataset key
#Treat (evolutionary history): M, monoculture-evolved; C, polyculture-evolved; An, Ancestral
#Res.type, evolutionary history of non-focal species are all coevolved in polyculture
#sp.fit = focal species growth rate
#Res.id = community line non-focal species originated from
#Sp = species. A = Achromobacter sp. O = Ochrobactrum sp. S = Stenotrophomonas sp. P = Pseudomonas sp. V = Variovorax sp. 

#   linear mixed effects model analysing growth rate against evolutionary history and species. Random effect of resident community evolutionary history

bm1 <- lmer(sp.fit ~ Treat * Sp + (1|Res.id), data = biotic_dat)
bm2 <- lmer(sp.fit ~ Treat + Sp + (1|Res.id), data = biotic_dat)
#model simplification
anova(bm1, bm2) #Interaction non-significant. Drop from the model

bm3 <- lmer(sp.fit ~ Sp + (1|Res.id), data = biotic_dat)
anova(bm2, bm3) #effect of treatment non-significant. Drop from the model

bm4 <- lmer(sp.fit ~ Treat + (1|Res.id), data = biotic_dat)
anova(bm2, bm4) #effect of species significant

bm5 <- lmer(sp.fit ~ 1 + (1|Res.id), data = biotic_dat)
anova(bm5, bm3) #Species remains as the only significant independent effect

#   Figure 2a

#make treatment labels more intuitive
biotic_dat$Treat <- as.character(biotic_dat$Treat)
biotic_dat$Treat[biotic_dat$Treat == "An"] <- "Anc."
biotic_dat$Treat[biotic_dat$Treat == "C"] <- "Poly."
biotic_dat$Treat[biotic_dat$Treat == "M"] <- "Mono."
names(biotic_dat)[3] <- "Population" #change column name of Treat to 'Population'

cols_1 <- c("black", "#4682b4", "#cd8500")

ggplot(biotic_dat, aes(x = Sp, y = sp.fit, col = Population)) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1), aes(alpha = 0.5)) +
  ylab(expression("Growth rate ("~italic(m)~")")) +
  xlab("Species") +
  theme(axis.title = element_text(size = 12, colour = "black"), axis.text = element_text(colour = "black", size = 11), axis.text.x =element_text(colour = "black"), strip.background = element_rect(fill="white")) +
  scale_colour_manual(values = cols_1) +
  scale_alpha(guide = 'none') +
  labs(title = element_text("(a) Biotic adaptation", size = 10))


### Do species populations differ in abiotic adaptation?

#read in dataset
abiotic_dat <- read.csv("abiotic_adaptation_060420.csv", header = T)

summary(abiotic_dat)
##dataset key
#Evo_line = evolution line from which species originated. For ancestral (AN number) clones, this  refers to treatment replicate
#Sp = focal species as previous
#Coev = evolutionary history of species' population. Anc = Ancestral; Com = Polyculture-evolved; Mono = Monoculture-evolved
#fit = growth rate in isolation

#Linear model - growth rate against interacting fixed effects of species and evolutionary history

am1 <- lm(fit ~ Coev * Sp, data = abiotic_dat)
am2 <- lm(fit ~ Coev + Sp, data = abiotic_dat)

anova(am1, am2, test = "F") #significant interaction
emmeans::emmeans(mod, pairwise ~ Coev | Sp) #pairwise comparisons of populations within species shows effect is driven by monoculture-evolved populations

#   Figure 2b

abiotic_dat$Coev <- as.character(abiotic_dat$Coev)
abiotic_dat$Coev[abiotic_dat$Coev == "Anc"] <- "Ancestral"
abiotic_dat$Coev[abiotic_dat$Coev == "Com"] <- "Polyculture-evolved"
abiotic_dat$Coev[abiotic_dat$Coev == "Mono"] <- "Monoculture-evolved"
names(abiotic_dat)[3] <- "Population"

sp.labs <- c(expression(italic("Achromobacter sp.")), expression(italic("Ochrobactrum sp.")), expression(italic("Pseudomonas sp.")), expression(italic("Stenotrophomonas sp.")), expression(italic("Variovorax sp.")))

ggplot(abiotic_dat, aes(x = Sp, y = fit, col = Population)) +
  theme_bw() +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1), aes(alpha = 0.5), size = 2.5) +
  ylab(expression("Growth rate ("~italic(m)~")")) +
  xlab("Species") +
  theme(axis.title = element_text(size = 14, colour = "black"), axis.text = element_text(colour = "black", size = 12), axis.text.x =element_text(colour = "black"), strip.background = element_rect(fill="white"), title = element_text(size = 14), legend.text = element_text(size = 12), legend.position = "bottom") +
  scale_colour_manual(values = cols_1) +
  scale_alpha(guide = 'none') +
  labs(title = element_text("(b) Abiotic adaptation")) +
  scale_x_discrete(labels = sp.labs)

#Table of pairwise comparisons - abiotic adaptation

abiotic_dat <- read.csv("abiotic_adaptation_060420.csv", header = T)
dat_grow <- abiotic_dat
dat_grow$Coev <- as.character(dat_grow$Coev)
dat_grow$Coev[dat_grow$Coev == "Com"] <- "Poly."
dat_grow$Coev[dat_grow$Coev == "Mono"] <- "Mono."
dat_grow$Coev[dat_grow$Coev == "Anc"] <- "Anc."
dat_grow$Sp <- as.character(dat_grow$Sp)
dat_grow$Sp <- substr(dat_grow$Sp, 1, 1)

mod_grow <- lm(fit ~ Coev * Sp, data = dat_grow)
grow_reads <- emmeans::emmeans(mod_grow, pairwise ~ Coev | Sp)
grow_tab <- data.frame(grow_reads)
grow_tab <- grow_tab[c(1:15), -c(1:7)]

grow_tab <- mutate(grow_tab,
                   contrasts.estimate = round(contrasts.estimate, digits = 3),
                   contrasts.SE = round(contrasts.SE, digits = 3),
                   contrasts.df = round(contrasts.df, digits = 0),
                   contrasts.t.ratio = round(contrasts.t.ratio, digits = 3),
                   contrasts.p.value = round(contrasts.p.value, digits = 3))

grow_tab$contrasts.p.value[grow_tab$contrasts.p.value == 0] <- '<0.001'
grow_tab$contrasts.estimate <- as.character(grow_tab$contrasts.estimate)
grow_tab$contrasts.SE <- as.character(grow_tab$contrasts.SE)
grow_tab$contrasts.df <- as.character(grow_tab$contrasts.df)
grow_tab$contrasts.p.value <- as.character(grow_tab$contrasts.p.value)

grow_flex <- flextable(grow_tab) %>%
  set_header_labels(contrasts.contrast = 'Contrast', contrasts.Sp = 'Sp.', contrasts.estimate = 'Estimate', contrasts.SE = 'SE', contrasts.df = 'df', contrasts.t.ratio = 't-ratio', contrasts.p.value = 'p-value') %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit() %>%
  bold(i = c(1, 13, 15)) %>%
  hline(i = c(3, 6, 9, 12, 15), border = officer::fp_border(color="black")) %>%
  align(j = c(2:7), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

grow_flex
