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

## Does evolutionary history impact community productivity?

#load dataset
prod <- read.csv("Com_prod_final.csv", header = T)

summary(prod)
##dataset key
#ID = community identity
#Coev = evolutionary history of community members. M = monoculture-evolved; N = polyculture-evolved from allopatric communities; Y = polyculture-evolved from sympatric communities

##    linear model testing community productivity against evolutionary history

prmod <- lm(prod ~ Coev, data = prod)
prodn <- lm(prod ~ 1, data = prod)
anova(prmod, prodn, test = "F") #non-significant effect of evolutionary history

#Figure 5a

lab <- c("Monoculture", expression(atop("Polyculture",
                                        "(Allopatric)")),
         expression(atop("Polyculture",
                         "(Sympatric)")))

cols_3 <- c("#4682b4", "#801061", "#cd8500") 

ggplot(prod, aes(x = Coev, y = prod, col = Coev), na.rm = T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.1), aes(alpha = 0.5), size = 3) +
  scale_x_discrete(labels = lab) +
  theme_bw() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14, colour = "black"), axis.text.x =element_text(size = 14, colour = "black"), legend.position = "none", title = element_text(size = 14)) +
  xlab("Evolutionary history of community members") +
  ylab("Productivity") +
  scale_colour_manual(values = cols_3) +
  labs(title = element_text("(a) Community productivity")) +
  scale_alpha(guide = 'none')


###   Do communities differ in structure (species proportions)?

#load dataset
Props <- read.csv("community_proportions_060420.csv", header = T)

summary(Props)
##dataset key:
#sp = species.  A = Achromobacter sp. O = Ochrobactrum sp. S = Stenotrophomonas sp. P = Pseudomonas sp. V = Variovorax sp. 
#comm = community evolution line and type. ".Ran" = community assembled from clones from random allopatric polyculture lines. ".Com" = community assembled from clones from the same evolution line (sympatric). ".Mono" = monoculture-based community. 
#count = number of colonies of focal species
#total = total number of colonies

#calculate proportion of each species within their community
Props <- mutate(Props, 
                prop = count/total)

#generalised linear mixed effects model testing proportion against species and community assembly
y <- cbind(Props$count, Props$total-Props$count)

com_props <- glm(y ~ sp * treat, data = Props, family = binomial)
summary(com_props) #heavily overdispersed

com_props <- glm(y ~ sp * treat, data = Props, family = quasibinomial) #switch to a quasibinomial error structure to deal with overdispersion
com_props2 <- glm(y ~ sp + treat, data = Props, family = quasibinomial)
anova(com_props2, com_props, test = "Chi") #sig interaction

emmeans::emmeans(com_props, pairwise ~ treat|sp, type = "response") #multiple pairwise comparisons comparing proportions of each evolved line within species

#Figure 5b

Props2 <- Props

Props2$Community <- as.character(Props2$treat)
Props2$Community[Props2$Community == "C"] <- "Poly. sym."
Props2$Community[Props2$Community == "M"] <- "Mono."
Props2$Community[Props2$Community == "R"] <- "Poly. allo."
Props2$sp <- as.character(Props2$sp)
Props2$sp <- substr(Props2$sp, 1, 1)

ggplot(Props2, aes(x = sp, y = prop, col = Community)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1), aes(alpha = 0.5)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14, colour = "black"), axis.text = element_text(size = 12, colour = "black"), axis.text.x =element_text(size = 12, colour = "black"), strip.background = element_rect(fill="white"), legend.position = "bottom", title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  ylab("Proportion") +
  xlab("Species") +
  scale_colour_manual(values = cols_3) +
  labs(title = element_text("(b) Community structure")) +
  scale_alpha(guide = 'none') +
  scale_x_discrete(labels = sp.labs)

#Table of proportion comparisons

Props2 <- Props
Props2$Community <- as.character(Props2$treat)
Props2$Community[Props2$Community == "C"] <- "Poly. sym."
Props2$Community[Props2$Community == "M"] <- "Mono."
Props2$Community[Props2$Community == "R"] <- "Poly. allo."
Props2$sp <- as.character(Props2$sp)
Props2$sp <- substr(Props2$sp, 1, 1)

com_props <- glm(y ~ sp * Community, data = Props2, family = quasibinomial)
prop_reads <- emmeans::emmeans(com_props, pairwise ~ Community|sp, type = "response")
prop_tab <- data.frame(prop_reads)
prop_tab <- prop_tab[c(1:15), -c(1:7, 12)]

prop_tab <- mutate(prop_tab,
                   contrasts.odds.ratio = round(contrasts.odds.ratio, digits = 3),
                   contrasts.SE = round(contrasts.SE, digits = 3),
                   contrasts.z.ratio = round(contrasts.z.ratio, digits = 3),
                   contrasts.p.value = round(contrasts.p.value, digits = 3))

prop_tab$contrasts.odds.ratio <- as.character(prop_tab$contrasts.odds.ratio)
prop_tab$contrasts.SE <- as.character(prop_tab$contrasts.SE)
prop_tab$contrasts.z.ratio <- as.character(prop_tab$contrasts.z.ratio)
prop_tab$contrasts.p.value <- as.character(prop_tab$contrasts.p.value)

prop_flex <- flextable(prop_tab) %>%
  set_header_labels(contrasts.odds.ratio = 'Odds-ratio', contrasts.sp = 'Sp.', contrasts.SE = 'SE', contrasts.z.ratio = 'z-ratio', contrasts.p.value = 'p-value', contrasts.contrast = "Contrast") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit() %>%
  bold(i = c(5, 9, 14, 15)) %>%
  hline(i = c(3, 6, 9, 12, 15), border = officer::fp_border(color="black")) %>%
  align(j = c(1:6), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

prop_flex