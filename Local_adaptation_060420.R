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

### Do species become locally adapted to their sympatric community?

#load in dataset
LA <- read.csv("local_adaptation_060420.csv", header = T)

summary(LA)
##dataset key:
#comm = evolution line from which clone originated
#clone = unique variable assigned to each species' clone
#Sp = species identity. A = Achromobacter sp. O = Ochrobactrum sp. P = Pseudomonas sp. S = Stenotrophomonas sp. V = Variovorax sp.
#Treat = whether growth rates was measured in the clones' sympatric (symp) or allopatric (allo) community
#fitness = growth rate

#linear mixed effects model testing growth rate against community type (symptric/allopatric). Random effect of clone to account for repeated measures design

LA_mod <- lmer(fitness ~ Treat * Sp + (1|clone), data = LA)
LA_mod2 <- lmer(fitness ~ Treat + Sp + (1|clone), data = LA)
anova(LA_mod, LA_mod2) #no interaction

LA_mod3 <- lmer(fitness ~ Sp + (1|clone), data = LA)
anova(LA_mod2, LA_mod3)

LA_mod4 <- lmer(fitness ~ Treat + (1|clone), data = LA)
anova(LA_mod2, LA_mod4) #both species and treatment (allopatric/sympatric) are significant

#multiple pairwise comparions contrasting growth rates with species across community types
emmeans::emmeans(LA_mod2, pairwise ~ Treat|Sp)

#see how many clones increased their growth rate in an allopatric community

#subset data
LA_sym <- filter(LA, Treat == "symp")
LA_allo <- filter(LA, Treat == "allo")
LA_sym$allo_fit <- LA_allo$fitness
#calculate difference in growth rate
LA_sym$diff <- LA_sym$allo_fit - LA_sym$fitness #If >0, fitness greater in allopatric community

LA_sym <- LA_sym[, c(3,7)]
LA_symAA <- filter(LA_sym, Sp == "A")
count(filter(LA_symAA, diff > 0)) #8 / 12 increase in allo
LA_symPC <- filter(LA_sym, Sp == "P")
count(filter(LA_symPC, diff > 0)) #11 / 12
LA_symVG <- filter(LA_sym, Sp == "V")
count(filter(LA_symVG, diff > 0)) #8 / 12
LA_symOD <- filter(LA_sym, Sp == "O")
count(filter(LA_symOD, diff > 0)) #9 / 12
LA_symSR <- filter(LA_sym, Sp == "S")
count(filter(LA_symSR, diff > 0)) #5 / 12

#Figure 3

label_facets <- function(string){
  len <- length(string)
  string = paste('(', letters[1:len], ') ', string, sep = '')
  return(string)
}

labs <- c('Allopatric', 'Sympatric')

sp.labs2 <- c("Achromobacter sp.", "Ochrobactrum sp.", "Pseudomonas sp.", "Stenotrophomonas sp.", "Variovorax sp.")
sp.labs2 <- as.character(sp.labs2)

LA$Sp <- as.factor(LA$Sp)
levels(LA$Sp) <- sp.labs2
LA$Sp <- as.character(LA$Sp)

LA_means <- data.frame(emmeans::emmeans(LA_mod2, pairwise ~ Treat|Sp))
names(LA_means)[1] <- "Treat"
names(LA_means)[3] <- "fitness"
names(LA_means)[2] <- "Sp"

LA_m2 <- group_by(LA, Treat, Sp) %>%
  summarise(., m.fit = mean(fitness))

LA$Sp2 <- LA$Sp
LA$letter <- letters[as.numeric(as.factor(LA$Sp2))]

LA_means$Sp2 <- LA_means$Sp
LA_means$letter <- letters[as.numeric(as.factor(LA_means$Sp2))]

ggplot() +
  geom_point(data = LA, aes(x = Treat, y = fitness, group = comm), colour = "lightgrey", position = position_dodge(0.2)) +
  geom_line(data = LA, aes(x = Treat, y = fitness, group = comm), colour = "lightgrey", position = position_dodge(0.2)) +
  geom_point(data = LA_means, aes(x = Treat, y = fitness), size = 3) +
  geom_line(data = LA_means, aes(x = Treat, y = fitness, group = Sp)) +
  theme_bw() +
  scale_x_discrete(labels = labs) +
  xlab("Community type") +
  ylab(expression("Growth rate ("~italic(m)~")")) +
  facet_grid(. ~ letter + Sp, labeller = label_bquote(cols = (.(letter))~italic(.(Sp)))) +
  theme(axis.title = element_text(size = 11, colour = "black"), axis.text = element_text(colour = "black", size = 10), strip.background = element_rect(fill="white", colour = "white"), legend.position = "none", strip.text = element_text(hjust = 0, colour = "black", size = 10))


#Table of pairwise comparisons

dat_LA <- LA
dat_LA$Treat <- as.character(dat_LA$Treat)
dat_LA$Treat[dat_LA$Treat == "symp"] <- "Sympatric"
dat_LA$Treat[dat_LA$Treat == "allo"] <- "Allopatric"

dat_LA$Sp <- as.character(dat_LA$Sp)
dat_LA$Sp <- substr(dat_LA$Sp, 1, 1)

LA_mod2 <- lmer(fitness ~ Treat + Sp + (1|clone), data = dat_LA)
LA_reads <- emmeans::emmeans(LA_mod2 , pairwise ~ Treat | Sp)
LA_tab <- data.frame(LA_reads)
LA_tab <- LA_tab[c(1:10), c(1:7)]

LA_tab <- mutate(LA_tab,
                 emmeans.emmean = round(emmeans.emmean, digits = 3),
                 emmeans.SE = round(emmeans.SE, digits = 3),
                 emmeans.df = round(emmeans.df, digits = 3),
                 emmeans.lower.CL = round(emmeans.lower.CL, digits = 3),
                 emmeans.upper.CL = round(emmeans.upper.CL, digits = 3))

LA_tab$emmeans.emmean <- as.character(LA_tab$emmeans.emmean)
LA_tab$emmeans.SE <- as.character(LA_tab$emmeans.SE)
LA_tab$emmeans.df <- as.character(LA_tab$emmeans.df)
LA_tab$emmeans.lower.CL <- as.character(LA_tab$emmeans.lower.CL)
LA_tab$emmeans.upper.CL <- as.character(LA_tab$emmeans.upper.CL)

LA_flex <- flextable(LA_tab) %>%
  set_header_labels(emmeans.emmean = 'Mean', emmeans.Sp = 'Sp.', emmeans.SE = 'SE', emmeans.df = 'df', emmeans.lower.CL = 'Lower CI', emmeans.upper.CL = 'Upper CI', emmeans.Treat = "Community type") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit() %>%
  hline(i = c(2, 4, 6, 8, 10), border = officer::fp_border(color="black")) %>%
  align(j = c(1:7), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

LA_flex