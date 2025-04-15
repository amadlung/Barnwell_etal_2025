###Script for Barnwell et al. 2025 Fig 5 (roots in R)

#load libraries we will need
library(ggplot2)
library (agricolae)
library(dplyr)
library(tidyverse)
library(devtools)
library(ggpubr)
library(stats)

my_data <- read.csv(file.choose(), na.strings=".")  # Barnwell Red Roots B1f - COMBO (1st plus redo) 121524
my_data

#run two-way anova with Genotype and Light_Condition and print results
anova_root = aov(Root_Length ~ Genotype*Light,data=my_data)
summary(anova_root)

###Evaluating Model Assumptions for taproot
lattice::densityplot(anova_root$residuals, auto.key=TRUE)  

#code to create a Q-Q plot
qqnorm(anova_root$residuals)
qqline(anova_root$residuals, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))

plot(anova_root$residuals~anova$fitted.values)
lines(lowess(anova_root$fitted.values,anova_root$residuals), col="blue")


#make interaction term to do Tukey's
Geno_Light_interaction = with(my_data, interaction(Genotype,Light))
interaction_anova = aov(Root_Length ~ Geno_Light_interaction,data=my_data)
summary(interaction_anova)
tukeys = HSD.test(interaction_anova,"Geno_Light_interaction", group = TRUE)
tukeys

#Recode factor levels
my_data$Genotype <- recode(my_data$Genotype, 
                           MM = 'WT',
                           B1  = 'phyB1',
                           E8  = 'phyE-8',
                           F11X = 'phyF-11X',
                           F413 = "phyF-413",
                           EF = "phyEF",
                           B2 = "phyB2",
                           B2E = "phyB2E",
                           B1B2 = "phyB1B2",
                           B1E = "phyB1E",
                           B1F = "phyB1F",
                           B2F = "phyB2F",
                           B1B2E = "phyB1B2E",
                           B1B2F = "phyB1B2F"
)


#to order the x-axis differently: 
my_data$Genotype <- factor(my_data$Genotype, levels=c( "WT", "phyB1", "phyE-8", "phyF-11X", "phyB1E","phyB1F", "phyB2", "phyB1B2", "phyB2E", "phyB2F"))
my_data

###Boxplot
p = #ggplot(germ_data, 
  #    aes(Germination, Genotype ) +
  #         fill=interaction(Germination, Hormone, sep="+", lex.order=TRUE)) + 
  ggplot(my_data, aes(Genotype, Root_Length, Light, colour = Light)) + 
  geom_boxplot(width=0.7) + 
  theme_classic() +
  geom_point(position=position_jitterdodge(0.1)) +
  ylab("Root length [cm]") +
  xlab("Genotype")+
  scale_color_manual(values = c ("black", "red")) +
  theme(legend.position="bottom")

p

p_plus <- p + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(fill="Light condition (D/R)") +
  theme(legend.position="bottom")
p_plus

#plot from above

#Now add the letters from the Tukey test
p_plus_annot <- p_plus + 
  annotate ("text", x = 0.75, y = 8, label = "ef")+  #NOTE: These letters use Tukey data from above
  annotate ("text", x = 1.25, y = 8, label = "a")+ #WT
  
  annotate ("text", x = 1.75, y = 8, label = "cd")+
  annotate ("text", x = 2.25, y = 8, label = "c")+ #B1
  
  annotate ("text", x = 3.75, y = 8, label = "de")+
  annotate ("text", x = 4.25, y = 8, label = "bc")+ #E8
  
  
  annotate ("text", x = 4.75, y = 8, label = "efg")+
  annotate ("text", x = 5.25, y = 8, label = "bc")+ #F11X
  
  
  annotate ("text", x = 5.75, y = 8, label = "cd")+
  annotate ("text", x = 6.25, y = 8, label = "abc")+ #B1E 
  
  annotate ("text", x = 6.75, y = 8, label = "fg")+
  annotate ("text", x = 7.25, y = 8, label = "efg")+ #B1F
  
  
  annotate ("text", x = 2.75, y = 8, label = "ef")+
  annotate ("text", x = 3.25, y = 8, label = "ab")+ #B2
  
  annotate ("text", x = 9.75, y = 8, label = "efg")+
  annotate ("text", x = 10.25, y = 8, label = "efg") +  #B1B2 
 
  annotate ("text", x = 7.75, y = 8, label = "g") +
  annotate ("text", x = 8.25, y = 8, label = "bc")+ #B2E
  
  annotate ("text", x = 8.75, y = 8, label = "efg")+
  annotate ("text", x = 9.25, y = 8, label = "ab") #B2F
  
 

# call annotated graph
p_plus_annot


