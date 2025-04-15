###Script for Barnwell et al., 2025 (Fig 3)

#load libraries 
library(ggplot2)
library (agricolae)
library(dplyr)
library(tidyverse)
library(devtools)
library(ggpubr)
library(stats)

my_data <- read.csv(file.choose(), na.strings=".")  # Barnwell phyE shoots Cr master sheet Sheet 2 121424 (note that B1 data is the same as from Balderrama, 2023 Fig 2)
my_data


# subset variables:
my_data <- subset(my_data, genotype=='MM' | genotype == 'B1' | genotype=='B2' | genotype == "B1B2" | genotype == "E8" | genotype=='F11' | genotype == "B1E" | genotype == "B1F"| genotype == "EF" | genotype == "B2E" | genotype == "F413" | genotype == "B2F" | genotype == "B1B2E" | genotype == "B1B2F")      
my_data

# Dodge width
pd = position_dodge(0.7)

my_data$genotype <- recode(my_data$genotype, 
                                     MM = 'WT',
                                     E8  = 'phyE-8',
                                     B1E = 'phyB1E',
                                     B1  = 'phyB1',
                                     F11X = 'phyF-11X',
                                     B2 = 'phyB2',
                                     B1B2 = 'phyB1B2',
                                     B2E = 'phyB2E',
                                     EF = 'phyEF',
                                     B1F = 'phyB1F')



#run two-way anova with Genotype and Light_Condition and print results
anova = aov(length ~ genotype*light,data=my_data)
summary(anova)

###Evaluating Model Assumptions 
lattice::densityplot(anova$residuals, auto.key=TRUE)  

#code to create a Q-Q plot
qqnorm(anova$residuals)
qqline(anova$residuals, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))

plot(anova$residuals~anova$fitted.values)
lines(lowess(anova$fitted.values,anova$residuals), col="blue")

#make interaction term to do Tukey's
Geno_Light_interaction = with(my_data, interaction(genotype,light))
interaction_anova = aov(length ~ Geno_Light_interaction,data=my_data)
summary(interaction_anova)
tukeys = HSD.test(interaction_anova,"Geno_Light_interaction", group = TRUE)
tukeys

#recode factor levels, needs dplyr already called above
my_data$genotype <- recode(my_data$genotype, 
                           MM = 'WT',
                           B1  = 'phyB1',
                           E8  = 'phyE-8',
                           F11 = 'phyF-11X',
                           F413 = "phyF-413",
                           EF = "phyEF",
                           B2 = "phyB2",
                           B1B2 = "phyB1B2",
                           B1E = "phyB1E",
                           B1F = "phyB1F",
                           B2F = "phyB2F",
                           B1B2E = "phyB1B2E",
                           B1B2F = "phyB1B2F"
                           )

#to order the x-axis differently: 
my_data$genotype <- factor(my_data$genotype, levels=c('WT', 'phyB1', "phyE-8", "phyF-11X", 'phyF-413', 'phyEF',  "phyB1E", "phyB1F",  "phyB2", "phyB1B2", "phyB2F",  'phyB1B2E', "phyB1B2F"))


###Boxplot
p = ggplot(my_data, aes(genotype, length, light, colour = light)) + 
  geom_boxplot(width=0.7) + 
  theme_classic() +
  geom_point(position=position_jitterdodge(0.1)) +
  ylab("Hypocotyl length ([cm]") +
  xlab("Genotype")+
  scale_color_manual(values = c ("grey30", "red")) +
  theme(legend.position="bottom")
p

p_plus <- p + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(fill="Light condition (D/R)") +
  theme(legend.position="bottom")
p_plus

#Now add the letters from the Tukey test
p_plus_annot <- p_plus + 
  annotate ("text", x = 0.75, y = 15, label = "cde")+  #NOTE: These letters use Tukey data calculated above
  annotate ("text", x = 1.25, y = 14, label = "j")+ #WT
  
  annotate ("text", x = 1.75, y = 15, label = "ab")+
  annotate ("text", x = 2.25, y = 14, label = "defg")+ #B1
  
  annotate ("text", x = 2.75, y = 15.2, label = "abc")+
  annotate ("text", x = 3.25, y = 14, label = "ghi")+ #E8
  
  annotate ("text", x = 3.75, y = 15, label = "abcd")+
  annotate ("text", x = 4.25, y = 14, label = "gh")+ #F11X
  
  annotate ("text", x = 4.75, y = 15, label = "ab")+
  annotate ("text", x = 5.25, y = 14, label = "ij")+ #F413
  
  annotate ("text", x = 5.75, y = 15, label = "abcde")+
  annotate ("text", x = 6.25, y = 14, label = "gh")+ #EF
  
  annotate ("text", x = 6.75, y = 15.2, label = "a")+
  annotate ("text", x = 7.25, y = 14.5, label = "a")+ #B1E
  
  annotate ("text", x = 7.75, y = 15, label = "cdef") +
  annotate ("text", x = 8.25, y = 14, label = "abcde")+ #B1F
  
  annotate ("text", x = 8.75, y = 15, label = "fg")+
  annotate ("text", x = 9.25, y = 14, label = "hij")+ #B2 

  annotate ("text", x = 9.75, y = 15, label = "abcde")+
  annotate ("text", x = 10.25, y = 14, label = "abcd")+  #B1B2
  
  annotate ("text", x = 10.75, y = 15, label = "a")+
  annotate ("text", x = 11.25, y = 14, label = "efg")+  #B2F
  
  annotate ("text", x = 11.75, y = 15, label = "fgh") +
  annotate ("text", x = 12.25, y = 14, label = "cdefg")+  #B1B2E
  
  annotate ("text", x = 12.75, y = 15, label = "bcde")+
  annotate ("text", x = 13.25, y = 14, label = "abc")  #B1B2F
  
# call annotated graph
p_plus_annot


