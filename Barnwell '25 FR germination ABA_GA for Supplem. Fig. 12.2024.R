### Barnwell et al. 2025 Supplemental Fig.

library(ggplot2)
library(agricolae)
library(dplyr)

#open csv file 
germ_data <- read.csv(file.choose(), na.strings=".") # Barnwell R_FR germination combo 122324
# look at data
germ_data       

# subset variables:
sub_germ_data_ABA <- subset(germ_data, Genotype=='MM' | Genotype == 'B1' | Genotype=='B2' | Genotype == "B1B2" | Genotype == "E8" | Genotype=='F11X' | Genotype == "B1E" | Genotype == "B1F"| Genotype == "EF" | Genotype == "B2E")      
sub_germ_data_ABA

# Dodge width
pd = position_dodge(0.7)

sub_germ_data_ABA$Genotype <- recode(sub_germ_data_ABA$Genotype, 
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

#to order the x-axis differently: 
sub_germ_data_ABA$Genotype <- factor(sub_germ_data_ABA$Genotype, levels=c( "WT", "phyB1", "phyE-8", "phyF-11X", "phyEF", "phyB1E","phyB1F","phyB2", "phyB1B2", "phyB2E"))
sub_germ_data_ABA

# Make plot
p = 
  ggplot(sub_germ_data_ABA, aes(Genotype, GermPercentFR, Hormone, colour = Hormone)) + 
  geom_boxplot(width=0.7, position=pd) + 
  theme_classic() +
  geom_point(position=position_jitterdodge(0.3)) +
  ylab("Germination rate") +
  scale_color_manual(values = c("black","#56B4E9", "#009E73" ), 
                   labels = c("none" = "untreated","with_ABA"  = "ABA treated", "with_GA" = "GA treated")   
)
p

p_plus <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(legend.position="bottom")

p_plus

#Run ANOVA and check assumptions 
combo_anova <- aov(GermPercentFR ~ Genotype*Hormone, data = sub_germ_data_ABA)
summary(combo_anova)

#Assumptions
layout(matrix(c(1,2,3,4),2,2)) # optional layout
plot(combo_anova) # diagnostic plots
layout(matrix(c(1),1)) # return to one figure per pane
# check density plot. Should be roughly following a bell-shaped curve.
plot(density(combo_anova$residuals))

# Run Tukey's Post-Hoc test to see significance between pairwise comparisons, print results
tukey <- TukeyHSD(combo_anova) 
tukey

# To get Tukey's post-hoc letter groups, make an interaction term (tx) variable
tx_r <- with(sub_germ_data_ABA, interaction(Genotype, Hormone))

#rerun anova with tx term
combo_anova_int <- aov(GermPercentFR ~ tx_r, data = sub_germ_data_ABA)

#get and print Tukey's groups (requires Agricolae)
library(agricolae)
letters <- HSD.test (combo_anova_int,"tx",group= TRUE)
letters
# check graph so far
p_plus

#Now add the letters from the Tukey test
p_plus_annot <- p_plus + 
  annotate ("text", x = 0.75, y =-0.03, label = "h")+  #NOTE: These letters use Tukey data from above
  annotate ("text", x = 1, y = -0.07, label = "h")+
  annotate ("text", x = 1.25, y = -0.03, label = "gh")+ #WT
  
  annotate ("text", x = 1.75, y = -0.03, label = "def")+
  annotate ("text", x = 2, y = -0.07, label = "fgh")+
  annotate ("text", x = 2.25, y = -0.03, label = "cd")+ #B1
  
  annotate ("text", x = 2.75, y = -0.03, label = "cd")+
  annotate ("text", x = 3, y = -0.07, label = "gh")+
  annotate ("text", x = 3.25, y = -0.03, label = "cd")+ #E8
  
  annotate ("text", x = 3.75, y = -0.03, label = "efgh")+ 
  annotate ("text", x = 4, y = -0.07, label = "h")+
  annotate ("text", x = 4.25, y = -0.03, label = "efgh")+ #F11X
  
  annotate ("text", x = 4.75, y = -0.03, label = "de")+
  annotate ("text", x = 5, y = -0.07, label = "h")+
  annotate ("text", x = 5.25, y = -0.03, label = "def")+ #EF
  
  annotate ("text", x = 5.75, y = -0.03, label = "a")+
  annotate ("text", x = 6, y = -0.07, label = "bc")+
  annotate ("text", x = 6.25, y = -0.03, label = "a")+ #B1E
    
  annotate ("text", x = 6.75, y = -0.03, label = "def")+
  annotate ("text", x = 7.25, y = -0.03, label = "h")+  #B1F
    
  annotate ("text", x = 7.75, y = -0.03, label = "h")+
  annotate ("text", x = 8, y = -0.07, label = "gh")+
  annotate ("text", x = 8.25, y = -0.03, label = "h")+ #B2
    
  annotate ("text", x = 8.75, y = -0.03, label = "b")+
  annotate ("text", x = 9, y = -0.07, label = "h")+
  annotate ("text", x = 9.25, y = -0.03, label = "defg")+ #B1B2

  annotate ("text", x = 9.75, y = -0.03, label = "fgh")+
  annotate ("text", x = 10, y = -0.07, label = "h")+
  annotate ("text", x = 10.25, y = -0.03, label = "gh") #B2E
  

# call annotated graph
p_plus_annot


letters$means
letters$groups

########################################################################

