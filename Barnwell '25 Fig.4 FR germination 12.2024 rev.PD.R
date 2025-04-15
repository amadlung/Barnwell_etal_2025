### Using R-studio for stats analysis
# 2-way ANOVA 

library(ggplot2)
library(agricolae)
library(dplyr)

#open csv file 
germ_data <- read.csv(file.choose(), na.strings=".") # Barnwell R_FR germination combo 122324
# look at data
germ_data       


#Subset:
sub_germ_data_ABA <- subset(germ_data, Genotype=='MM' | Genotype == 'B1' | Genotype=='B2' | Genotype == "B1B2" | Genotype == "E8" | Genotype=='F11X' | Genotype == "B1E" | Genotype == "B1F"| Genotype == "EF" | Genotype == "B2E" | Hormone == 'ABA' | Hormone == 'none')      
sub_germ_data_ABA <- subset(
  germ_data, 
  (Genotype == 'MM' | Genotype == 'B1' | Genotype == 'B2' | Genotype == "B1B2" | Genotype == "E8" | Genotype == 'F11X' | Genotype == "B1E" | Genotype == "B1F" | Genotype == "EF" | Genotype == "B2E") & 
    (Hormone == 'with_ABA' | Hormone == 'none')
)

sub_germ_data_ABA

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
p_ABA = #ggplot(germ_data, 
  #    aes(Germination, Genotype ) +
  #         fill=interaction(Germination, Hormone, sep="+", lex.order=TRUE)) + 
  ggplot(sub_germ_data_ABA, aes(Genotype, GermPercentFR, Hormone, colour = Hormone)) + 
  geom_boxplot(width=0.7) + 
  theme_classic() +
  geom_point(position=position_jitterdodge(0.3)) +
  ylab("Germination rate") +
  scale_color_manual(values = c("black","#0072B2"),
                     labels = c("none" = "untreated","with_ABA"  = "ABA treated")   
  )

p_ABA

p_plus_ABA <- p_ABA + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme(legend.position="bottom")

p_plus_ABA

#Run ANOVA and check assumptions 

combo_anova_ABA <- aov(GermPercentFR ~ Genotype*Hormone, data = sub_germ_data_ABA)
summary(combo_anova_ABA)

#Assumptions
layout(matrix(c(1,2,3,4),2,2)) # optional layout
plot(combo_anova_ABA) # diagnostic plots
layout(matrix(c(1),1)) # return to one figure per pane
# check density plot. Should be roughly following a bell-shaped curve.
plot(density(combo_anova_ABA$residuals))

# Run Tukey's Post-Hoc test to see significance between pairwise comparisons, print results
tukey_ABA <- TukeyHSD(combo_anova_ABA) 
tukey_ABA

# To get Tukey's post-hoc letter groups, make an interaction term (tx) variable
tx_r_ABA <- with(sub_germ_data_ABA, interaction(Genotype, Hormone))

#rerun anova with tx term
combo_anova_int_ABA <- aov(GermPercentFR ~ tx_r_ABA, data = sub_germ_data_ABA)

#get and print Tukey's groups (requires Agricolae)
library(agricolae)
letters_ABA <- HSD.test (combo_anova_int_ABA,"tx",group= TRUE)
letters_ABA


#Now add the letters from the Tukey test
p_plus_annot_ABA <- p_plus_ABA + 
  annotate ("text", x = 0.75, y =-0.02, label = "g")+  #NOTE: These letters use Tukey data from above
  annotate ("text", x = 1.25, y =-0.02, label = "g")+
  
  annotate ("text", x = 1.75, y =-0.02, label = "def")+
  annotate ("text", x = 2.25, y =-0.02, label = "fg")+
  
  annotate ("text", x = 2.75, y =-0.02, label = "cd")+
  annotate ("text", x = 3.25, y =-0.02, label = "g")+
 
  annotate ("text", x = 3.75, y =-0.02, label = "efg")+
  annotate ("text", x = 4.25, y =-0.02, label = "g")+
 
  annotate ("text", x = 4.75, y =-0.02, label = "de")+
  annotate ("text", x = 5.25, y =-0.02, label = "g")+
 
  annotate ("text", x = 5.75, y =-0.02, label = "a")+ 
  annotate ("text", x = 6.25, y =-0.02, label = "bc")+
 
  annotate ("text", x = 6.75, y =-0.02, label = "def")+
  annotate ("text", x = 7.25, y =-0.02, label = "g")+
  
  annotate ("text", x = 7.75, y =-0.02, label = "g")+
  annotate ("text", x = 8.25, y =-0.02, label = "g")+
  
  annotate ("text", x = 8.75, y =-0.02, label = "b")+
  annotate ("text", x = 9.25, y =-0.02, label = "g")+
  
  annotate ("text", x = 9.75, y =-0.02, label = "g")+
  annotate ("text", x = 10.25, y =-0.02, label = "g")
 

# call annotated graph
p_plus_annot_ABA

# call letters and means
letters_ABA$means
letters_ABA$groups




