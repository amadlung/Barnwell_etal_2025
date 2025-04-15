### Barnwell et al. 2025 Fig.6

library(ggplot2)
library(agricolae)

#open csv file 
ABA_data <- read.csv(file.choose(), na.strings=".") # Barnwell ABA_redroots ALL ABA RED ROOTS 121324
# look at data
ABA_data       

#subset 
ABA_data <- subset(ABA_data, Genotype %in% c("MM", "B1", "B1E", "E8"))
ABA_data

ABA_data$Genotype <- recode(ABA_data$Genotype, 
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

ABA_data


#Change axis order
ABA_data$Genotype <- factor(ABA_data$Genotype,c("WT", "phyB1", "phyE-8", "phyB1E")) 
ABA_data_sub <- ABA_data


##Boxplot

p = ggplot(ABA_data, 
       aes(x = Genotype, y = Root_Length)) + 
  geom_boxplot(aes(fill = interaction(Light, ABA)), 
               position = position_dodge(width = 0.8), 
               width = 0.7, 
               outlier.shape = NA, 
               color = "black") + # 16 boxes grouped by genotype
  geom_point(aes(color = interaction(Light, ABA)), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.5) + # Points jittered and grouped by Genotype and Light × ABA
  scale_fill_manual(values = c("white", "white", "grey90", "grey90")) + # Fill color for Light × ABA combinations
  scale_color_manual(values = c("black", "red", "darkgrey", "darkred")) + # Matching point colors
  ylab("Root Length (cm)") +
  xlab("Genotype") +
  labs(fill = "ABA (yes/no) - Light condition (D/R)", 
       color = "ABA (yes/no) - Light condition (D/R)") + # Legend for color and fill
  theme_classic() + # Clean white background
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        legend.position = "bottom")

p


#Run ANOVA and check assumptions 


ABA_anova <- aov(Root_Length ~ Genotype*ABA*Light, data = ABA_data)
summary(ABA_anova)

layout(matrix(c(1,2,3,4),2,2)) # optional layout
plot(ABA_anova) # diagnostic plots
layout(matrix(c(1),1)) # return to one figure per pane
# check density plot. Should be roughly following a bell-shaped curve.
plot(density(ABA_anova$residuals))

# Run Tukey's Post-Hoc test to see significance between pairwise comparisons, print results
tukey <- TukeyHSD(ABA_anova) # for Lengths, no output yet
tukey

# To get Tukey's post-hoc letter groups, make an interaction term (tx) variable
tx_r <- with(ABA_data, interaction(Genotype, Light, ABA))

#rerun anova with tx term
ABA_anova_int <- aov(Root_Length ~ tx_r, data = ABA_data)

#get and print Tukey's groups (requires Agricolae)
library(agricolae)
letters <- HSD.test (ABA_anova_int,"tx",group= TRUE)
letters

#Now add the letters to the boxplot

#plot from above
p_plus <- p + labs(fill="ABA (yes/no) - Light condition (D/R)") +
  theme(legend.position="bottom")
p_plus

#Now add the letters from the Tukey test
p_plus_annot <- p_plus + 
  annotate ("text", x = 0.75, y = 8, label = "e")+  #NOTE: These letters use Tukey data from above
  annotate ("text", x = 0.9, y = 8, label = "a")+
  annotate ("text", x = 1.1, y = 8, label = "ef")+
  annotate ("text", x = 1.25, y = 8, label = "bc")+
  
  annotate ("text", x = 1.75, y = 8, label = "bc")+
  annotate ("text", x = 1.9, y = 8, label = "bc")+
  annotate ("text", x = 2.1, y = 8, label = "ef")+
  annotate ("text", x = 2.25, y = 8, label = "b")+

  annotate ("text", x = 2.75, y = 8, label = "de")+
  annotate ("text", x = 2.9, y = 8, label = "bc")+
  annotate ("text", x = 3.1, y = 8, label = "f")+
  annotate ("text", x = 3.25, y = 8, label = "cd")+

  annotate ("text", x = 3.75, y = 8, label = "c")+
  annotate ("text", x = 3.9, y = 8, label = "bc")+
  annotate ("text", x = 4.1, y = 8, label = "ef")+
  annotate ("text", x = 4.25, y = 8, label = "efg")


# call annotated graph
p_plus_annot


