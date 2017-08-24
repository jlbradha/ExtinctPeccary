# Created by Jen Bradham on January 12, 2017
# Script to plot and statistically evaluate extinct peccary microwear and isotopes and compare them with modern peccary microwear and isotops

rm(list=ls())

library(tidyverse)
library(RColorBrewer)
library(dunn.test)


setwd("~/Dropbox/PhD_Bradham/ExtinctPeccaries/")

#isotopes <- read.csv('PeccaryIsotopes_R.csv', header = TRUE)
isotopes <- read.csv('PeccaryIsotopes__LifeStatusR.csv', header = TRUE)
#microwear <- read.csv('PeccaryMicrowear_R.csv', header = TRUE)
microwear <- read.csv('PeccaryMicrowear_R_LifeStatus.csv', header = TRUE)
micro_isotopes <- read.csv('PeccaryMicroIsotopes_R.csv', header = TRUE)
#names(microwear) 

# Test for Normality
isotopes <- isotopes[! (isotopes$Genus == "Collared"),] #remove collared peccaries
shapiro.test(isotopes$Delta13C)
qqnorm(isotopes$Delta13C)

normal_mylo <- filter(isotopes, Genus == 'Mylohyus')
normal_platy <- filter(isotopes, Genus == 'Platygonus')
normal_catag <- filter(isotopes, Genus == 'Catagonus')
normal_WLP <- filter(isotopes, Genus == 'WLP')

shapiro.test(normal_mylo$Delta13C)
shapiro.test(normal_platy$Delta13C)
shapiro.test(normal_catag$Delta13C)
shapiro.test(normal_WLP$Delta13C)

qqnorm(normal_mylo$Delta13C)
qqnorm(normal_platy$Delta13C)
qqnorm(normal_catag$Delta13C)
qqnorm(normal_WLP$Delta13C)


shapiro.test(microwear$Asfc)
shapiro.test(microwear$epLsar)
shapiro.test(microwear$Tfv)
shapiro.test(microwear$X3x3HAsfc)
shapiro.test(microwear$X9x9HAsfc)


# Define Functions
MannWhitney <- function (data1, data2) {
  wilcox.test(data1, data2, paired = FALSE)$p.value
}

## ---------------- Isotopes ------------------ ##
isotopes <- isotopes[! (isotopes$Genus == "Collared"),] #remove collared peccaries

### CARBON BOXPLOTS
#### Boxplot of delta13C per EXTINCT species
# # Species on y axis
# isotopes_extinct <- filter (isotopes, Status == 'Extinct')
# Extinct_carbon <- ggplot(isotopes_extinct, aes(x=Genus, y=Delta13C)) + geom_boxplot(aes(fill=NALMA)) + coord_flip() 
# Extinct_carbon + ylab(delta^13~C) + xlab("Species") 
# 
# # Time period on y axis
# isotopes_extinct <- filter (isotopes, Status == 'Extinct') 
# Extinct_carbon <- ggplot(isotopes_extinct, aes(x=NALMA, y=Delta13C)) + geom_boxplot(aes(fill=Genus)) + coord_flip() 
# Extinct_carbon + ylab(delta^13~C) + xlab("NALMA") 

#Time period on y axis in correct temporal order for EXTINCT peccaries
# isotopes_extinct <- filter (isotopes, Status == 'Extinct')
# isotopes_extinct$NALMA <- factor(isotopes_extinct$NALMA,
#                                  levels = c('Hemphillian', 'Blancan', 'Irvingtonian',                                   'Rancholabrean'), ordered = TRUE)
# Extinct_carbon <- ggplot(isotopes_extinct, aes(x=NALMA, y=Delta13C)) + geom_boxplot(aes(fill=Genus)) + coord_flip() 
# Extinct_carbon + ylab(delta^13~C) + xlab("NALMA") 
# 



##### ~%~%~%~ ##### ~%~%~%~ ##### ~%~%~%~ ##### ~%~%~%~ ##### ~%~%~%~ ##### ~%~%~%~
##### ~%~%~%~ ##### ~%~%~%~ ##### ~%~%~%~ ##### ~%~%~%~ ##### ~%~%~%~ ##### ~%~%~%~

#### BOXPLOT OF DELTA13C PER NALMA FOR ALL SPECIES 

isotopes$NALMA <- factor(isotopes$NALMA,
                                 levels = c('Hemphillian', 'Blancan', 'Irvingtonian',                                   'Rancholabrean', 'Modern'), ordered = TRUE)
species_carbon <- ggplot(isotopes, aes(x=NALMA, y=Delta13C_Seuss)) + geom_boxplot(aes(fill=Genus)) + coord_flip()
species_carbon +  ylab(delta^13~C) + theme_classic() + scale_y_continuous(limits=c(-16, 0)) +
  scale_fill_manual(values = c("#FFFFFF", "#666666", "#000000", "#CCCCCC"),
                    labels = c("Catagonus sp.", "Mylohyus sp.", "Platygonus sp.", "Tayassu pecari")) +
  theme(plot.background = element_blank(), #remove the gray background and grid
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(), #remove tick marks on axes
        axis.ticks.y=element_blank(), 
        axis.title.y=element_blank()) +
  theme(panel.border= element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 1), #drawing the x and y axis lines
        axis.line.y = element_line(color="black", size = 1)) +
  theme(legend.text = element_text(face="italic")) + #Italisize legend words
  theme (axis.text = element_text(size = 15), 
         axis.title = element_text(size = 17, face = "bold"), 
         legend.title = element_text(size = 15), 
         legend.text = element_text(size = 14), 
         legend.key.size = unit(1.5, 'lines')) 


#Colors: FFFFFF=white , 666666=darker gray, 000000"=black, CCCCCC = light gray








#### Carbon Boxplot for Mylo and Catag in Hemph
Hemph_carbon <- filter(isotopes, NALMA == 'Hemphillian')
Hemph_carbon_plot <- ggplot(Hemph_carbon, aes(x=Genus, y=Delta13C)) + geom_boxplot() + coord_flip() + ggtitle('Hemphilian: Mylohyus and Catagonus Carbon isotopes') 
Hemph_carbon_plot

Mylo_Hemph_isotopes <- filter(Hemph_carbon, Genus == 'Mylohyus')
Catag_Hemph_isotopes <- filter(Hemph_carbon, Genus == 'Catagonus')
wilcox.test(Mylo_Hemph_isotopes$Delta13C,Catag_Hemph_isotopes$Delta13C,paired=FALSE,alternative="greater")

#### Carbon Boxplot for Mylo and Platygonus in Irv (only NALMA co-occur with large n)
Irv_carbon <- filter(isotopes, NALMA == 'Irvingtonian')
Irv_carbon_plot <- ggplot(Irv_carbon, aes(x=Genus, y=Delta13C)) + geom_boxplot() + coord_flip() + ggtitle('Irvingtonian: Mylohyus and Platygonus Carbon isotopes') 
Irv_carbon_plot

######## ----------- %%%%%%%%%% ---------- ######### 
# P VALUES ARE DIFFERENT DEPENDING ON WHAT VERSTION OF WILCOX TEST YOU DO: 
# see: http://www.statmethods.net/stats/nonparametric.html
# Wilcox or Mann Whitney? paired = false (like we did in the function above) indicates a Mann Whitney test

# independent but specifying one tail
wilcox.test(Mylo_Irv_isotopes$Delta13C,Platy_Irv_isotopes$Delta13C,paired=FALSE,alternative="greater") #alternative = greater specifies one tailed test
# paired = false says the two are independent

#independent 2-group Mann-Whitney U test
wilcox.test(Mylo_Irv_isotopes$Delta13C,Platy_Irv_isotopes$Delta13C)
## Malu - do 2 tailed test. If you expect one should have lower c values than other, then do a 1 tail. She says do 2-tailed for everything because this is more conservative. One tailed already says we have a reason to believe these are different. Two-tailed is the way to go becuase we aren't pre-expecting a difference. If figure looks like there is a difference but two-tailed doesn't show, then run 1 and if it's significant then maybe there is a slight difference. 
######## ----------- %%%%%%%%%% ---------- #########


### CARBON STATS


# WLPC_seuss <- isotopes %>% filter(Genus == 'WLP') %>% transmute(Delta13C, Delta13C = Delta13C + 1.5) 
# WLPC_seuss <- list(WLPC_seuss)
# 
# isotopes <- isotopes %>% mutate(Delta13C = 
#   if (isotopes$Genus == 'WLP') {
#   isotopes$Delta13C[i] <- WLPC_seuss[[i]]
# })
# 
# #replacing the original WLP Carbon column with the new values calculated according to the Seuss effect above



# Change the column headings so the column in the spreadsheet labeled "delta13C_seuss" is now the delta13C column (I already had the code written with "delta13C" so I didn't want to rewrite that to add 'seuss' to it every time. 
colnames(isotopes) <- c('Number', 'Genus', 'Site', 'NALMA', 'Date', 'Delta13C_original', 'Delta18O', 'Status', 'Delta13C')
  
### Organizing isotope date into groups

mylo_iso <- filter(isotopes, Genus == 'Mylohyus') 
platy_iso <- filter(isotopes, Genus == 'Platygonus') 
catag_iso <- filter(isotopes, Genus == 'Catagonus') 
wlp_iso <- filter(isotopes, Genus == 'WLP') 

mylo_iso_hemph <- filter(mylo_iso, NALMA == 'Hemphillian') 
mylo_iso_irv <- filter(mylo_iso, NALMA == 'Irvingtonian') 
mylo_iso_blanc <- filter(mylo_iso, NALMA == 'Blancan') 
mylo_iso_ranch <- filter(mylo_iso, NALMA == 'Rancholabrean') 

platy_iso_irv <- filter(platy_iso, NALMA == 'Irvingtonian')
platy_iso_blanc <- filter(platy_iso, NALMA == 'Blancan')
platy_iso_ranch <- filter(platy_iso, NALMA == 'Rancholabrean')

# iso_hemph <- filter(isotopes, NALMA == 'Hemphillian')
# iso_irv <- filter(isotopes, NALMA == 'Irvingtonian')
# iso_blanc <- filter(isotopes, NALMA == 'Blancan')
# iso_ranch <- filter(isotopes, NALMA == 'Rancholabrean')


# MannWhitney tests only work for comparing two species. If comparing more species or if comparing one species through time (which counts as a comparison of more than 2 entities), then need to do Kruskal Wallis test using Dunn's Test. Also known as a Kruskal-Wallis rank sum test?



## Sum Stats per Genus isotopes
SumStats_Iso.Genus <- isotopes %>% group_by(Genus) %>% summarise(
  avg_Delta13C = mean(Delta13C),
  med_Delta13C = median(Delta13C),
  sd_Delta13C = sd(Delta13C),
  max_Delta13C = max(Delta13C),
  min_Delta13C = min(Delta13C),
  range_Delta13C = max_Delta13C - min_Delta13C
)

SumStats_Iso.Genus <- t(SumStats_Iso.Genus)
write.csv(SumStats_Iso.Genus, 'Isotopes_SumStats_Genus.csv')

## Sum Stats per Genus per NALMA
SumStats_Iso.GenusNALMA <- isotopes %>% group_by(Genus, NALMA) %>% summarise(
  avg_Delta13C = mean(Delta13C),
  med_Delta13C = median(Delta13C),
  sd_Delta13C = sd(Delta13C),
  max_Delta13C = max(Delta13C),
  min_Delta13C = min(Delta13C),
  range_Delta13C = max_Delta13C - min_Delta13C
)

SumStats_Iso.GenusNALMA <- t(SumStats_Iso.GenusNALMA)
write.csv(SumStats_Iso.GenusNALMA, 'Isotopes_SumStats_Genus-NALMA.csv')

## Sum Stats per NALMA (clumped Genera)
SumStats_Iso.NALMA <- isotopes %>% group_by(NALMA) %>% summarise(
  avg_Delta13C = mean(Delta13C),
  med_Delta13C = median(Delta13C),
  sd_Delta13C = sd(Delta13C),
  max_Delta13C = max(Delta13C),
  min_Delta13C = min(Delta13C),
  range_Delta13C = max_Delta13C - min_Delta13C
)

SumStats_Iso.NALMA <- t(SumStats_Iso.NALMA)
write.csv(SumStats_Iso.NALMA, 'Isotopes_SumStats_NALMA.csv')


## Calculating P values between Genera for each isotope analysis

# # Isotope p-values between genera (not divided by time period)
# Isotopes_Genus <- matrix(NA, nrow=2, ncol=2)
# rownames(Isotopes_Genus) <- c('Mylohyus', 'Platygonus')
# colnames(Isotopes_Genus) <- c('Platygonus', 'Catagonus')
# 
# Isotopes_Genus[1,1] <- MannWhitney(mylo_iso$Delta13C, platy_iso$Delta13C)
# Isotopes_Genus[1,2] <- MannWhitney(mylo_iso$Delta13C, catag_iso$Delta13C)
# Isotopes_Genus[2,2] <- MannWhitney(platy_iso$Delta13C, catag_iso$Delta13C)
# 
# write.csv(Isotopes_Genus, 'Isotopes-Genera.csv')

# # Isotope p-values between time periods (not divided by genera)
# Isotope_NALMAS <- matrix(NA, nrow=3, ncol=3)
# rownames(Isotope_NALMAS) <- c('Hemphillian', 'Irvingtonian', 'Blancan')
# colnames(Isotope_NALMAS) <- c('Irvingtonian', 'Blancan', 'Rancholabrean')
#   
# Isotope_NALMAS[1,1]<- MannWhitney(iso_hemph$Delta13C, iso_irv$Delta13C)
# Isotope_NALMAS[1,2]<- MannWhitney(iso_hemph$Delta13C, iso_blanc$Delta13C)
# Isotope_NALMAS[1,3]<- MannWhitney(iso_hemph$Delta13C, iso_ranch$Delta13C)
# Isotope_NALMAS[2,2]<- MannWhitney(iso_irv$Delta13C, iso_blanc$Delta13C)
# Isotope_NALMAS[2,3]<- MannWhitney(iso_irv$Delta13C, iso_ranch$Delta13C)
# Isotope_NALMAS[3,3]<- MannWhitney(iso_blanc$Delta13C, iso_ranch$Delta13C)
# 
# write.csv(Isotope_NALMAS, 'Isotopes-NALMAS.csv')


# #Isotope p-values Mylohyus through time

#Kruskal Wallis with Dunn test evaluating Mylohyus delta13c through time
isotope_Mylo.time <- dunn.test(list(mylo_iso_ranch$Delta13C, mylo_iso_irv$Delta13C, mylo_iso_blanc$Delta13C, mylo_iso_hemph$Delta13C), method = "none", kw=TRUE, label = TRUE)
#Group 1 = Mylo rancholabrean; group 2=mylo irvingtonian; group 3=mylo blancan; group 4 = mylo hemphillian

# Matrix_Mylo_Iso <- matrix(NA, nrow=3, ncol=3)
# rownames(Matrix_Mylo_Iso) <- c('Hemphillian', 'Irvingotnian', 'Blancan')
# colnames(Matrix_Mylo_Iso) <- c('Irvingtonian', 'Blancan', 'Rancholabrean')
# 
# Matrix_Mylo_Iso[1,1]<- MannWhitney(mylo_iso_hemph$Delta13C, mylo_iso_irv$Delta13C)
# Matrix_Mylo_Iso[1,2]<- MannWhitney(mylo_iso_hemph$Delta13C, mylo_iso_blanc$Delta13C)
# Matrix_Mylo_Iso[1,3]<- MannWhitney(mylo_iso_hemph$Delta13C, mylo_iso_ranch$Delta13C)
# Matrix_Mylo_Iso[2,2]<- MannWhitney(mylo_iso_irv$Delta13C, mylo_iso_blanc$Delta13C)
# Matrix_Mylo_Iso[2,3]<- MannWhitney(mylo_iso_irv$Delta13C, mylo_iso_ranch$Delta13C)
# Matrix_Mylo_Iso[3,3]<- MannWhitney(mylo_iso_blanc$Delta13C, mylo_iso_ranch$Delta13C)
# 
# write.csv(Matrix_Mylo_Iso, 'Isotopes_Mylo_throughtime.csv')


# #Isotope p-values Platygonus through time with Dunn test

isotope_platy.time <- dunn.test(list(platy_iso_ranch$Delta13C, platy_iso_irv$Delta13C, platy_iso_blanc$Delta13C), method = "none", kw=TRUE, label = TRUE, table = TRUE)
#Group 1 = platy rancho; group2=platy irv; group3=platy blanc


# Matrix_platy_iso <- matrix(NA, nrow=2, ncol=2)
# colnames(Matrix_platy_iso) <- c("Blancan", "Rancholabrean")
# rownames(Matrix_platy_iso) <- c("Irvingtonian", "Blancan")
# 
# Matrix_platy_iso[1,1] <- MannWhitney(platy_iso_irv$Delta13C, platy_iso_blanc$Delta13C)
# Matrix_platy_iso[1,2] <- MannWhitney(platy_iso_irv$Delta13C, platy_iso_ranch$Delta13C)
# Matrix_platy_iso[2,2] <- MannWhitney(platy_iso_blanc$Delta13C, platy_iso_ranch$Delta13C)
# 
# write.csv(Matrix_platy_iso, 'Isotopes_Platy_throughtime.csv')


# Comparison of all peccaries exclusive of NALMA with Dunn test
isotope_mylo.platy.catag <- dunn.test(list(mylo_iso$Delta13C, platy_iso$Delta13C, catag_iso$Delta13C, wlp_iso$Delta13C), method="none", kw=TRUE, label = TRUE)
#group1 = Mylohyus, group2=platygonus, group3=catagonus, group 4 = WLP



#Isotope p-values Mylo vs Catag in Hemphillian
MannWhitney(mylo_iso_hemph$Delta13C, catag_iso$Delta13C) # p= 0.03212329

#Isotope p-values Mylo vs Platy through time
Iso_mylo.platy <- matrix(NA, nrow=1, ncol=3)
colnames(Iso_mylo.platy) <- c('Irvingtonian', 'Blancan', 'Rancholabrean')

Iso_mylo.platy[,1] <- MannWhitney(mylo_iso_irv$Delta13C, platy_iso_irv$Delta13C) 
#0.5174734
Iso_mylo.platy[,2] <- MannWhitney(mylo_iso_blanc$Delta13C, platy_iso_blanc$Delta13C) 
#0.1835734
Iso_mylo.platy[,3] <- MannWhitney(mylo_iso_ranch$Delta13C, platy_iso_ranch$Delta13C) 
#0.002204394 --> only significant one

write.csv(Iso_mylo.platy, 'Isotopes_mylo-platy_NALMA.csv')

# Isotope p-vales WLP vs each of the extinct species
## NO SHOULD BE DUNN TEST
Iso_WLP.extinct <- matrix(NA, nrow=1, ncol=3)
colnames(Iso_WLP.extinct) <- c('Catagonus', 'Mylohyus', 'Platygonus')
rownames(Iso_WLP.extinct) <- c('WLP')

Iso_WLP.extinct [1,1] <- MannWhitney(wlp_iso$Delta13C, catag_iso$Delta13C)
Iso_WLP.extinct [1,2] <- MannWhitney(wlp_iso$Delta13C, mylo_iso$Delta13C)
Iso_WLP.extinct [1,3] <- MannWhitney(wlp_iso$Delta13C, platy_iso$Delta13C)

write.csv(Iso_WLP.extinct, 'Isotopes_WLP-Extinct.csv')







## ---------------- Microwear ------------------ ##

microwear <- microwear[! (microwear$Genus == "Collared"),] #remove collared peccaries

###############


micro_scatter <- ggplot(data=microwear, aes(x=Asfc, y=epLsar))
micro_scatter + geom_point(aes(shape = Genus, color = Genus), size = 4) +  #4 is size of points
  theme_classic() +
  scale_color_manual(values = c("#CCCCCC", "#666666", "#000000", "#999999")) + 
  theme(plot.background = element_blank(), #remove the gray background and grid
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(), #remove tick marks on axes
        axis.ticks.y=element_blank()) +
  theme(panel.border= element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 1), #drawing the x and y axis lines
        axis.line.y = element_line(color="black", size = 1)) +
  theme(legend.text = element_text(face="italic")) + #Italisize legend words
  theme (axis.text = element_text(size = 15), 
         axis.title = element_text(size = 17, face = "bold"), 
         legend.title = element_text(size = 15), 
         legend.text = element_text(size = 14), 
         legend.key.size = unit(1.5, 'lines')) 


#Colors: FFFFFF=white , 666666=darker gray, 000000"=black, CCCCCC = light gray
                                                                                                                   
                                                                                                                   
#Creating Polygons around all_species scatterplot...with Jonathan's Help!
hulls <- na.omit(microwear) %>% 
  group_by(Genus) %>% 
  do(hull = chull(.$Asfc, .$epLsar)) %>%
  ungroup()

micro_hulls = na.omit(microwear) %>% 
  left_join(hulls, by = 'Genus') %>%
  group_by(Genus) %>% 
  mutate(index = 1:n()) %>%
  ungroup() %>%
  rowwise() %>%
  filter(index %in% hull) %>%
  mutate(h_index = which(index == hull)) %>%
  ungroup() %>%
  select(Genus, index, h_index, Asfc, epLsar, hull) %>%
  arrange(Genus, h_index)

Micro_species <- ggplot(micro_hulls, aes(Asfc, epLsar, color = Genus, fill = Genus)) +
  geom_polygon(data = filter(micro_hulls, !is.na(hull)), alpha = 0.5) +
  geom_point() 
Micro_species + ggtitle('Asfc vs epLsar for all peccaries')


# Scatterplot of Asfc-epLsar per EXTINCT species
Extinct_microwear <- filter(microwear, Status == 'Extinct')
Extinct_micro <- qplot(Asfc, epLsar, data=Extinct_microwear, color=Genus)
Extinct_micro

###############





######### %%%%%%% ######### %%%%%%%%% ############

# Grouping Data by Genus
micro_mylo <- filter(microwear, Genus == 'Mylohyus')
micro_platy <- filter(microwear, Genus == 'Platygonus')
micro_catag <- filter(microwear, Genus == 'Catagonus')
micro_WLP <- filter(microwear, Genus == 'T. pecari')

# # Grouping data by Time period
# micro_Rancho <- filter(microwear, NALMA == 'Rancholabrean')
# micro_Irv <- filter(microwear, NALMA == 'Irvingtonian')
# micro_Blanc <- filter(microwear, NALMA == 'Blancan') 
# micro_Hemph <- filter(microwear, NALMA == 'Hemphillian') 

# # Grouping data by Genus and Time Period
mylo_hemph <- filter(micro_mylo, NALMA == 'Hemphillian')
mylo_blanc <- filter(micro_mylo, NALMA == 'Blancan')
mylo_irv <- filter(micro_mylo, NALMA == 'Irvingtonian')
mylo_ranch <- filter(micro_mylo, NALMA == 'Rancholabrean')

platy_blanc <- filter(micro_platy, NALMA == 'Blancan')
platy_irv <- filter(micro_platy, NALMA == 'Irvingtonian')
platy_ranch <- filter(micro_platy, NALMA == 'Rancholabrean')

catag_hemph <- filter(micro_catag, NALMA == 'Hemphillian')

######### %%%%%%% ######### %%%%%%%%% ############


# Microwear between all peccaries exclusive of NALMA with Dunn Test

All_Asfc <- dunn.test(list(micro_mylo$Asfc, micro_platy$Asfc, micro_catag$Asfc, micro_WLP$Asfc), method = "none", kw=TRUE, label = TRUE)
All_epLsar <- dunn.test(list(micro_mylo$epLsar, micro_platy$epLsar, micro_catag$epLsar, micro_WLP$epLsar), method = "none", kw=TRUE)
All_Tfv <- dunn.test(list(micro_mylo$Tfv, micro_platy$Tfv, micro_catag$Tfv, micro_WLP$Tfv), method = "none", kw=TRUE)
All_X3x3HAsfc <- dunn.test(list(micro_mylo$X3x3HAsfc, micro_platy$X3x3HAsfc, micro_catag$X3x3HAsfc, micro_WLP$X3x3HAsfc), method = "none", kw=TRUE)
All_X9x9HAsfc <- dunn.test(list(micro_mylo$X9x9HAsfc, micro_platy$X9x9HAsfc, micro_catag$X9x9HAsfc, micro_WLP$X9x9HAsfc), method = "none", kw=TRUE)
## 1 = mylohyus, 2= platygonus, 3= catagonus, 4= wlp


#Kruskal Wallis with Dunn test evaluating Mylohyus microwear attributes through time

Asfc_Mylo.time <- dunn.test(list(mylo_ranch$Asfc, mylo_irv$Asfc, mylo_blanc$Asfc, mylo_hemph$Asfc), method = "none", kw=TRUE, label = TRUE)
epLsar_Mylo.time <- dunn.test(list(mylo_ranch$epLsar, mylo_irv$epLsar, mylo_blanc$epLsar, mylo_hemph$epLsar), method = "none", kw=TRUE, label = TRUE)
Tfv_Mylo.time <- dunn.test(list(mylo_ranch$Tfv, mylo_irv$Tfv, mylo_blanc$Tfv, mylo_hemph$Tfv), method = "none", kw=TRUE, label = TRUE)
X3x3HAsfc_Mylo.time <- dunn.test(list(mylo_ranch$X3x3HAsfc, mylo_irv$X3x3HAsfc, mylo_blanc$X3x3HAsfc, mylo_hemph$X3x3HAsfc), method = "none", kw=TRUE, label = TRUE)
X9x9HAsfc_Mylo.time <- dunn.test(list(mylo_ranch$X9x9HAsfc, mylo_irv$X9x9HAsfc, mylo_blanc$X9x9HAsfc, mylo_hemph$X9x9HAsfc), method = "none", kw=TRUE, label = TRUE)

Asfc_allpecc_TEST <- dunn.test(list(micro_mylo$Asfc, micro_platy$Asfc, micro_catag$Asfc, micro_WLP$Asfc), method = "none", kw=TRUE, label = TRUE)


# Microwear over time Platygonus with Dunn Test

Asfc_platy.time <- dunn.test(list(platy_ranch$Asfc, platy_irv$Asfc, platy_blanc$Asfc), method = "none", kw=TRUE, label = TRUE)
epLsar_platy.time <- dunn.test(list(platy_ranch$epLsar, platy_irv$epLsar, platy_blanc$epLsar), method = "none", kw=TRUE, label = TRUE)
Tfv_platy.time <- dunn.test(list(platy_ranch$Tfv, platy_irv$Tfv, platy_blanc$Tfv), method = "none", kw=TRUE, label = TRUE)
X3x3HAsfc_platy.time <- dunn.test(list(platy_ranch$X3x3HAsfc, platy_irv$X3x3HAsfc, platy_blanc$X3x3HAsfc), method = "none", kw=TRUE, label = TRUE)
X9x9HAsfc_platy.time <- dunn.test(list(platy_ranch$X9x9HAsfc, platy_irv$X9x9HAsfc, platy_blanc$X9x9HAsfc), method = "none", kw=TRUE, label = TRUE)







# # Calculating summary statistics per Genus %%% Method 1
# trial <- micro_WLP %>% select(Asfc, epLsar, Tfv)
# apply(trial, 2, mean)
# # gives the mean for each column as an array with headers labeled!

# Big_Function <- function (Genus) {
#   mean = apply(Genus, 2, mean)
#   median = apply(Genus, 2, median)
#   sd = apply(Genus, 2, sd)
#   max = apply(Genus, 2, max)
#   min = apply(Genus, 2, min)
#   range = max - min
#   Matrix <- matrix(NA, nrow= 6, ncol = 5)
#   colnames()
# }
#  # Put in matrix and return every variable
# trial_mylo <- filter(microwear, Genus == 'Mylohyus') %>% select(Asfc, epLsar, Tfv)
# 
# Big_Function(trial_mylo)
# 



################ %%%% ################ %%%% ###################
# Calculating summary statistics per Genus
SumStats <- microwear %>% group_by(Genus) %>% summarise(
  avg_Asfc = mean(Asfc),
  med_Asfc = median(Asfc),
  sd_Asfc = sd(Asfc),
  max_Asfc = max(Asfc),
  min_Asfc = min(Asfc),
  range_Asfc = max_Asfc - min_Asfc,
  
  avg_epLsar = mean(epLsar),
  med_epLsar = median(epLsar),
  sd_epLsar = sd(epLsar),
  max_epLsar = max(epLsar),
  min_epLsar = min(epLsar),
  range_epLsar = max_epLsar - min_epLsar,
  
  avg_Tfv = mean(Tfv),
  med_Tfv = median(Tfv),
  sd_Tfv = sd(Tfv),
  max_Tfv = max(Tfv),
  min_Tfv = min(Tfv),
  range_Tfv = max_Tfv - min_Tfv,
  
  avg_X3x3HAsfc = mean(X3x3HAsfc),
  med_X3x3HAsfc = median(X3x3HAsfc),
  sd_X3x3HAsfc = sd(X3x3HAsfc),
  max_X3x3HAsfc = max(X3x3HAsfc),
  min_X3x3HAsfc = min(X3x3HAsfc),
  range_X3x3HAsfc = max_X3x3HAsfc - min_X3x3HAsfc,
  
  avg_X9x9HAsfc = mean(X9x9HAsfc),
  med_X9x9HAsfc = median(X9x9HAsfc),
  sd_X9x9HAsfc = sd(X9x9HAsfc),
  max_X9x9HAsfc = max(X9x9HAsfc),
  min_X9x9HAsfc = min(X9x9HAsfc),
  range_X9x9HAsfc = max_X9x9HAsfc - min_X9x9HAsfc
  )

SumStats <- t(SumStats)
write.csv(SumStats, "Genus_SummaryStats.csv")



## Calculating P values between Genera for each microwear analysis 

## WRONG == SHOULD BE A DUNN TEST 
# Mann-Whitney for Genus$Asfc
Asfc_Matrix <- matrix(NA, nrow= 3, ncol = 4)
rownames(Asfc_Matrix) <- c("Mylohyus", "Platygonus", "Catagonus")
colnames(Asfc_Matrix) <- c("Mylohyus", "Platygonus", "Catagonus", "WLP")

Asfc_Matrix [1,2] <- MannWhitney(micro_mylo$Asfc, micro_platy$Asfc)
Asfc_Matrix [1,3] <- MannWhitney(micro_mylo$Asfc, micro_catag$Asfc)
Asfc_Matrix [1,4] <- MannWhitney(micro_mylo$Asfc, micro_WLP$Asfc)
Asfc_Matrix [2,3] <- MannWhitney(micro_platy$Asfc, micro_catag$Asfc)
Asfc_Matrix [2,4] <- MannWhitney(micro_platy$Asfc, micro_WLP$Asfc)
Asfc_Matrix [3,4] <- MannWhitney(micro_catag$Asfc, micro_WLP$Asfc)

table(microwear$Genus) #counts the sample size in per Genus
write.csv(Asfc_Matrix, 'Genus~Asfc_MannWhitney.csv')

# # gets mean of Asfc for each peccary and returns a vector
# x <- list(micro_mylo$Asfc, micro_platy$Asfc, micro_catag$Asfc, micro_WLP$Asfc)
# sapply(x, mean)


# Mann-Whitney for Genus$epLsar
epLsar_Matrix <- matrix(NA, nrow= 3, ncol = 4)
rownames(epLsar_Matrix) <- c("Mylohyus", "Platygonus", "Catagonus")
colnames(epLsar_Matrix) <- c("Mylohyus", "Platygonus", "Catagonus", "WLP")

epLsar_Matrix [1,2] <- MannWhitney(micro_mylo$epLsar, micro_platy$epLsar)
epLsar_Matrix [1,3] <- MannWhitney(micro_mylo$epLsar, micro_catag$epLsar)
epLsar_Matrix [1,4] <- MannWhitney(micro_mylo$epLsar, micro_WLP$epLsar)
epLsar_Matrix [2,3] <- MannWhitney(micro_platy$epLsar, micro_catag$epLsar)
epLsar_Matrix [2,4] <- MannWhitney(micro_platy$epLsar, micro_WLP$epLsar)
epLsar_Matrix [3,4] <- MannWhitney(micro_catag$epLsar, micro_WLP$epLsar)

write.csv(epLsar_Matrix, 'Genus~epLsar_MannWhitney.csv')

# Mann-Whitney for Genus$Tfv
Tfv_Matrix <- matrix(NA, nrow= 3, ncol = 4)
rownames(Tfv_Matrix) <- c("Mylohyus", "Platygonus", "Catagonus")
colnames(Tfv_Matrix) <- c("Mylohyus", "Platygonus", "Catagonus", "WLP")

Tfv_Matrix [1,2] <- MannWhitney(micro_mylo$Tfv, micro_platy$Tfv)
Tfv_Matrix [1,3] <- MannWhitney(micro_mylo$Tfv, micro_catag$Tfv)
Tfv_Matrix [1,4] <- MannWhitney(micro_mylo$Tfv, micro_WLP$Tfv)
Tfv_Matrix [2,3] <- MannWhitney(micro_platy$Tfv, micro_catag$Tfv)
Tfv_Matrix [2,4] <- MannWhitney(micro_platy$Tfv, micro_WLP$Tfv)
Tfv_Matrix [3,4] <- MannWhitney(micro_catag$Tfv, micro_WLP$Tfv)

write.csv(Tfv_Matrix, 'Genus~Tfv_MannWhitney.csv')

# Mann-Whitney for Genus$HAsfc3x3
X3x3HAsfc_Matrix <- matrix(NA, nrow= 3, ncol = 4)
rownames(X3x3HAsfc_Matrix) <- c("Mylohyus", "Platygonus", "Catagonus")
colnames(X3x3HAsfc_Matrix) <- c("Mylohyus", "Platygonus", "Catagonus", "WLP")

X3x3HAsfc_Matrix [1,2] <- MannWhitney(micro_mylo$X3x3HAsfc, micro_platy$X3x3HAsfc)
X3x3HAsfc_Matrix [1,3] <- MannWhitney(micro_mylo$X3x3HAsfc, micro_catag$X3x3HAsfc)
X3x3HAsfc_Matrix [1,4] <- MannWhitney(micro_mylo$X3x3HAsfc, micro_WLP$X3x3HAsfc)
X3x3HAsfc_Matrix [2,3] <- MannWhitney(micro_platy$X3x3HAsfc, micro_catag$X3x3HAsfc)
X3x3HAsfc_Matrix [2,4] <- MannWhitney(micro_platy$X3x3HAsfc, micro_WLP$X3x3HAsfc)
X3x3HAsfc_Matrix [3,4] <- MannWhitney(micro_catag$X3x3HAsfc, micro_WLP$X3x3HAsfc)

write.csv(X3x3HAsfc_Matrix, 'Genus~X3x3HAsfc_MannWhitney.csv')

# Mann-Whitney for Genus$HAsfc9x9
X9x9HAsfc_Matrix <- matrix(NA, nrow= 3, ncol = 4)
rownames(X9x9HAsfc_Matrix) <- c("Mylohyus", "Platygonus", "Catagonus")
colnames(X9x9HAsfc_Matrix) <- c("Mylohyus", "Platygonus", "Catagonus", "WLP")

X9x9HAsfc_Matrix [1,2] <- MannWhitney(micro_mylo$X9x9HAsfc, micro_platy$X9x9HAsfc)
X9x9HAsfc_Matrix [1,3] <- MannWhitney(micro_mylo$X9x9HAsfc, micro_catag$X9x9HAsfc)
X9x9HAsfc_Matrix [1,4] <- MannWhitney(micro_mylo$X9x9HAsfc, micro_WLP$X9x9HAsfc)
X9x9HAsfc_Matrix [2,3] <- MannWhitney(micro_platy$X9x9HAsfc, micro_catag$X9x9HAsfc)
X9x9HAsfc_Matrix [2,4] <- MannWhitney(micro_platy$X9x9HAsfc, micro_WLP$X9x9HAsfc)
X9x9HAsfc_Matrix [3,4] <- MannWhitney(micro_catag$X9x9HAsfc, micro_WLP$X9x9HAsfc)

write.csv(X9x9HAsfc_Matrix, 'Genus~X9x9HAsfc_MannWhitney.csv')






################ %%%% ################ %%%% ###################

## Microwear through time for each genus
SumStats_genus.time <- microwear %>% group_by(NALMA, Genus) %>% summarise(
  avg_Asfc = mean(Asfc),
  med_Asfc = median(Asfc),
  sd_Asfc = sd(Asfc),
  max_Asfc = max(Asfc),
  min_Asfc = min(Asfc),
  range_Asfc = max_Asfc - min_Asfc,
  
  avg_epLsar = mean(epLsar),
  med_epLsar = median(epLsar),
  sd_epLsar = sd(epLsar),
  max_epLsar = max(epLsar),
  min_epLsar = min(epLsar),
  range_epLsar = max_epLsar - min_epLsar,
  
  avg_Tfv = mean(Tfv),
  med_Tfv = median(Tfv),
  sd_Tfv = sd(Tfv),
  max_Tfv = max(Tfv),
  min_Tfv = min(Tfv),
  range_Tfv = max_Tfv - min_Tfv,
  
  avg_X3x3HAsfc = mean(X3x3HAsfc),
  med_X3x3HAsfc = median(X3x3HAsfc),
  sd_X3x3HAsfc = sd(X3x3HAsfc),
  max_X3x3HAsfc = max(X3x3HAsfc),
  min_X3x3HAsfc = min(X3x3HAsfc),
  range_X3x3HAsfc = max_X3x3HAsfc - min_X3x3HAsfc,
  
  avg_X9x9HAsfc = mean(X9x9HAsfc),
  med_X9x9HAsfc = median(X9x9HAsfc),
  sd_X9x9HAsfc = sd(X9x9HAsfc),
  max_X9x9HAsfc = max(X9x9HAsfc),
  min_X9x9HAsfc = min(X9x9HAsfc),
  range_X9x9HAsfc = max_X9x9HAsfc - min_X9x9HAsfc
)

SumStats_genus.time <- t(SumStats_genus.time) #transpose data
write.csv(SumStats_genus.time, "Genus-Time_SummaryStats.csv")




######## %%%%% ######
## MYLOHYUS -- Calculating P values ~ time period | Microwear analysis 

## WRONG THIS SHOULD BE A DUNN TEST

# Mann-Whitney for Genus$Asfc
Mylo_Matrix <- matrix(NA, nrow= 15, ncol = 4)
rownames(Mylo_Matrix) <- c("Hemphilian", "Irvingtonian", "Blancan", "Hemphilian", "Irvingtonian", "Blancan","Hemphilian", "Irvingtonian", "Blancan","Hemphilian", "Irvingtonian", "Blancan","Hemphilian", "Irvingtonian", "Blancan")
colnames(Mylo_Matrix) <- c("DMTA feature", "Irvingtonian", "Blancan", "Rancholabrean")

# Mylohyus Asfc
Mylo_Matrix [1,1] <- 'Asfc'
Mylo_Matrix [1,2] <- MannWhitney(mylo_hemph$Asfc, mylo_irv$Asfc)
Mylo_Matrix [1,3] <- MannWhitney(mylo_hemph$Asfc, mylo_blanc$Asfc)
Mylo_Matrix [1,4] <- MannWhitney(mylo_hemph$Asfc, mylo_ranch$Asfc)
Mylo_Matrix [2,3] <- MannWhitney(mylo_irv$Asfc, mylo_blanc$Asfc)
Mylo_Matrix [2,4] <- MannWhitney(mylo_irv$Asfc, mylo_ranch$Asfc)
Mylo_Matrix [3,4] <- MannWhitney(mylo_blanc$Asfc, mylo_ranch$Asfc)

#Mylohyus epLsar
Mylo_Matrix [4,1] <- 'epLsar'
Mylo_Matrix [4,2] <- MannWhitney(mylo_hemph$epLsar, mylo_irv$epLsar)
Mylo_Matrix [4,3] <- MannWhitney(mylo_hemph$epLsar, mylo_blanc$epLsar)
Mylo_Matrix [4,4] <- MannWhitney(mylo_hemph$epLsar, mylo_ranch$epLsar)
Mylo_Matrix [5,3] <- MannWhitney(mylo_irv$epLsar, mylo_blanc$epLsar)
Mylo_Matrix [5,4] <- MannWhitney(mylo_irv$epLsar, mylo_ranch$epLsar)
Mylo_Matrix [6,4] <- MannWhitney(mylo_blanc$epLsar, mylo_ranch$epLsar)

#Mylohyus Tfv
Mylo_Matrix [7,1] <- 'Tfv'
Mylo_Matrix [7,2] <- MannWhitney(mylo_hemph$Tfv, mylo_irv$Tfv)
Mylo_Matrix [7,3] <- MannWhitney(mylo_hemph$Tfv, mylo_blanc$Tfv)
Mylo_Matrix [7,4] <- MannWhitney(mylo_hemph$Tfv, mylo_ranch$Tfv)
Mylo_Matrix [8,3] <- MannWhitney(mylo_irv$Tfv, mylo_blanc$Tfv)
Mylo_Matrix [8,4] <- MannWhitney(mylo_irv$Tfv, mylo_ranch$Tfv)
Mylo_Matrix [9,4] <- MannWhitney(mylo_blanc$Tfv, mylo_ranch$Tfv)

#Mylohyus HAsfc(3x3)
Mylo_Matrix [10,1] <- 'HAsfc3x3'
Mylo_Matrix [10,2] <- MannWhitney(mylo_hemph$X3x3HAsfc, mylo_irv$X3x3HAsfc)
Mylo_Matrix [10,3] <- MannWhitney(mylo_hemph$X3x3HAsfc, mylo_blanc$X3x3HAsfc)
Mylo_Matrix [10,4] <- MannWhitney(mylo_hemph$X3x3HAsfc, mylo_ranch$X3x3HAsfc)
Mylo_Matrix [11,3] <- MannWhitney(mylo_irv$X3x3HAsfc, mylo_blanc$X3x3HAsfc)
Mylo_Matrix [11,4] <- MannWhitney(mylo_irv$X3x3HAsfc, mylo_ranch$X3x3HAsfc)
Mylo_Matrix [12,4] <- MannWhitney(mylo_blanc$X3x3HAsfc, mylo_ranch$X3x3HAsfc)

#Mylohyus HAsfc(9x9)
Mylo_Matrix [13,1] <- 'HAsfc9x9'
Mylo_Matrix [13,2] <- MannWhitney(mylo_hemph$X9x9HAsfc, mylo_irv$X9x9HAsfc)
Mylo_Matrix [13,3] <- MannWhitney(mylo_hemph$X9x9HAsfc, mylo_blanc$X9x9HAsfc)
Mylo_Matrix [13,4] <- MannWhitney(mylo_hemph$X9x9HAsfc, mylo_ranch$X9x9HAsfc)
Mylo_Matrix [14,3] <- MannWhitney(mylo_irv$X9x9HAsfc, mylo_blanc$X9x9HAsfc)
Mylo_Matrix [14,4] <- MannWhitney(mylo_irv$X9x9HAsfc, mylo_ranch$X9x9HAsfc)
Mylo_Matrix [15,4] <- MannWhitney(mylo_blanc$X9x9HAsfc, mylo_ranch$X9x9HAsfc)

# Getting sample sizes of Mylohyus per NALMA
nrow(mylo_hemph)
nrow(mylo_irv)
nrow(mylo_blanc)
nrow(mylo_ranch)

write.csv(Mylo_Matrix, 'Mylohyus-NALMA_Pvalues.csv')



######## %%%%% ######
## PLATYGONUS -- Calculating P values ~ time period | Microwear analysis 

## WRONG -> THIS SHOULD BE A DUNN TEST

# Mann-Whitney for Genus$Asfc
Platy_Matrix <- matrix(NA, nrow= 10, ncol = 3)
rownames(Platy_Matrix) <- c("Irvingtonian", "Blancan", "Irvingtonian", "Blancan", "Irvingtonian", "Blancan", "Irvingtonian", "Blancan", "Irvingtonian", "Blancan")
colnames(Platy_Matrix) <- c("DMTA Feature", "Blancan", "Rancholabrean")


#Platygonus Asfc
Platy_Matrix [1,1] <- 'Asfc'
Platy_Matrix [1,2] <- MannWhitney(platy_irv$Asfc, platy_blanc$Asfc)
Platy_Matrix [1,3] <- MannWhitney(platy_irv$Asfc, platy_ranch$Asfc)
Platy_Matrix [2,3] <- MannWhitney(platy_blanc$Asfc, platy_ranch$Asfc)

#Platygonus epLsar
Platy_Matrix [3,1] <- 'epLsar'
Platy_Matrix [3,2] <- MannWhitney(platy_irv$epLsar, platy_blanc$epLsar)
Platy_Matrix [3,3] <- MannWhitney(platy_irv$epLsar, platy_ranch$epLsar)
Platy_Matrix [4,3] <- MannWhitney(platy_blanc$epLsar, platy_ranch$epLsar)

#Platygonus Tfv
Platy_Matrix [5,1] <- 'Tfv'
Platy_Matrix [5,2] <- MannWhitney(platy_irv$Tfv, platy_blanc$Tfv)
Platy_Matrix [5,3] <- MannWhitney(platy_irv$Tfv, platy_ranch$Tfv)
Platy_Matrix [6,3] <- MannWhitney(platy_blanc$Tfv, platy_ranch$Tfv)

#Platygonus HAsfc(3x3)
Platy_Matrix [7,1] <- 'HAsfc3x3'
Platy_Matrix [7,2] <- MannWhitney(platy_irv$X3x3HAsfc, platy_blanc$X3x3HAsfc)
Platy_Matrix [7,3] <- MannWhitney(platy_irv$X3x3HAsfc, platy_ranch$X3x3HAsfc)
Platy_Matrix [8,3] <- MannWhitney(platy_blanc$X3x3HAsfc, platy_ranch$X3x3HAsfc)

#Platygonus HAsfc(9x9)
Platy_Matrix [9,1] <- 'HAsfc9x9'
Platy_Matrix [9,2] <- MannWhitney(platy_irv$X9x9HAsfc, platy_blanc$X9x9HAsfc)
Platy_Matrix [9,3] <- MannWhitney(platy_irv$X9x9HAsfc, platy_ranch$X9x9HAsfc)
Platy_Matrix [10,3] <- MannWhitney(platy_blanc$X9x9HAsfc, platy_ranch$X9x9HAsfc)

# Getting sample sizes of Platygonus per NALMA
nrow(platy_irv)
nrow(platy_blanc)
nrow(platy_ranch)

write.csv(Platy_Matrix, 'Platygonus-NALMA_Pvalues.csv')



######## %%%%% ######
## MYLOHYUS vs PLATYGONUS -- Calculating P values between Genera ~ time period | Microwear analysis 

Matrix_platy.mylo <- matrix(NA, nrow=5, ncol=3)
colnames(Matrix_platy.mylo) <- c('Irvingtonian', 'Blancan', 'Rancholabrean')
rownames(Matrix_platy.mylo) <- c('Asfc', 'epLsar', 'Tfv', 'HAsfc3x3', 'HAsfc9x9')

#Mylohyus vs Platygonus __ Irvingtonian
Matrix_platy.mylo [1,1] <- MannWhitney(mylo_irv$Asfc, platy_irv$Asfc)
Matrix_platy.mylo [2,1] <- MannWhitney(mylo_irv$epLsar, platy_irv$epLsar)
Matrix_platy.mylo [3,1] <- MannWhitney(mylo_irv$Tfv, platy_irv$Tfv)
Matrix_platy.mylo [4,1] <- MannWhitney(mylo_irv$X3x3HAsfc, platy_irv$X3x3HAsfc)
Matrix_platy.mylo [5,1] <- MannWhitney(mylo_irv$X9x9HAsfc, platy_irv$X9x9HAsfc)

#Mylohyus vs Platygonus __ Blancan
Matrix_platy.mylo [1,2] <- MannWhitney(mylo_blanc$Asfc, platy_blanc$Asfc)
Matrix_platy.mylo [2,2] <- MannWhitney(mylo_blanc$epLsar, platy_blanc$epLsar)
Matrix_platy.mylo [3,2] <- MannWhitney(mylo_blanc$Tfv, platy_blanc$Tfv)
Matrix_platy.mylo [4,2] <- MannWhitney(mylo_blanc$X3x3HAsfc, platy_blanc$X3x3HAsfc)
Matrix_platy.mylo [5,2] <- MannWhitney(mylo_blanc$X9x9HAsfc, platy_blanc$X9x9HAsfc)

#Mylohyus vs Platygonus __ Rancholabrean
Matrix_platy.mylo [1,3] <- MannWhitney(mylo_ranch$Asfc, platy_ranch$Asfc)
Matrix_platy.mylo [2,3] <- MannWhitney(mylo_ranch$epLsar, platy_ranch$epLsar)
Matrix_platy.mylo [3,3] <- MannWhitney(mylo_ranch$Tfv, platy_ranch$Tfv)
Matrix_platy.mylo [4,3] <- MannWhitney(mylo_ranch$X3x3HAsfc, platy_ranch$X3x3HAsfc)
Matrix_platy.mylo [5,3] <- MannWhitney(mylo_ranch$X9x9HAsfc, platy_ranch$X9x9HAsfc)

write.csv(Matrix_platy.mylo, 'NALMA_platy-mylo.csv')




######## %%%%% ######
## MYLOHYUS vs Catagonus -- Calculating P values between Genera ~ time period | Microwear analysis 

Matrix_catag.mylo <- matrix(NA, nrow=5, ncol=1)
colnames(Matrix_catag.mylo) <- c('Hemphillian')
rownames(Matrix_catag.mylo) <- c('Asfc', 'epLsar', 'Tfv', 'HAsfc3x3', 'HAsfc9x9')

#Mylohyus vs Platygonus __ Irvingtonian
Matrix_catag.mylo [1,1] <- MannWhitney(mylo_hemph$Asfc, catag_hemph$Asfc)
Matrix_catag.mylo [2,1] <- MannWhitney(mylo_hemph$epLsar, catag_hemph$epLsar)
Matrix_catag.mylo [3,1] <- MannWhitney(mylo_hemph$Tfv, catag_hemph$Tfv)
Matrix_catag.mylo [4,1] <- MannWhitney(mylo_hemph$X3x3HAsfc, catag_hemph$X3x3HAsfc)
Matrix_catag.mylo [5,1] <- MannWhitney(mylo_hemph$X9x9HAsfc, catag_hemph$X9x9HAsfc)

write.csv(Matrix_catag.mylo, 'NALMA_catag-mylo.csv')


## WLP vs each extinct species
micro_wlp.extinct <- matrix(NA, nrow=5, ncol=3)
rownames(micro_wlp.extinct) <- c('Asfc', 'epLsar', 'Tfv', 'HAsfc3x3', 'HAsfc9x9')
colnames(micro_wlp.extinct) <- c('Catagonus','Mylohyus','Platygonus')

micro_wlp.extinct [1,1] <- MannWhitney(micro_WLP$Asfc, micro_catag$Asfc)
micro_wlp.extinct [1,2] <- MannWhitney(micro_WLP$Asfc, micro_mylo$Asfc)
micro_wlp.extinct [1,3] <- MannWhitney(micro_WLP$Asfc, micro_platy$Asfc)

micro_wlp.extinct [2,1] <- MannWhitney(micro_WLP$epLsar, micro_catag$epLsar)
micro_wlp.extinct [2,2] <- MannWhitney(micro_WLP$epLsar, micro_mylo$epLsar)
micro_wlp.extinct [2,3] <- MannWhitney(micro_WLP$epLsar, micro_platy$epLsar)

micro_wlp.extinct [3,1] <- MannWhitney(micro_WLP$Tfv, micro_catag$Tfv)
micro_wlp.extinct [3,2] <- MannWhitney(micro_WLP$Tfv, micro_mylo$Tfv)
micro_wlp.extinct [3,3] <- MannWhitney(micro_WLP$Tfv, micro_platy$Tfv)

micro_wlp.extinct [4,1] <- MannWhitney(micro_WLP$X3x3HAsfc, micro_catag$X3x3HAsfc)
micro_wlp.extinct [4,2] <- MannWhitney(micro_WLP$X3x3HAsfc, micro_mylo$X3x3HAsfc)
micro_wlp.extinct [4,3] <- MannWhitney(micro_WLP$X3x3HAsfc, micro_platy$X3x3HAsfc)

micro_wlp.extinct [5,1] <- MannWhitney(micro_WLP$X9x9HAsfc, micro_catag$X9x9HAsfc)
micro_wlp.extinct [5,2] <- MannWhitney(micro_WLP$X9x9HAsfc, micro_mylo$X9x9HAsfc)
micro_wlp.extinct [5,3] <- MannWhitney(micro_WLP$X9x9HAsfc, micro_platy$X9x9HAsfc)

write.csv(micro_wlp.extinct, 'micro_wlp-extinct_Pvalues.csv')





## ---------------- Microwear vs Carbon ------------------ ##

qplot(Delta13C, Asfc, data=  micro_isotopes, color = Genus, main = "Extinct Peccaries: Asfc vs Delta13C")
qplot(Delta13C, epLsar, data=  micro_isotopes, color = Genus, main = 'Extinct Peccaries: epLsar vs Delta13C')






