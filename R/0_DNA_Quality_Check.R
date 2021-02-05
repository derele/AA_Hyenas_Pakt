library("ggplot2")
library(gridExtra)
library(grid)
library(tidyverse)
library("frequency")
library(httr)

Hyena_O <-read.csv("~/GitProjects/AA_Hyenas_Pakt/Data/Hyena_DNA_Concentration.csv") ##Original concentrations
Hyena_O$Sample<- gsub(pattern = " ", replacement = "", x= Hyena_O$Sample)

Hyena_O[duplicated(Hyena_O$Sample), ]
Hyena_O %>% distinct(Sample, .keep_all= T)-> Hyena_O

Hyena_F <-read.csv("~/GitProjects/AA_Hyenas_Pakt/Data/Hyena_DNA_Concentration_Speedvac.csv") ##After speedvac concentrations
Hyena_F$Sample<- gsub(pattern = " ", replacement = "", x= Hyena_F$Sample)
Hyena_F[duplicated(Hyena_F$Sample), ]
Hyena_F %>% distinct(Sample, .keep_all= T)-> Hyena_F

###Check of Original samples 
#How many samples have "good quality" DNA?

##Strict conditions
table(Hyena_O$DNA_conc_O >= 50 & Hyena_O$P260_280_O >= 1.8 & Hyena_O$P260_230_O >= 2.0)

##Relaxed conditions
table(Hyena_O$DNA_conc_O >= 50 & Hyena_O$P260_280_O >= 1.8 & Hyena_O$P260_230_O >= 1.5)

##Just by amount of DNA
table(Hyena_O$DNA_conc_O >= 50)

##Potentialy to be concentrated
table(Hyena_O$DNA_conc_O >= 30 & Hyena_O$DNA_conc_O < 50)


##Actions to take for DNAs 
Hyena_O <- mutate(Hyena_O, DNA_Quality = ifelse(DNA_conc_O >= 50.0, "Optimal", "Low"))
Hyena_O <- mutate(Hyena_O, Purity_1 = ifelse(P260_280_O >= 1.80, "Optimal", "Low"))
Hyena_O <- mutate(Hyena_O, Purity_2 = ifelse(P260_230_O >= 2.00, "Optimal", "Low"))

Hyena_O %>% mutate(Action = case_when(DNA_Quality == "Low" & Purity_1 == "Low" & Purity_2 == "Low"  ~ "Reextraction",
                                      DNA_Quality == "Low" & Purity_1 == "Optimal"  ~ "Concentrate",
                                      DNA_Quality == "Optimal" & Purity_1 == "Optimal" ~ "Metabarcoding")) -> Hyena_O

##Now merged Original concentration with after speedvac concentration

Hyena<- plyr::join(Hyena_O, Hyena_F, by= "Sample")

a <- ggplot(Hyena, aes(x=P260_280_F, y=DNA_conc_F, color= P260_230_F)) + 
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=45, linetype="dashed", color = "red")+
  scale_color_gradient(low="red", high="green", name= "Purity\n(260/230 ratio)")+
  xlab("Purity (260/280 ratio)")+
  ylab("DNA concentration ng/µL")+
  labs(tag = "A)")

b <-ggplot(Hyena, aes(x=P260_230_F, y=DNA_conc_F, color= P260_280_F)) + 
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=45, linetype="dashed", color = "red")+
  scale_color_gradient(low="red", high="green", name= "Purity\n(260/280 ratio)")+
  xlab("Purity (260/230 ratio)")+
  ylab("DNA concentration ng/µL")+
  labs(tag = "B)")

c <- ggplot(Hyena, aes(x=P260_230_F, y=P260_280_F, color=DNA_conc_F)) + 
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=1.8, linetype="dashed", color = "red")+
  geom_vline(xintercept =2.0, linetype="dashed", color = "red")+
  scale_color_gradient(low="red", high="green", name="DNA\nconcentration\nng/µL")+
  xlab("Purity (260/230 ratio)")+
  ylab("Purity (260/280 ratio)")+
  labs(tag = "C)")

pdf("~/GitProjects/AA_Hyenas_Pakt/Figures/Hyena_DNA_Concentration.pdf", width=10,height=10)
gridExtra::grid.arrange(a, b, c, ncol=3)
dev.off()

###Selection of samples
Hyena %>% mutate(Final_selection = 
                   case_when(DNA_conc_F >= 45.0 ~ "Selected")) -> Hyena

table(Hyena$Final_selection == "Selected")

#Hyena %>% 
#  add_row(Sample= 203:240) -> Hyena
#Hyena[203:240,]<- "Negative" ##Adding negative controls 
#Chip <- as.data.frame(split(Hyena$Sample, sample(5))) 
#colnames(Chip)<- c("Chip_1", "Chip_2", "Chip_3", "Chip_4", "Chip_5")

#for (i in 1:ncol(Chip)) {
#  Chip[,i]<- sample(Chip[,i])
#}

#write.csv(Chip, "~/GitProjects/AA_Hyenas_Pakt/Data/Sample_distribution.csv")

##Add index information 
Index_1 <-read.csv("~/AA_Hyena_Pakt/Index_Pool_1.csv") ##Pool 1 Dec 2019
Index_2 <-read.csv("~/AA_Hyena_Pakt/Index_Pool_2.csv") ##Pool 2 Jan 2020
Index<- merge(Index_1, Index_2, by= "Sample")

rm(Hyena_F, Hyena_O, Index_1, Index_2)

###Final data set with all information including negatives 
Hyena<- full_join(Hyena, Index, by= "Sample")
