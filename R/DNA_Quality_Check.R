library("ggplot2")
library(gridExtra)
library(grid)
library(dplyr)
library("frequency")
library(httr)

Hyena <-read.csv("~/GitProject/AA_Hyenas_Pakt/Data/Hyena_DNA_Concentration.csv")

a <- ggplot(Hyena, aes(x=X260_280, y=DNA.conc, color= X260_230)) + 
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=50, linetype="dashed", color = "red")+
  scale_color_gradient(low="red", high="green", name= "Purity\n(260/230 ratio)")+
  xlab("Purity (260/280 ratio)")+
  ylab("DNA concentration ng/µL")
  
b <-ggplot(Hyena, aes(x=X260_230, y=DNA.conc, color= X260_280)) + 
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=50, linetype="dashed", color = "red")+
  scale_color_gradient(low="red", high="green", name= "Purity\n(260/280 ratio)")+
  xlab("Purity (260/230 ratio)")+
  ylab("DNA concentration ng/µL")

c <- ggplot(Hyena, aes(x=X260_230, y=X260_280, color=DNA.conc)) + 
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=1.8, linetype="dashed", color = "red")+
  geom_vline(xintercept =2.0, linetype="dashed", color = "red")+
  scale_color_gradient(low="red", high="green", name="DNA\nconcentration\nng/µL")+
  xlab("Purity (260/230 ratio)")+
  ylab("Purity (260/280 ratio)")

grid.arrange(a, b, c, ncol=3)

#How many samples have "good quality" DNA?

##Strict conditions
table(Hyena$DNA.conc >= 50 & Hyena$X260_280 >= 1.8 & Hyena$X260_230 >= 2.0)

##Relaxed conditions
table(Hyena$DNA.conc >= 50 & Hyena$X260_280 >= 1.8 & Hyena$X260_230 >= 1.5)

##Just by amount of DNA
table(Hyena$DNA.conc >= 50)

##Potentialy to be concentrated
table(Hyena$DNA.conc >= 30 & Hyena$DNA.conc < 50)
