### FOR COMPOSITION WE USE proportional transformation

library(microbiome)

############ Preparing datasets for ordination
prune_prune <- function (ps, samples=0, taxa=10) {
    s <- prune_samples(sample_sums(ps) > samples, ps)
    S <- prune_samples(!grepl("Negative", sample_names(s)), s)
    prune_taxa(taxa_sums(S) > taxa, S)
}
pseq.compositional <- microbiome::transform(pseq, "compositional")

## trying the species glomed dataset
cBac <- microbiome::transform(
                        prune_prune(PBac, 20, 10),
                        "compositional")

cEuk <- microbiome::transform(
                        prune_prune(PEuk, 20, 10),
                        "compositional")

    

library(vegan)

######### Prepare environmental data
EDatrawBac <- sample_data(cBac)
class(EDatrawBac) <- "data.frame"

EDatrawEuk <- sample_data(cEuk)
class(EDatrawEuk) <- "data.frame"

immuno <- c("cortisol", "neopterin", "lysozyme", "IgG", "IgA", "mucin")

para <- c("Ancylostoma_egg_load", "Cystoisospora_oocyst_load",
          "Spirometra_egg_load", "Dipylidium_egg_load",
          "Trichuris_egg_load", "Taeniidae_egg_load",
          "Spirurida_egg_load", "Mesocestoides_egg_load")

lhist <- c("season", "sampling_month",
           "prey_level_imputed", # "prey_level",
           "CSocialRank",
           "age_sampling", "sex", "clan",
           "hyena_ID")

EDatBac <- EDatrawBac[, c(immuno, para, lhist)]
EDatEuk <- EDatrawEuk[, c(immuno, para, lhist)]


### Bacterial nMDS ##############################################
BacData <- otu_table(cBac)
colnames(BacData) <- as.vector(tax_table(cBac)[, "species"])

mdsBac <- vegan::metaMDS(BacData, try=350, trymax=350, k=3)

BacEFit <- envfit(mdsBac, EDatBac, na.rm=TRUE)
BacEFit

### Eukaryote nMDS ###########################################

EukData <- otu_table(cEuk)
colnames(EukData) <- as.vector(tax_table(cEuk)[, "genus"])

mdsEuk <- vegan::metaMDS(EukData, try=350, trymax=350, k=3,
                         noshare=0.3)

EukEFit <- envfit(mdsEuk, EDatEuk, na.rm=TRUE)
EukEFit


### exporting the data for plotting with ggplot

data.scoresEuk <-  as.data.frame(scores(mdsEuk))
data.scoresEuk <- cbind(data.scoresEuk, EDatEuk)

EnCoordEuk <-
    as.data.frame(
        rbind(scores(EukEFit, "vectors") *
              ordiArrowMul(EukEFit),
              scores(EukEFit, "factors") *
              ordiArrowMul(EukEFit))
    )

EnCoordEuk$Cat <- ifelse(rownames(EnCoordEuk)%in%immuno, "Immune",
                  ifelse(rownames(EnCoordEuk)%in%para, "Parasite", "other"))


data.scoresBac <-  as.data.frame(scores(mdsBac))
data.scoresBac <- cbind(data.scoresBac, EDatBac)

EnCoordBac <-
    as.data.frame(
        rbind(scores(BacEFit, "vectors") * ordiArrowMul(BacEFit),
              scores(BacEFit, "factors") * ordiArrowMul(BacEFit))
    )

EnCoordBac$Cat <- ifelse(rownames(EnCoordBac)%in%immuno, "Immune",
                  ifelse(rownames(EnCoordBac)%in%para, "Parasite", "other"))


### the basic plots
ggBac <-  ggplot(data = data.scoresBac, aes(x = NMDS1, y = NMDS2)) +
    geom_point(data = data.scoresBac, aes(colour = age_sampling), size = 3, alpha = 0.5) +
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey30"),
          axis.ticks = element_blank(), axis.text = element_blank(),
          legend.key = element_blank(),
          legend.title = element_text(size = 10, face = "bold", colour = "grey30"),
          legend.text = element_text(size = 9, colour = "grey30"))


pdf("Figures/Age_Bac_ordi_EFit.pdf")
plot(ggBac)
dev.off()

ggBacIm <- ggBac +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 data = subset(EnCoordBac, Cat%in%"Immune"), size =1, alpha = 0.5) +
    geom_point(data = subset(EnCoordBac, Cat%in%"Immune"), aes(x = NMDS1, y = NMDS2),
               shape = "diamond", size = 4, alpha = 0.6) +    
    geom_text(data = subset(EnCoordBac, Cat%in%"Immune"), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(subset(EnCoordBac, Cat%in%"Immune")), colour = "red", fontface = "bold")

pdf("Figures/Age_Bac_ordi_EFit_Immune.pdf")
ggBacIm
dev.off()

ggBacPara <- ggBac +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 data = subset(EnCoordBac, Cat%in%"Parasite"), size =1, alpha = 0.5) +
    geom_point(data = subset(EnCoordBac, Cat%in%"Parasite"), aes(x = NMDS1, y = NMDS2),
               shape = "diamond", size = 4, alpha = 0.6) +    
    geom_text(data = subset(EnCoordBac, Cat%in%"Parasite"), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(subset(EnCoordBac, Cat%in%"Parasite")), colour = "red", fontface = "bold") 

pdf("Figures/Age_Bac_ordi_EFit_Parasite.pdf")
ggBacPara
dev.off()

ggBacOther <- ggBac +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 data = subset(EnCoordBac, Cat%in%"other"), size =1, alpha = 0.5) +
    geom_point(data = subset(EnCoordBac, Cat%in%"other"), aes(x = NMDS1, y = NMDS2),
               shape = "diamond", size = 4, alpha = 0.6) +    
    geom_text(data = subset(EnCoordBac, Cat%in%"other"), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(subset(EnCoordBac, Cat%in%"other")), colour = "red", fontface = "bold")

pdf("Figures/Age_Bac_ordi_EFit_other.pdf")
ggBacOther
dev.off()


BacSig <- names(BacEFit$vectors$pvals)[BacEFit$vectors$pvals<.05]


ggBacSig <- ggBac +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 data = EnCoordBac[BacSig,],
                 size =1, alpha = 0.5) +
    geom_point(data = EnCoordBac[BacSig,],
               aes(x = NMDS1, y = NMDS2),
               shape = "diamond", size = 4, alpha = 0.6) +    
    geom_text(data = EnCoordBac[BacSig,], 
              aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(EnCoordBac[BacSig,]), 
              colour = "red", fontface = "bold")

pdf("Figures/Age_Bac_ordi_EFit_sig.pdf")
ggBacSig
dev.off()



ggEuk <-  ggplot(data = data.scoresEuk, aes(x = NMDS2, y = NMDS3)) +
    geom_point(data = data.scoresEuk, aes(colour = age_sampling, shape=season),
               size = 3, alpha = 0.5) +
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey30"),
          axis.ticks = element_blank(), axis.text = element_blank(),
          legend.key = element_blank(),
          legend.title = element_text(size = 10, face = "bold", colour = "grey30"),
          legend.text = element_text(size = 9, colour = "grey30"))


pdf("Figures/Age_Euk_ordi_EFit.pdf")
plot(ggEuk)
dev.off()



pdf("Figures/Age_Euk_ordi_EFit_Immune.pdf")
ggEuk +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 data = subset(EnCoordEuk, Cat%in%"Immune"), size =1, alpha = 0.5) +
    geom_point(data = subset(EnCoordEuk, Cat%in%"Immune"), aes(x = NMDS1, y = NMDS2),
               shape = "diamond", size = 4, alpha = 0.6) +    
    geom_text(data = subset(EnCoordEuk, Cat%in%"Immune"), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(subset(EnCoordEuk, Cat%in%"Immune")), colour = "red", fontface = "bold") 
dev.off()


pdf("Figures/Age_Euk_ordi_EFit_Parasite.pdf")
ggEuk +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 data = subset(EnCoordEuk, Cat%in%"Parasite"), size =1, alpha = 0.5) +
    geom_point(data = subset(EnCoordEuk, Cat%in%"Parasite"), aes(x = NMDS1, y = NMDS2),
               shape = "diamond", size = 4, alpha = 0.6) +    
    geom_text(data = subset(EnCoordEuk, Cat%in%"Parasite"), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(subset(EnCoordEuk, Cat%in%"Parasite")), colour = "red", fontface = "bold") 
dev.off()


pdf("Figures/Age_Euk_ordi_EFit_other.pdf")
ggEuk +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 data = subset(EnCoordEuk, Cat%in%"other"), size =1, alpha = 0.5) +
    geom_point(data = subset(EnCoordEuk, Cat%in%"other"), aes(x = NMDS1, y = NMDS2),
               shape = "diamond", size = 4, alpha = 0.6) +    
    geom_text(data = subset(EnCoordEuk, Cat%in%"other"), aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(subset(EnCoordEuk, Cat%in%"other")), colour = "red", fontface = "bold") 
dev.off()



EukSig <- names(EukEFit$vectors$pvals)[EukEFit$vectors$pvals<.05]


ggEukSig <- ggEuk +
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 data = EnCoordEuk[EukSig,],
                 size =1, alpha = 0.5) +
    geom_point(data = EnCoordEuk[EukSig,],
               aes(x = NMDS1, y = NMDS2),
               shape = "diamond", size = 4, alpha = 0.6) +    
    geom_text(data = EnCoordEuk[EukSig,], 
              aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(EnCoordEuk[EukSig,]), 
              colour = "red", fontface = "bold")

pdf("Figures/Age_Euk_ordi_EFit_sig.pdf")
ggEukSig
dev.off()



############## PERMANOVA ###


EDatNABac <- na.omit(EDatBac,
                     cols=c("cortisol", "mucin", "age_sampling",
                            "Ancylostoma_egg_load", "CSocialRank", "IgA",
                            "season", "sex"))

BacDataNA <- BacData[rownames(EDatNABac),]


EDatNAEuk <- na.omit(EDatEuk,
                     cols=c("mucin", "age_sampling",
                            "Ancylostoma_egg_load", "CSocialRank", 
                            "hyena_ID"))

EukDataNA <- EukData[rownames(EDatNAEuk),]



EukAdonis <- adonis2(EukDataNA ~ mucin + age_sampling + Ancylostoma_egg_load,
                     data=EDatNAEuk, na.action = na.omit, by="margin")


write.csv(round(EukAdonis, 4), "EukAdonis.csv")


BakAdonis <- adonis2(BacDataNA ~ age_sampling + season + sex +
                         mucin + cortisol + CSocialRank +
                         Ancylostoma_egg_load,
        data=EDatNABac, na.action = na.omit, by="margin")

write.csv(round(EukAdonis, 4), "BacAdonis.csv")


## EDatNA_juvenile <- EDatNA[EDatNA$age_sampling<730, ]
## EDatNA_adult <-  EDatNA[EDatNA$age_sampling>=730, ]


## EukDataNA_adult <- EukDataNA[rownames(EDatNA_adult),]
## EukDataNA_juvenile <- EukDataNA[rownames(EDatNA_juvenile),]

## EDatNA_juvenile <- EDatNA[EDatNA$age_sampling<730, ]
## EDatNA_adult <-  EDatNA[EDatNA$age_sampling>=730, ]


## BacDataNA_adult <- BacDataNA[rownames(EDatNA_adult),]
## BacDataNA_juvenile <- BacDataNA[rownames(EDatNA_juvenile),]


## adonis2(BacDataNA_juvenile ~ age_sampling + season +
##             cortisol + CSocialRank + mucin + 
##             neopterin + lysozyme + IgG + IgA + Trichuris_egg_load +
##             Taeniidae_egg_load +Ancylostoma_egg_load  + Cystoisospora_oocyst_load,
##         data=EDatNA_juvenile, na.action = na.omit, by="margin")


## adonis2(BacDataNA_adult ~ age_sampling + season +
##             cortisol + CSocialRank +  mucin + 
##             neopterin + lysozyme + IgG + IgA + Trichuris_egg_load +
##             Taeniidae_egg_load +Ancylostoma_egg_load + Cystoisospora_oocyst_load,
##         data=EDatNA_adult, na.action = na.omit, by="margin")


## adonis2(EukDataNA_juvenile ~ age_sampling + season +
##             season + cortisol + CSocialRank +  mucin + 
##             neopterin + lysozyme + IgG + IgA + Trichuris_egg_load +
##             Taeniidae_egg_load +Ancylostoma_egg_load + Cystoisospora_oocyst_load,
##         data=EDatNA_juvenile, na.action = na.omit, by="margin")


## adonis2(EukDataNA_adult ~ age_sampling + season +
##             cortisol + CSocialRank +  mucin + 
##             neopterin + lysozyme + IgG + IgA + Trichuris_egg_load +
##             Taeniidae_egg_load +Ancylostoma_egg_load + Cystoisospora_oocyst_load,
##         data=EDatNA_adult, na.action = na.omit, by="margin")

