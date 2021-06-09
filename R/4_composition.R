### FOR COMPOSITION WE USE NO NORMAIZATION, but we log10 transform
### counts!!!


############ Preparing datasets for ordination
prune_prune <- function (ps, samples=0, taxa=10) {
    s <- prune_samples(sample_sums(ps) > samples, ps)
    S <- prune_samples(!grepl("Negative", sample_names(s)), s)
    prune_taxa(taxa_sums(S) > taxa, S)
}

## ## genus glomed data wont converge for bacteria...
## logBAC <- transform_sample_counts(
##     prune_prune(PBac),
##     function(x) log10(1+x))

## trying the species glomed dataset
logBac <- transform_sample_counts(
    prune_prune(PBac, -1, 50),
    function(x) log10(1+x))

logEuk <- transform_sample_counts(
    prune_prune(PEuk, 0, 50),
    function(x) log10(1+x))

library(vegan)

######### Prepare environmental data
EDatraw <- sample_data(logBac)
class(EDatraw) <- "data.frame"

immuno <- c("cortisol", "neopterin", "lysozyme", "IgG", "IgA", "mucin")

para <- c("Ancylostoma_egg_load", "Cystoisospora_oocyst_load",
          "Spirometra_egg_load", "Dipylidium_egg_load",
          "Trichuris_egg_load", "Taeniidae_egg_load",
          "Spirurida_egg_load", "Mesocestoides_egg_load")

lhist <- c("season", "sampling_month",
           "prey_level_imputed", # "prey_level",
           "CSocialRank",
           "age_sampling", "sex", "clan")

EDat <- EDatraw[, c(immuno, para, lhist)]


### Bacterial nMDS ##############################################
BacData <- otu_table(logBAC)
colnames(BacData) <- as.vector(tax_table(logBAC)[, "species"])

BacDist <- distance(BacData, method="bray")

mdsBac <- vegan::metaMDS(BacDist, try=150, trymax=150, k=3)

BacEFit <- envfit(mdsBac, EDat, na.rm=TRUE)
BacEFit

### Eukaryote nMDS ###########################################

EukData <- otu_table(logEuk)
colnames(EukData) <- as.vector(tax_table(logEuk)[, "genus"])

EukDist <- distance(EukData, method="bray")

mdsEuk <- vegan::metaMDS(EukDist, try=150, trymax=150, k=3)

EukEFit <- envfit(mdsEuk, EDat, na.rm=TRUE)
EukEFit


### exporting the data for plotting with ggplot

data.scoresEuk <-  as.data.frame(scores(mdsEuk))
data.scoresEuk <- cbind(data.scoresEuk, EDat)

EnCoordEuk <-
    as.data.frame(
        rbind(scores(EukEFit, "vectors") * ordiArrowMul(EukEFit),
              scores(EukEFit, "factors") * ordiArrowMul(EukEFit))
    )

EnCoordEuk$Cat <- ifelse(rownames(EnCoordEuk)%in%immuno, "Immune",
                  ifelse(rownames(EnCoordEuk)%in%para, "Parasite", "other"))


data.scoresBac <-  as.data.frame(scores(mdsBac))
data.scoresBac <- cbind(data.scoresBac, EDat)

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

pdf("Figures/Age_Bac_ordi_EFit_cut.pdf")
plot(ggBac) + xlim( -0.0826, 0.01)
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

pdf("Figures/Age_Bac_ordi_EFit_Immune_cut.pdf")
ggBacIm +
        xlim( -0.0826, 0.01)
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

pdf("Figures/Age_Bac_ordi_EFit_Parasite_cut.pdf")
ggBacPara +
    xlim( -0.0826, 0.01)
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

pdf("Figures/Age_Bac_ordi_EFit_other_cut.pdf")
ggBacOther +
    xlim( -0.0826, 0.01)
dev.off()


ggEuk <-  ggplot(data = data.scoresEuk, aes(x = NMDS1, y = NMDS2)) +
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





############## PERMANOVA ###


EDatNA <- na.omit(EDat, cols=c("cortisol", "mucin", "age_sampling", "CSocialRank", "season"))

EukDataNA <- EukData[rownames(EDatNA),]

BacDataNA <- BacData[rownames(EDatNA),]


adonis2(EukDataNA ~ age_sampling + season + cortisol + CSocialRank + mucin + 
            neopterin + lysozyme + IgG + IgA + mucin + Trichuris_egg_load +
            Taeniidae_egg_load +Ancylostoma_egg_load + Cystoisospora_oocyst_load,
        data=EDatNA, na.action = na.omit)


adonis2(EukDataNA ~ age_sampling + season + cortisol + CSocialRank + mucin + 
            neopterin + lysozyme + IgG + IgA + mucin + Trichuris_egg_load +
            Taeniidae_egg_load +Ancylostoma_egg_load + Cystoisospora_oocyst_load,
        data=EDatNA, na.action = na.omit, by="margin")


adonis2(BacDataNA ~ age_sampling + season + cortisol + CSocialRank + mucin + 
            neopterin + lysozyme + IgG + IgA + mucin + Trichuris_egg_load +
            Taeniidae_egg_load +Ancylostoma_egg_load + Cystoisospora_oocyst_load,
        data=EDatNA, na.action = na.omit)


adonis2(BacDataNA ~ age_sampling + season + cortisol + CSocialRank + mucin + 
            neopterin + lysozyme + IgG + IgA + mucin + Trichuris_egg_load +
            Taeniidae_egg_load +Ancylostoma_egg_load + Cystoisospora_oocyst_load,
        data=EDatNA, na.action = na.omit, by="margin")


EDatNA_juvenile <- EDatNA[EDatNA$age_sampling<730, ]
EDatNA_adult <-  EDatNA[EDatNA$age_sampling>=730, ]


EukDataNA_adult <- EukDataNA[rownames(EDatNA_adult),]
EukDataNA_juvenile <- EukDataNA[rownames(EDatNA_juvenile),]

EDatNA_juvenile <- EDatNA[EDatNA$age_sampling<730, ]
EDatNA_adult <-  EDatNA[EDatNA$age_sampling>=730, ]


BacDataNA_adult <- BacDataNA[rownames(EDatNA_adult),]
BacDataNA_juvenile <- BacDataNA[rownames(EDatNA_juvenile),]



adonis2(BacDataNA_juvenile ~ age_sampling + season +
            cortisol + CSocialRank + mucin + 
            neopterin + lysozyme + IgG + IgA + Trichuris_egg_load +
            Taeniidae_egg_load +Ancylostoma_egg_load  + Cystoisospora_oocyst_load,
        data=EDatNA_juvenile, na.action = na.omit, by="margin")


adonis2(BacDataNA_adult ~ age_sampling + season +
            cortisol + CSocialRank +  mucin + 
            neopterin + lysozyme + IgG + IgA + Trichuris_egg_load +
            Taeniidae_egg_load +Ancylostoma_egg_load + Cystoisospora_oocyst_load,
        data=EDatNA_adult, na.action = na.omit, by="margin")


adonis2(EukDataNA_juvenile ~ age_sampling + season +
            season + cortisol + CSocialRank +  mucin + 
            neopterin + lysozyme + IgG + IgA + Trichuris_egg_load +
            Taeniidae_egg_load +Ancylostoma_egg_load + Cystoisospora_oocyst_load,
        data=EDatNA_juvenile, na.action = na.omit, by="margin")


adonis2(EukDataNA_adult ~ age_sampling + season +
            cortisol + CSocialRank +  mucin + 
            neopterin + lysozyme + IgG + IgA + Trichuris_egg_load +
            Taeniidae_egg_load +Ancylostoma_egg_load + Cystoisospora_oocyst_load,
        data=EDatNA_adult, na.action = na.omit, by="margin")

