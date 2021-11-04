library(microbiome)

if(!exists("PBac")){
    source("R/1_HyenaPakt_MA.R")
}

############ Preparing datasets for ordination
prune_prune <- function (ps, samples=0, taxa=10) {
    s <- prune_samples(sample_sums(ps) > samples, ps)
    S <- prune_samples(!grepl("Negative", sample_names(s)), s)
    prune_taxa(taxa_sums(S) > taxa, S)
}

### FOR COMPOSITION WE USE proportional transformation
## trying the species glomed dataset
cBac <- microbiome::transform(
                        prune_prune(PBac, 20, 10),
                        "compositional")


## cEuk <- subset_samples(PEuk, !sample_names(cEuk)%in%c("B4285", "B9157"))


cEuk <- microbiome::transform(
                        prune_prune(PEuk, 40, 20),
                        "compositional")

## cEuk <- subset_taxa(PEuk,
##                      !phylum%in%c("Nematoda",
##                                   "Apicomplexa",
##                                   "Platyhelminthes"))


## cEuk <- subset_taxa(cEuk, !(is.na(phylum)))
    



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

fitn <- c("AFR", "LRS", "age_death", "Survival_Ad")

EDatBac <- EDatrawBac[, c(immuno, para, lhist, fitn)]
EDatEuk <- EDatrawEuk[, c(immuno, para, lhist, fitn)]


### Bacterial nMDS ##############################################
BacData <- otu_table(cBac)
colnames(BacData) <- as.vector(tax_table(cBac)[, "species"])

mdsBac <- vegan::metaMDS(BacData, try=350, trymax=350, k=3)

BacEFit <- envfit(mdsBac, EDatBac, na.rm=TRUE)
BacEFit

### Eukaryote nMDS ###########################################

### Removing ASVs from the target taxa... to not associate those with
### themselves


EukData <- otu_table(cEuk)

mdsEuk <- vegan::metaMDS(EukData, try=250, trymax=250, k=3)

EukEFit <- envfit(mdsEuk, EDatEuk, na.rm=TRUE)
EukEFit


### exporting the data for plotting with ggplot

data.scoresEuk <-  as.data.frame(scores(mdsEuk))
data.scoresEuk <- cbind(data.scoresEuk, EDatEuk[rownames(data.scoresEuk),])

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

data.scoresEuk$lab <- rownames(data.scoresEuk)

ggEuk <-  ggplot(data = data.scoresEuk, aes(x = NMDS1, y = NMDS2)) +
    geom_label(data = data.scoresEuk, aes(colour = age_sampling, 
                                          label=lab),
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

EukSig

EDatNAEuk <- na.omit(EDatEuk,
                     cols=c("mucin", "neopterin", "IgA", "age_sampling", "CSocialRank", 
                            "hyena_ID"))

EukDataNA <- EukData[rownames(EDatNAEuk),]



EukAdonis <- adonis2(EukDataNA ~ mucin + neopterin + IgA +
                         age_sampling + CSocialRank,
                     data=EDatNAEuk, na.action = na.omit, by="margin")


write.csv(round(EukAdonis, 4), "EukAdonis.csv")


BacAdonis <- adonis2(BacDataNA ~ age_sampling + season + sex +
                         mucin + cortisol + CSocialRank +
                         Ancylostoma_egg_load,
        data=EDatNABac, na.action = na.omit, by="margin")

write.csv(round(BacAdonis, 4), "BacAdonis.csv")


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



## partial least squares regression

library(caret)
library(pls)

foo <- cbind(EukData, sample_data(P)[rownames(EukData), "LRS.y"])
foo <- foo[!is.na(foo$LRS.y), ]

set.seed(123)
training.samples <- sample(rownames(foo), size=0.8*nrow(foo))


train.data  <- foo[training.samples, ]
test.data <- foo[!rownames(foo)%in%training.samples, ]


set.seed(123)
model <- train(
    LRS.y~., data = train.data, method = "pls",
    scale = FALSE,
    trControl = trainControl("cv", number = 10),
    tuneLength = 10
)
## Plot model RMSE vs different values of components
plot(model)
## Print the best tuning parameter ncomp that
## minimize the cross-validation error 

summary(model$finalModel)

## Make predictions
predictions <- model %>% predict(test.data)
## Model performance metrics

data.frame(
    RMSE = caret::RMSE(predictions, test.data$LRS.y),
    Rsquare = caret::R2(predictions, test.data$LRS.y)
)

ggplot(data.frame(predictions, test.data$LRS.y),
       aes(predictions, test.data.LRS.y)) +
    geom_point()
       


foo <- cbind(EukData, sample_data(P)[rownames(EukData), "Survival_Ad.y"])
foo <- foo[!is.na(foo$Survival_Ad.y), ]

set.seed(123)
training.samples <- sample(rownames(foo), size=0.8*nrow(foo))


train.data  <- foo[training.samples, ]
test.data <- foo[!rownames(foo)%in%training.samples, ]


set.seed(123)
model <- train(
    Survival_Ad.y~., data = train.data, method = "pls",
    scale = FALSE,
    preProc= "center",
    trControl = trainControl("cv", number = 10),
    tuneLength = 10
)
## Plot model RMSE vs different values of components
plot(model)
## Print the best tuning parameter ncomp that
## minimize the cross-validation error 

summary(model$finalModel)
