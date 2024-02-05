library(ggplot2)
library(phyloseq)
library(vegan)
#library(tidyverse)
library(cowplot)
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(vegan)
library(phyloseq)
library(igraph)
library(Hmisc)
library(Matrix)
library(SpiecEasi)
library(microbiome)
library(brms)

PMS <- readRDS("/SAN/Susanas_den/gitProj/AA_Hyenas_Pakt/tmp/fPMS.rds")
PMS <- prune_taxa(taxa_sums(PMS)>0, PMS)

do_DataPrep <- FALSE

if (do_DataPrep){## we want to add genetic mum and surrogate mum, in the future this
    
    ## will be done directly in MA. For now, we do it here
    mytable <- read.csv("Data/Covariates_int_biomes_23-06-29.csv")
    metadt <- sample_data(PMS)
    mytable <- mytable[match(metadt$hyena_ID, mytable$hyena_ID),]
    all(metadt$hyena_ID==mytable$hyena_ID) # looks good
    names(mytable)
    metadt$ID_sib <- mytable$ID_sib
    metadt$genetic_mum <- mytable$genetic_mum
    metadt$surrogate_mum <- mytable$surrogate_mum

    ## now some more data cleanig and prep we have some NA's at the
    ## social rank, for juveniles whose genetic mums were dead at hte
    ## time of sampling. We will give them the social rank of their
    ## surrogate mums
    metadt[is.na(metadt$CSocialRank),]
    ## B8423 and B8180 have 2 surrogate mums. Let's give them a mean
    ## social rank.  B8427, B8738 and B8833 have 1 surrogate mum and
    ## we will give them that social rank.
    metadt[metadt$Sample%in%c("B8427", "B8738", "B8833"),
           "CSocialRank"] <-
        metadt[metadt$Sample%in%c("B8427", "B8738", "B8833"),
               "social_rank_surrogate_mum1"]
    metadt[metadt$Sample%in%c("B8423"),]$CSocialRank <-
        mean(c(metadt[metadt$Sample%in%c("B8423"),]$social_rank_surrogate_mum1,
               metadt[metadt$Sample%in%c("B8423"),]$social_rank_surrogate_mum2))
    metadt[metadt$Sample%in%c("B8180"),]$CSocialRank <-
        mean(c(metadt[metadt$Sample%in%c("B8180"),]$social_rank_surrogate_mum1,
               metadt[metadt$Sample%in%c("B8180"),]$social_rank_surrogate_mum2))
    
    ## now let's see if we can complete immune and parasite data with
    ## Susana's dataset
    comp.I <- read.csv("Data/Sample_database_Soares_20230525.csv")
    id_NA <- metadt[is.na(metadt$IgA),]$ hyena_ID
    tb <- subset(comp.I, comp.I$ID%in%id_NA)
    tb <- tb[1:19,]
    tb[, c("ID", "sample_date", "frozen_sample", "formalin_sample")]
    metadt[is.na(metadt$IgA),c("hyena_ID", "date_sampling")]
    #B6493 instead of B6492 (one day later)
    
    metadt[metadt$Sample=="B6492",
           c("Ancylostoma_egg_load",
             "Cystoisospora_oocyst_load",
             "Spirometra_egg_load",
             "Dipylidium_egg_load",
             "Trichuris_egg_load",
             "Taeniidae_egg_load",
             "Spirurida_egg_load",
             "Mesocestoides_egg_load",
             "mucin",
             "IgA",
             "IgG",
             "lysozyme",
             "neopterin")] <- unique(
        tb[tb$formalin_sample=="B6493",
           c("Ancylostoma_egg_load",
             "Cystoisospora_oocyst_load",
             "Spirometra_egg_load",
             "Dipylidium_egg_load",
             "Trichuris_egg_load",
             "Taeniidae_egg_load",
             "Spirurida_egg_load",
             "Mesocestoides_egg_load",
             "mucin",
             "IgA",
             "IgG",
             "lysozyme",
             "neopterin")])

#ES2736 instead of B9482 (8 days later)
    metadt[metadt$Sample=="B9482",
           c("Ancylostoma_egg_load",
             "Cystoisospora_oocyst_load",
             "Spirometra_egg_load",
             "Dipylidium_egg_load",
             "Trichuris_egg_load",
             "Taeniidae_egg_load",
             "Spirurida_egg_load",
             "Mesocestoides_egg_load",
             "mucin",
             "IgA",
             "IgG",
             "lysozyme",
             "neopterin")] <- tb[tb$frozen_sample=="ES2736",
                                 c("Ancylostoma_egg_load",
                                   "Cystoisospora_oocyst_load",
                                   "Spirometra_egg_load",
                                   "Dipylidium_egg_load",
                                   "Trichuris_egg_load",
                                   "Taeniidae_egg_load",
                                   "Spirurida_egg_load",
                                   "Mesocestoides_egg_load", "mucin",
                                   "IgA", "IgG", "lysozyme",
                                   "neopterin")]
#ES887 instead of X6643 (one day later)
    metadt[metadt$Sample=="X6643",
           c("Ancylostoma_egg_load",
             "Cystoisospora_oocyst_load",
             "Spirometra_egg_load",
             "Dipylidium_egg_load",
             "Trichuris_egg_load",
             "Taeniidae_egg_load",
             "Spirurida_egg_load",
             "Mesocestoides_egg_load",
             "mucin",
             "IgA",
             "IgG",
             "lysozyme",
             "neopterin")] <- tb[tb$frozen_sample=="ES887",
                                 c("Ancylostoma_egg_load",
                                   "Cystoisospora_oocyst_load",
                                   "Spirometra_egg_load",
                                   "Dipylidium_egg_load",
                                   "Trichuris_egg_load",
                                   "Taeniidae_egg_load",
                                   "Spirurida_egg_load",
                                   "Mesocestoides_egg_load",
                                   "mucin",
                                   "IgA",
                                   "IgG",
                                   "lysozyme",
                                   "neopterin")]
    
    sample_data(PMS)$AFR[sample_data(PMS)$Survival_Ad=="No"] # this makes sense
    PMS@ sam_data <- metadt
    
    #let's input what we can
    library(mice)
    pMiss <- function(x){sum(is.na(x))/length(x)*100}
    apply(PMS@sam_data,2,pMiss)
    
    # to input: mucin
    df <- data.frame(Sample=PMS@sam_data$Sample, mucin=PMS@sam_data$mucin)
    Mucin <- mice(df, m=5, maxit=50, meth="pmm", seed=500)
    Mucin <- complete(Mucin, 1)
    PMS@sam_data$mucin_inputed<- Mucin$mucin
    
    # to input: IgA
    df <- data.frame(Sample=PMS@sam_data$Sample, IgA=PMS@sam_data$IgA)
    IGA <- mice(df, m=5, maxit=50, meth="pmm", seed=500)
    IGA <- complete(IGA, 1)
    PMS@sam_data$IgA_inputed<- IGA$IgA

    # to input: IgG
    df <- data.frame(Sample=PMS@sam_data$Sample, IgG=PMS@sam_data$IgG)
    IGG <- mice(df, m=5, maxit=50, meth="pmm", seed=500)
    IGG <- complete(IGG, 1)
    PMS@sam_data$IgG_inputed<- IGG$IgG
    
    #defining Clan
    PMS@sam_data$Clan <- substr(PMS@sam_data$hyena_ID, 1,1)
    saveRDS(PMS, "tmp/PMS_imputed.rds")
} else
    PMS <- readRDS("tmp/PMS_imputed.rds")

# let's see which parasites we've got
Apicomplexa <- subset_taxa(PMS, Phylum%in%"Apicomplexa")
Nematodes <- subset_taxa(PMS, Phylum%in%"Nematozoa")
Helminthes <- subset_taxa(PMS, Phylum%in%"Platyhelminthes")

get_taxa_unique(Apicomplexa, "Genus")
get_taxa_unique(Nematodes, "Genus")
get_taxa_unique(Helminthes, "Genus")

#Cysto <- subset_taxa(PMS, Genus%in%c("Eimeriorina", "Toxoplasma", "Eimeria", "Cystoisospora"))
#Cysto <- subset_taxa(PMS, Genus%in%c("Eimeriorina", "Eimeria", "Sarcocystis"))
Cysto <- subset_taxa(PMS, Genus%in%c("Sarcocystis"))
Spiru <- subset_taxa(PMS, Genus%in%c("Spirurida"))
Ancy <- subset_taxa(PMS, Genus%in%c("Rhabditida"))
Diphylo <- subset_taxa(PMS, Genus%in%"Diphyllobothriidea")
Dipy <- subset_taxa(PMS, Genus%in%"Cyclophyllidea")

# let's see the correlation with egg counts
Cysto.df <-
    na.omit(data.frame(Cyso_ASV=sample_sums(tax_glom(Cysto, "Genus")),
                       Cysto_OPG=Cysto@sam_data$
                           Cystoisospora_oocyst_load))

#na.omit(data.frame(Cyso_ASV=sample_sums(tax_glom(Cysto, "Genus")),
#                   Cysto_OPG=Cysto@sam_data$Cystoisospora_oocyst_load))

#na.omit(data.frame(Cyso_ASV=sample_sums(tax_glom(Cysto, "Genus")),
#                   Cysto_OPG=Cysto@sam_data$Cystoisospora_oocyst_load))

#na.omitA(data.frame(Cyso_ASV=sample_sums(tax_glom(Cysto, "Genus")),
#                    Cysto_OPG=Cysto@sam_data$Cystoisospora_oocyst_load))

#na.omit(data.frame(Cyso_ASV=sample_sums(tax_glom(Cysto, "Genus")),
#                   Cysto_OPG=Cysto@sam_data$Cystoisospora_oocyst_load))

Spiru.df <- na.omit(data.frame(ASV=sample_sums(tax_glom(Spiru, "Genus")),
                               OPG=Spiru@sam_data$Spirurida_egg_load))

#na.omit(data.frame(ASV=sample_sums(tax_glom(Spiru, "Genus")),
#                   OPG=Spiru@sam_data$Spirurida_egg_load))

Ancy.df <- na.omit(data.frame(ASV=sample_sums(tax_glom(Ancy, "Genus")),
                              OPG=Spiru@sam_data$Ancylostoma_egg_load))
#na.omit(data.frame(ASV=sample_sums(tax_glom(Ancy, "Genus")),
#                   OPG=Spiru@sam_data$Ancylostoma_egg_load))

Diphylo.df <-na.omit(data.frame(ASV=sample_sums(tax_glom(Diphylo, "Genus")),
                                OPG=Spiru@sam_data$Spirometra_egg_load))
#na.omit(data.frame(ASV=sample_sums(tax_glom(Diphylo, "Genus")),
#                   OPG=Spiru@sam_data$Spirometra_egg_load))

Dipy.df <- na.omit(data.frame(ASV=sample_sums(tax_glom(Dipy, "Genus")),
                              OPG=Dipy@sam_data$Dipylidium_egg_load,
                              OPG_Q=Dipy@sam_data$Dipylidium_rich..includes.proglotids.visualization.info.))
#na.omit(data.frame(ASV=sample_sums(tax_glom(Dipy, "Genus")),
#                   OPG=Dipy@sam_data$Dipylidium_egg_load,
#                   OPG_Q=Dipy@sam_data$Dipylidium_rich..includes.proglotids.visualization.info.))

# Dipylidium is unreliable in egg counts. Let's just see which method is more sensitve
Dipy.df$ASV_Q <- 3
Dipy.df$ASV_Q[Dipy.df$ASV>0] <- 1
Dipy.df$ASV_Q[Dipy.df$ASV==0] <- 0
table(Dipy.df$ASV_Q, Dipy.df$OPG_Q)
# now the correlations
cor.test(Dipy.df$ASV, Dipy.df$OPG)
cor.test(Cysto.df$Cyso_ASV, Cysto.df$Cysto_OPG)
cor.test(Spiru.df$ASV, Spiru.df$OPG)
cor.test(Ancy.df$ASV, Ancy.df$OPG)
cor.test(Diphylo.df$ASV, Diphylo.df$OPG)

#ggplot(Diphylo.df, aes(x=ASV, y=OPG))+
#    geom_point()


# We consider parasites of the gastrointestinal tract of hyenas:
Parasite <- subset_taxa(PMS, Genus%in%c("Sarcocystis", "Spirurida", "Rhabditida", "Diphyllobothriidea", "Cyclophyllidea", "Cystoisospora", "Cryptosporidium", "Ascaridida"))

##  only bacteria
Bacteria <- subset_taxa(PMS, Kingdom %in%"Bacteria") 
Bacteria

# only fungi
Fungi <- subset_taxa(PMS, Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Neocallimastigomycota"))

# and merge per genus
Parasite2 <- tax_glom(Parasite, "Genus")
Parasite2

do_Models <- TRUE

if (do_Models){## we want to add genetic mum and surrogate mum, in the future this
#############################################################################
############# let's do dyad comparisons now
    sample_data(PMS)$key <- paste(sample_data(PMS)$hyena_ID,
                              sample_data(PMS)$Sample, sep="_")
key <- data.frame(ID=sample_data(PMS)$key)
metadt <- sample_data(PMS)
metadt$IgAP <- log(metadt$IgA_inputed)
metadt$MucinP <- log(metadt$mucin_inputed)
# renaming this now
sample_names(PMS) <- key$ID

####################
## 1) Jaccard distance
JACM <- as.matrix(phyloseq::distance(PMS,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
jac<-c(as.dist(JACM))

## 2) Aitchison distance
AITM <- as.matrix(vegan::vegdist(PMS@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
AITM <- 1-AITM
ait<-c(as.dist(AITM))

## 3) Jaccard distance Parasites
JACP <- as.matrix(phyloseq::distance(Parasite,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
JACP[is.na(JACP)]<- 0 # defining those as 0 distances
JACP <- 1-JACP
jacP<-c(as.dist(JACP))

## 4) Aitchison distance Parasites
AITP <- as.matrix(vegan::vegdist(Parasite@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
AITP[is.na(AITP)]<- 0 # defining those as 0 distances
AITP <- 1-AITP
aitP<-c(as.dist(AITP))

## 5) Jaccard distance Bacteria
JACB <- as.matrix(phyloseq::distance(Bacteria,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
JACB[is.na(JACB)]<- 0 # defining those as 0 distances
JACB <- 1-JACB
jacB<-c(as.dist(JACB))

## 6) Aitchison distance Bacteria
AITB <- as.matrix(vegan::vegdist(Bacteria@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
AITB[is.na(AITB)]<- 0 # defining those as 0 distances
AITB <- 1-AITB
aitB<-c(as.dist(AITB))

## 7) Jaccard distance Fungi
JACF <- as.matrix(phyloseq::distance(Fungi,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
JACF[is.na(JACF)]<- 0 # defining those as 0 distances
JACF <- 1-JACF
jacF<-c(as.dist(JACF))

## 8) Aitchison distance Fungi
AITF <- as.matrix(vegan::vegdist(Fungi@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
AITF[is.na(AITF)]<- 0 # defining those as 0 distances
AITF <- 1-AITF
aitF<-c(as.dist(AITF))

## 9) Sex pairs
Sex_frame<-metadt[,c("key","sex")]
Sex_frame$key<-as.character(Sex_frame$key)
Sex_frame$sex<-as.character(Sex_frame$sex)
#Create an empty character matrix to fill with characters
SEXM<-array(as.character(NA),c(nrow(Sex_frame),nrow(Sex_frame)))
for(i in 1:nrow(Sex_frame)){
    for(j in 1:nrow(Sex_frame)){
        if(Sex_frame$sex[i]=="F" & Sex_frame$sex[i]==Sex_frame$sex[j]){
            SEXM[i,j]= "FF"}
        if(Sex_frame$sex[i]=="M" & Sex_frame$sex[i]==Sex_frame$sex[j]){
           SEXM[i,j]= "MM"}
        if( Sex_frame$sex[i]!=Sex_frame$sex[j]){
            SEXM[i,j]= "FM"}
    }
}
sex<-c(SEXM[lower.tri(SEXM)])

# shared clan
ClanM<-array(0,c(nrow(metadt),nrow(metadt)))
for(i in 1:nrow(metadt)){
    for(j in 1:nrow(metadt)){
        if(metadt$Clan[i]==metadt$Clan[j]){
            ClanM[i,j]= "1"
        } else{
            ClanM[i,j]= "0"
        }
    }
}
ClanM<-c(ClanM[lower.tri(ClanM)])  

## now other distances
AGEM <- c(dist(metadt[,c("age_sampling")]))
metadt$date_sampling <- as.Date(metadt$date_sampling,
                                format='%m/%d/%Y')
SamplingM <- c(dist(metadt$date_sampling))
RANKM<-c(dist(metadt[,c("CSocialRank")]))
IGAM<-c(dist(metadt[,c("IgA_inputed")]))
IGAPM<-c(dist(metadt[,c("IgAP")]))
MUCM<-c(dist(metadt[,c("mucin_inputed")]))
MUCPM<-c(dist(metadt[,c("MucinP")]))

# genetic mum
MUM_frame<-metadt[,c("key","genetic_mum")]
MUMM<-array(0,c(nrow(MUM_frame),nrow(MUM_frame)))
for(i in 1:nrow(MUM_frame)){
    for(j in 1:nrow(MUM_frame)){
        if(MUM_frame$genetic_mum[i]==MUM_frame$genetic_mum[j]){
            MUMM[i,j]= "1"
        } else{
            MUMM[i,j]= "0"
        }
    }
}
MUMM <- as.character(c(as.dist(MUMM)))

#Combine these vectors into a data frame
data.dyad<-data.frame(MS_J=jac,
                      MS_A=ait,
                      Parasite_J=jacP,
                      Parasite_A=aitP,
                      Bacteria_J=jacB,
                      Bacteria_A=aitB,
                      Fungi_J=jacF,
                      Fungi_A=aitF,
                      Age=AGEM,
                      Sex=sex,
                      IgA=IGAM,
                      IgAP=IGAPM,
                      Mucin=MUCM,
                      MucinP=MUCPM,
                      Rank=RANKM,
                      Gen_mum=MUMM,
                      Temporal=SamplingM,
                      Clan=ClanM)

list<-expand.grid(key$ID, key$ID)
# This created individual-to-same-individual pairs as well. Get rid of these:
list<-list[which(list$Var1!=list$Var2),]
# this still has both quantiles in--> add 'unique' key
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))
# sanity check that the Individual name combinations are in the same exact order as the lower quantile value vector of the matrices
i=nrow(key)
JACM[which(rownames(JACM)==list$Var1[i]),which(colnames(JACM)==list$Var2[i])]==jac[i]
# add the names of both individuals participating in each dyad into the data frame
data.dyad$IDA_s<-list$Var2
data.dyad$IDB_s<-list$Var1
data.dyad$IDA <- gsub("_.*", "", data.dyad$IDA_s)
data.dyad$IDB <- gsub("_.*", "", data.dyad$IDB_s)

# We keep intra individual comparisons (not the same sample!!)
#data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]

######################### Now we model the data ####################
#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x, min.use, max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("Rank",
             "Age",
             "IgA",
             "IgAP",
             "MucinP",
             "Parasite_A",
             "Bacteria_A",
             "Parasite_A",
             "Fungi_A",
             "MS_A",
             "MS_J",
             "Temporal")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-
        range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
}
    
### Multivariate model for Parasite, Fungi and Bacteria
bformJ <- bf(mvbind(Fungi_J, Parasite_J, Bacteria_J)~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)))+ set_rescor(TRUE)  
fitJ <- brm(bformJ, data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)
saveRDS(fitJ, "tmp/modelJ_FPB.rds")
#
bformA <- bf(mvbind(Fungi_A, Parasite_A, Bacteria_A)~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)))+ set_rescor(TRUE)  
fitA <- brm(bformA, data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)
saveRDS(fitA, "tmp/modelA_FPB.rds")

#### model
modelJ<-brm(MS_J~ 1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 3000,
            cores = 20, chains = 4,
            inits=0)
saveRDS(modelJ, "tmp/modelJ.rds")
#
modelA<-brm(MS_A~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
saveRDS(modelA, "tmp/modelA.rds")
#
#### individual models for plotting
modelFJ <- brm(Fungi_J~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan +Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
saveRDS(modelFJ, "tmp/modelFJ.rds")
#
modelFA <- brm(Fungi_A ~1 + IgAP+ MucinP + Age+  Rank+ Gen_mum+
               Temporal+ Clan+ Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
saveRDS(modelFA, "tmp/modelFA.rds")
#
modelBJ <- brm(Bacteria_J~1+IgAP+ MucinP + Age+  Rank+ Gen_mum+
                   Temporal+ Clan+Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
saveRDS(modelBJ, "tmp/modelBJ.rds")
#
modelBA <- brm(Bacteria_A~1+IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan +Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
saveRDS(modelBA, "tmp/modelBA.rds")
#
modelPJ <- brm(Parasite_J~1+ IgAP+ MucinP + Age+  Rank+ Gen_mum+
                Temporal+ Clan +Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
saveRDS(modelPJ, "tmp/modelPJ.rds")
#
modelPA <- brm(Parasite_A ~1 + IgAP+ MucinP + Age+  Rank+ Gen_mum+
               Temporal+ Clan+ Age:IgAP+Age:MucinP+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
saveRDS(modelPA, "tmp/modelPA.rds")


    
saveRDS(data.dyad, "tmp/data.dyad.rds")
} else
    data.dyad <- readRDS("tmp/data.dyad.rds")

modelJ_FPB <- readRDS("tmp/modelJ_FPB.rds")
modelA_FPB <- readRDS("tmp/modelA_FPB.rds")

modelA <- readRDS("tmp/modelA.rds")
modelJ <- readRDS("tmp/modelJ.rds")
#model_intJ <- readRDS("tmp/model_intJ.rds")
#model_intA <- readRDS("tmp/model_intA.rds")
modelFA <- readRDS("tmp/modelFA.rds")
modelFJ <- readRDS("tmp/modelFJ.rds")
modelBJ <- readRDS("tmp/modelBJ.rds")
modelBA <- readRDS("tmp/modelBA.rds")
modelPJ <- readRDS("tmp/modelPJ.rds")
modelPA <- readRDS("tmp/modelPA.rds")

print(summary(modelJ), digits=3)

print(summary(modelA), digits=3)

print(summary(modelFJ), digits=3)

print(summary(modelPJ), digits=3)

print(summary(modelBJ), digits=3)

print(summary(modelFA), digits=3)

print(summary(modelPA), digits=3)

print(summary(modelBA), digits=3)

print(summary(modelJ_FPB), digits=3) # same as individual models but with autcorr between dependent variables
print(summary(modelA_FPB), digits=3) # same as individual models but with autcorr between dependent variables

############################# plotting
#stil figuring out which is the best vizualizing approach

Parasite
Fungi
Bacteria
PMS


para <- summary(modelJ)$fixed

para

name <- "Overall"

res.fun <- function(modelJ, name){
    para <- summary(modelJ)$fixed
    data.frame(Domain=rep(name, 10),
               Effect=rownames(para),
               Estimate=para$Estimate,
               lCI=para$'l-95% CI',
               uCI=para$'u-95% CI')
}

    
res <- res.fun(modelJ, "Overall")    
res <- rbind(res, res.fun(modelPJ, "Parasites"))
res <- rbind(res, res.fun(modelFJ, "Fungi"))
res <- rbind(res, res.fun(modelBJ, "Bacteria"))

res$Domain <- factor(res$Domain, levels=c("Bacteria", "Fungi", "Parasites", "Overall"))

res <- res[res$Effect%in%c("Age", "IgAP", "MucinP"),]

coul=c("#F8B195","#F67280", "#6C5B7B", "#355C7D")

ggplot(res, aes(x = Estimate, y = Effect, fill = Effect)) +
    geom_errorbar(aes(xmin=lCI, xmax=uCI, colour=Effect), size=1, width=0.4)+
    geom_point(shape = 21, size=3) +
    scale_fill_manual(values = coul) +
    scale_colour_manual(values = coul) +
    xlab("Parameter estimate") +
    ylab("") +
    scale_y_discrete(labels = c("Age dist", "f-IgA dist", "f-mucin dist")) +
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    facet_grid(Domain ~ ., scales = "free_y", space = "free_y")+
    theme_classic()+
    theme(legend.position="none")

resdf.fun<- function(modelJ, name){
    para <- summary(modelJ)$fixed
    data.frame(Domain=name,
               IgA_Estimate=para[rownames(para)=="IgAP", "Estimate"],
               IgA_lCI=para[rownames(para)=="IgAP", "l-95% CI"],
               IgA_uCI=para[rownames(para)=="IgAP", "u-95% CI"],
               Mucin_Estimate=para[rownames(para)=="MucinP", "Estimate"],
               Mucin_lCI=para[rownames(para)=="MucinP", "l-95% CI"],
               Mucin_uCI=para[rownames(para)=="MucinP", "u-95% CI"],
               Age_Estimate=para$Estimate[rownames(para)=="Age"],
               Age_lCI=para[rownames(para)=="Age", "l-95% CI"],
               Age_uCI=para[rownames(para)=="Age", "u-95% CI"],
               Rank_Estimate=para$Estimate[rownames(para)=="Rank"],
               Rank_lCI=para[rownames(para)=="Rank", "l-95% CI"],
               Rank_uCI=para[rownames(para)=="Rank", "u-95% CI"],
               Gen_mum1_Estimate=para$Estimate[rownames(para)=="Gen_mum1"],
               Gen_mum1_lCI=para[rownames(para)=="Gen_mum1", "l-95% CI"],
               Gen_mum1_uCI=para[rownames(para)=="Gen_mum1", "u-95% CI"],
               Temporal_Estimate=para$Estimate[rownames(para)=="Temporal"],
               Temporal_lCI=para[rownames(para)=="Temporal", "l-95% CI"],
               Temporal_uCI=para[rownames(para)=="Temporal", "u-95% CI"],
               Clan1_Estimate=para$Estimate[rownames(para)=="Clan1"],
               Clan1_lCI=para[rownames(para)=="Clan1", "l-95% CI"],
               Clan1_uCI=para[rownames(para)=="Clan1", "u-95% CI"],
               IgA_Age_Estimate=para$Estimate[rownames(para)=="IgAP:Age"],
               IgA_Age_lCI=para[rownames(para)=="IgAP:Age", "l-95% CI"],
               IgA_Age_uCI=para[rownames(para)=="IgAP:Age", "u-95% CI"],
               Mucin_Age_Estimate=para$Estimate[rownames(para)=="MucinP:Age"],
               Mucin_Age_lCI=para[rownames(para)=="MucinP:Age", "l-95% CI"],
               Mucin_Age_uCI=para[rownames(para)=="MucinP:Age", "u-95% CI"]
               )
}

res.df <-resdf.fun(modelJ, "Full")
res.df <- rbind(res.df, resdf.fun(modelBJ, "Bacteria"))
res.df <- rbind(res.df, resdf.fun(modelFJ, "Fungi"))
res.df <- rbind(res.df, resdf.fun(modelPJ, "Parasite"))

res.dfA <-resdf.fun(modelA, "Full")
res.dfA <- rbind(res.dfA, resdf.fun(modelBA, "Bacteria"))
res.dfA <- rbind(res.dfA, resdf.fun(modelFA, "Fungi"))
res.dfA <- rbind(res.dfA, resdf.fun(modelPA, "Parasite"))

res.df$Domain <- factor(res.df$Domain, levels=c("Full", "Bacteria", "Fungi", "Parasite"))
res.dfA$Domain <- factor(res.dfA$Domain, levels=c("Full", "Bacteria", "Fungi", "Parasite"))

coul <- c("#744700", "#d63232", "#ffbf00", "#293885")

IgA<-ggplot(res.df, aes(x=IgA_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=IgA_lCI, xmax=IgA_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="f-IgA distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

IgAA<-ggplot(res.dfA, aes(x=IgA_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=IgA_lCI, xmax=IgA_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="f-IgA distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

Mucin<-ggplot(res.df, aes(x=Mucin_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=Mucin_lCI, xmax=Mucin_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="f-mucin distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

MucinA<-ggplot(res.dfA, aes(x=Mucin_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=Mucin_lCI, xmax=Mucin_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="f-mucin distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

Age<-ggplot(res.df, aes(x=Age_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=Age_lCI, xmax=Age_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Age distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")
AgeA<-ggplot(res.dfA, aes(x=Age_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=Age_lCI, xmax=Age_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Age distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

Rank<-ggplot(res.df, aes(x=Rank_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=Rank_lCI, xmax=Rank_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Social rank distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

Gen_Mum<-ggplot(res.df, aes(x=Gen_mum1_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=Gen_mum1_lCI, xmax=Gen_mum1_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Same genetic mother", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

Temporal<-ggplot(res.df, aes(x=Temporal_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=Temporal_lCI, xmax=Temporal_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Temporal distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

Clan1<-ggplot(res.df, aes(x=Clan1_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=Clan1_lCI, xmax=Clan1_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Same clan", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

Figure2 <- plot_grid(IgA, IgAA, Mucin, MucinA, Age, AgeA, labels="auto", ncol=2)
ggplot2::ggsave(file="Figures/Figure2.pdf", Figure2, width = 190, height = 150, dpi = 300, units="mm")

# interaction
library(tidyr)
library(modelr)
library(tidybayes)

######################## Figure 3: interaction
newdata1 <- data.frame(IgAP=seq_range(0:1, n=51),
                       MucinP=rep(median(data.dyad$MucinP), n=51),
                       Age=rep(1, n=51),
                       Clan=rep(0, n=51),
                       IDA=rep("M668", 51),
                       IDB=rep("M604", 51),
                       Gen_mum=rep(0, n=51),
                       Temporal=rep(0, n=51),
                       Rank=rep(median(data.dyad$Rank)))
newdata0 <- data.frame(IgAP=seq_range(0:1, n=51),
                       MucinP=rep(median(data.dyad$MucinP), n=51),
                       Age=rep(0, n=51),
                       Clan=rep(0, n=51),
                       IDA=rep("M668", 51),
                       IDB=rep("M604", 51),
                       Gen_mum=rep(0, n=51),
                       Temporal=rep(0, n=51),
                       Rank=rep(median(data.dyad$Rank)))
pred.df0 <- add_epred_draws(newdata0, modelPJ)
pred.df1 <- add_epred_draws(newdata1, modelPJ)


pred.df <- rbind(pred.df0, pred.df1)
pred.df$Age <- as.factor(pred.df$Age)

IgA_Parasite <-    ggplot(pred.df, aes(x=IgAP, y=.epred, group=Age))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5)+
    ylab("Parasite community similarity")+
    xlab("faecal IgA distances")+
    theme_classic()

pred.df0B <- add_epred_draws(newdata0, modelBJ)
pred.df1B <- add_epred_draws(newdata1, modelBJ)

pred.dfB <- rbind(pred.df0B, pred.df1B)
pred.dfB$Age <- as.factor(pred.dfB$Age)

IgA_Bacteria <-    ggplot(pred.dfB, aes(x=IgAP, y=.epred, group=Age))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5)+
    ylab("Bacteria community similarity")+
    xlab("faecal IgA distances")+
    theme_classic()


ab <- plot_grid(IgA_Parasite, IgA_Bacteria, labels="auto")

newdata1M <- data.frame(MucinP=seq_range(0:1, n=51),
                       IgAP=rep(median(data.dyad$IgAP), n=51),
                       Age=rep(1, n=51),
                       Clan=rep(0, n=51),
                       IDA=rep("M668", 51),
                       IDB=rep("M604", 51),
                       Gen_mum=rep(0, n=51),
                       Temporal=rep(0, n=51),
                       Rank=rep(median(data.dyad$Rank)))
newdata0M <- data.frame(MucinP=seq_range(0:1, n=51),
                       IgAP=rep(median(data.dyad$IgAP), n=51),
                       Age=rep(0, n=51),
                       Clan=rep(0, n=51),
                       IDA=rep("M668", 51),
                       IDB=rep("M604", 51),
                       Gen_mum=rep(0, n=51),
                       Temporal=rep(0, n=51),
                       Rank=rep(median(data.dyad$Rank)))

pred.df0M <- add_epred_draws(newdata0M, modelPA)
pred.df1M <- add_epred_draws(newdata1M, modelPA)
pred.dfM <- rbind(pred.df0M, pred.df1M)
pred.dfM$Age <- as.factor(pred.dfM$Age)


Mucin_Parasite <- ggplot(pred.dfM, aes(x=MucinP, y=.epred, group=Age))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5)+
    ylab("Parasite community similarity")+
    xlab("faecal mucin distances")+
    theme_classic()

pred.df0MF <- add_epred_draws(newdata0M, modelFA)
pred.df1MF <- add_epred_draws(newdata1M, modelFA)
pred.dfMF <- rbind(pred.df0MF, pred.df1MF)
pred.dfMF$Age <- as.factor(pred.dfMF$Age)

Mucin_Fungi <- ggplot(pred.dfMF, aes(x=MucinP, y=.epred, group=Age))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5)+
    ylab("Fungi community similarity")+
    xlab("faecal mucin distances")+
    theme_classic()

cd <- plot_grid(Mucin_Parasite, Mucin_Fungi, labels=c("c", "d"))

Fig3 <- plot_grid(ab, cd, nrow=1)

ggplot2::ggsave(file="Figures/Figure3.pdf", Fig3, width = 190, height = 140, dpi = 300, units="mm")

#############################
# random forest regression

library(caret)
library(ranger)
library(pdp)
library(patchwork)

PMS.t <- transform(PMS, "compositional") # transform to make all taxa from 0 to 1
otu <- PMS.t@otu_table

colnames(otu) <- paste("ASV", seq(1, length(colnames(otu))), PMS@tax_table[,6], sep="_")
df <- data.frame(otu, IgA=PMS.t@sam_data$IgA_inputed)
df$IgA <- log(df$IgA)

set.seed(56)
trainIndex <- createDataPartition(df$IgA, p=0.8, list=FALSE, times=1)
IgA_df_train <- df[trainIndex,]
IgA_df_test <- df[-trainIndex,]

set.seed(12)
fitControl <- trainControl(method="repeatedcv", number=10, repeats=10)

rfFit1 <- train(IgA~., data=IgA_df_train,
                method="ranger",
                trControl=fitControl,
                importance="permutation")


print(rfFit1)
print(rfFit1$finalModel)


test_predictions <- predict(rfFit1, newdata=IgA_df_test)
print(postResample(test_predictions, IgA_df_test$IgA))

pred_obs <- data.frame(predicted=test_predictions, observed=IgA_df_test$IgA)

corr <- ggplot(pred_obs, aes(x=predicted, y=observed))+
    geom_point(size=3, color="orange")+
    geom_abline(linetype=5, color="blue", size=1)+
    theme_classic()

imp <- varImp(rfFit1)
imp <- imp$importance
imp$taxa <- rownames(imp)
rownames(imp) <- NULL
imp <- imp[order(-imp$Overall),]
imp20 <- imp[1:20,]
imp20 <- droplevels(imp20)
imp20$taxa <- with(imp20, reorder(taxa, Overall))

topImp <- ggplot(imp20, aes(y=taxa, x=Overall))+
    geom_segment( aes(yend=taxa, xend=0)) +
    geom_point(size=4, color="orange")+
    labs(x="importance", y="")+
    theme_classic()

Fig4 <- plot_grid(corr, topImp, labels="auto", rel_widths=c(0.6, 1))

top_features <- imp20$taxa

pd_plots <- list(NULL)

top_features <- as.character(top_features)

for (a in 1:length(top_features)) {
    pd_plots[[a]] <-pdp::partial(rfFit1, pred.var=top_features[a], rug=TRUE)%>%
        autoplot()+
        geom_hline(yintercept = mean(IgA_df_train$IgA), linetype = 2, color = "gray") + 
#            scale_y_continuous(limits=c(1.5,2.3)) + # Harmonize the scale of yhat on all plots
        theme(panel.border = element_rect(colour = "black", fill = NA),
                      panel.background = element_blank())
    print(paste0("Partial dependence of ", top_features[a]))
}

fig4_2 <- wrap_plots(pd_plots, ncol=4)
fig4 <- plot_grid(Fig4, fig4_2, nrow=2, rel_heights=c(0.5, 1))
ggplot2::ggsave(file="Figures/Figure4.pdf", fig4, width = 200, height = 250, dpi = 300, units="mm")


#plot(varImp(rfFit1))


# for mucin
df <- data.frame(otu, Mucin=PMS.t@sam_data$mucin_inputed)
df$Mucin <- log(df$Mucin)

set.seed(57)
trainIndex <- createDataPartition(df$Mucin, p=0.8, list=FALSE, times=1)
Mucin_df_train <- df[trainIndex,]
Mucin_df_test <- df[-trainIndex,]

set.seed(13)
fitControl <- trainControl(method="repeatedcv", number=10, repeats=10)

rfFitM <- train(Mucin~., data=Mucin_df_train,
                method="ranger",
                trControl=fitControl,
                importance="permutation")


print(rfFitM)
print(rfFitM$finalModel)


test_predictions <- predict(rfFitM, newdata=Mucin_df_test)
print(postResample(test_predictions, Mucin_df_test$Mucin))

pred_obs <- data.frame(predicted=test_predictions, observed=Mucin_df_test$Mucin)

corrM <- ggplot(pred_obs, aes(x=predicted, y=observed))+
    geom_point(size=3, color="orange")+
    geom_abline(linetype=5, color="blue", size=1)+
    theme_classic()

impM <- varImp(rfFitM)
impM <- impM$importance
impM$taxa <- rownames(impM)
rownames(impM) <- NULL
impM <- impM[order(-impM$Overall),]
impM20 <- impM[1:20,]
impM20 <- droplevels(impM20)
impM20$taxa <- with(impM20, reorder(taxa, Overall))

topMImp <- ggplot(impM20, aes(y=taxa, x=Overall))+
    geom_segment( aes(yend=taxa, xend=0)) +
    geom_point(size=4, color="orange")+
    labs(x="importance", y="")+
    theme_classic()

Fig5 <- plot_grid(corrM, topMImp, labels="auto", rel_widths=c(0.6, 1))

top_features <- impM20$taxa

pd_plots <- list(NULL)

top_features <- as.character(top_features)

for (a in 1:length(top_features)) {
    pd_plots[[a]] <-pdp::partial(rfFitM, pred.var=top_features[a], rug=TRUE)%>%
        autoplot()+
        geom_hline(yintercept = mean(Mucin_df_train$Mucin), linetype = 2, color = "gray") + 
#            scale_y_continuous(limits=c(1.5,2.3)) + # Harmonize the scale of yhat on all plots
        theme(panel.border = element_rect(colour = "black", fill = NA),
                      panel.background = element_blank())
    print(paste0("Partial dependence of ", top_features[a]))
}

fig5_2 <- wrap_plots(pd_plots, ncol=4)
fig5 <- plot_grid(Fig5, fig5_2, nrow=2, rel_heights=c(0.5, 1))
ggplot2::ggsave(file="Figures/Figure5.pdf", fig5, width = 200, height = 250, dpi = 300, units="mm")



#####################################################################
###### network
Euk <- merge_phyloseq(Parasite, Fungi)
Euk
Bacteria

## prevalebce filtering of 10%
KeepTaxap <- microbiome::prevalence(Euk)>0.1
Euk<- phyloseq::prune_taxa(KeepTaxap, Euk)
Euk

KeepTaxap <- microbiome::prevalence(Bacteria)>0.1
Bac<- phyloseq::prune_taxa(KeepTaxap, Bacteria)
Bac

#### spiec easi
pargs <- list(rep.num=1000, seed=10010, ncores=90, thresh=0.05)

## mb
t1 <- Sys.time()
se.net <- spiec.easi(list(Bac, Euk), method="mb", pulsar.params=pargs)
t2 <- Sys.time()
t2-t1
saveRDS(se.net, "tmp/se.fnet.rds")
se.net <- readRDS("tmp/se.fnet.rds")

se.net$select$stars$summary # lambda path


se.net$select

#coding bacteria/eukaryote nodes
dtype <- c(rep(1,ntaxa(Bac)), rep(2,ntaxa(Euk)))

bac.ids=taxa_names(Bac)
euk.ids= taxa_names(Euk)
net.ids <- c(bac.ids,euk.ids)

#### plotting
bm=symBeta(getOptBeta(se.net), mode="maxabs")
diag(bm) <- 0
#weights <- Matrix::summary(t(bm))[,3] # includes negative weights
weights <- (1-Matrix::summary(t(bm))[,3])/2 # ort
net <- SpiecEasi::adj2igraph(Matrix::drop0(getRefit(se.net)),
                             edge.attr=list(weight=weights),
                             vertex.attr = list(name=net.ids))

E(net)
betaMat=as.matrix(symBeta(getOptBeta(se.net)))
# we want positive edges to be green and negative to be red
edges <- E(net)
edge.colors=c()
for (e.index in 1:length(edges)){
    adj.nodes=ends(net, edges[e.index])
    xindex=which(net.ids==adj.nodes[1])
    yindex=which(net.ids==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta>0){
        edge.colors=append(edge.colors, "#1B7837")
    }else if(beta<0){
        edge.colors=append(edge.colors, "#762A83")
    }
}
E(net)$color=edge.colors

### defining attributes
V(net)$type=c(rep("Bacteria", length(taxa_names(Bac))), rep("Eukaryote", length(taxa_names(Euk))))
V(net)$genus=c(Bac@tax_table[,6], Euk@tax_table[,6])
V(net)$family=c(Bac@tax_table[,5], Euk@tax_table[,5])
V(net)$phylum=c(Bac@tax_table[,2], Euk@tax_table[,2])
V(net)$domain=c(Bac@tax_table[,1], Euk@tax_table[,1])

# curating names
V(net)$phylum[V(net)$phylum=="Unknown_phylum"] <- "Unknown"
V(net)$genus <- gsub("Unknown_genus_in_", "", V(net)$genus)
V(net)$genus <- gsub("_uncultered", "", V(net)$genus)


deg <- igraph::degree(net)
label <- V(net)$genus
label[deg<4] <- ""

close <- closeness(net)
betw <- betweenness(net)

vip <- data.frame(degree=deg, closeness=close, betweness=betw, genus=V(net)$genus, phyla=V(net)$phylum)

rownames(vip) <- NULL

V(net)$stype <- c(rep("circle",ntaxa(Bac)), rep("square",ntaxa(Euk)))

#V(net)$lab.hub <- ""
#V(net)$lab.hub <- V(net)$genus
#V(net)$label.cex <- 0.5
#V(net)$label.dist <- 0
#V(net)$label.degree <- pi/2

# we also want the node color to code for phylum
nb.col <- length(levels(as.factor(V(net)$phylum)))
coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)
mc <- coul[as.numeric(as.factor(V(net)$phylum))]


pdf("Figures/Network.pdf",
    width =10, height = 10)
set.seed(1002)
plot(net,
     layout=layout_with_fr(net),
     vertex.shape=V(net)$stype,
     vertex.label="",
     vertex.label.dist=0.4,
     vertex.label.degree=-pi/2,
     vertex.size=deg+3,
     vertex.color=adjustcolor(mc,0.8),
     edge.width=as.integer(cut(E(net)$weight, breaks=6))/3,
     margin=c(0,1,0,0))
legend(x=-2, y=1, legend=levels(as.factor(V(net)$phylum)), col=coul, bty="n",x.intersp=0.25,text.width=0.045, pch=20, pt.cex=1.5)
dev.off()

