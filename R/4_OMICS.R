library(ggplot2)
library(phyloseq)
library(vegan)
#library(tidyverse)

PMS <- readRDS("/SAN/Susanas_den/gitProj/AA_Hyenas_Pakt/tmp/PMS_imputed.rds")
metab <- read.csv("/SAN/Susanas_den/gitProj/AA_Hyenas_Pakt/Data/microbiome_samples_tagged_life_history_fitness_immunity_metabolomics.csv")

# joining ID and sample date for merging
PMS@sam_data$ID_age <- paste(PMS@sam_data$hyena_ID, PMS@sam_data$age_sampling, sep="_")
metab$ID_age <- paste(metab$hyena_ID, metab$age_sampling_days, sep="_")
# subsetting
PMS_met <- subset_samples(PMS, ID_age%in%metab$ID_age)
#subseting
metab <- metab[metab$ID_age%in%PMS_met@sam_data$ID_age,]
#ordering
metab <- metab[match(PMS_met@sam_data$ID_age, metab$ID_age),]
# sanity check
all(metab$ID_age==PMS_met@sam_data$ID_age)
# now only metabolites from metab table
metabolite <- metab[,71:168]

all(is.na(metab)) # sanity check

#### procrustes analysis

aitMi <- vegan::vegdist(PMS_met@otu_table,
                                 method="aitchison",
                                 pseudocount=1)

aitMe <- vegan::vegdist(metabolite,
                        method="aitchison")


library(ade4)

add <-  !(is.euclid(aitMi))
pcoa.mi <- cmdscale(aitMi, k = nrow(PMS_met@otu_table)-1, eig = TRUE, add = add)
add <-  !(is.euclid(aitMi))
pcoa.me <- cmdscale(aitMe, k = nrow(metabolite)-1, eig = TRUE, add = add)

pro <- procrustes(aitMi, aitMe, symmetric=TRUE)

pdf("Figures/Procrustes.pdf",width = 4, height = 4)
plot(pro, kind=1, type="p")
dev.off()

plot(pro, kind=2)

protest(X = aitMi, Y = aitMe, scores = "sites", permutations = 999)

## distances
aitMic <- as.matrix(vegan::vegdist(PMS_met@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
aitMic <- 1-aitMic
aitMic<-c(as.dist(aitMic))

aitMet <- as.matrix(vegan::vegdist(metabolite,
                                 method="aitchison"))
aitMet <- 1-aitMet
aitMet<-c(as.dist(aitMet))



# We consider parasites of the gastrointestinal tract of hyenas:
Parasite <- subset_taxa(PMS_met, Genus%in%c("Sarcocystis", "Spirurida", "Rhabditida", "Diphyllobothriidea", "Cyclophyllidea", "Cystoisospora", "Cryptosporidium", "Ascaridida"))

aitP <- as.matrix(vegan::vegdist(Parasite@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
aitP[is.na(aitP)]<- 0 # defining those as 0 distances
aitP <- 1-aitP
aitP<-c(as.dist(aitP))


##  only bacteria
Bacteria <- subset_taxa(PMS_met, Kingdom %in%"Bacteria")
aitB <- as.matrix(vegan::vegdist(Bacteria@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
aitB[is.na(aitB)]<- 0 # defining those as 0 distances
aitB <- 1-aitB
aitB<-c(as.dist(aitB))

# only fungi
Fungi <- subset_taxa(PMS_met, Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Neocallimastigomycota"))
aitF <- as.matrix(vegan::vegdist(Fungi@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
aitF[is.na(aitF)]<- 0 # defining those as 0 distances
aitF <- 1-aitF
aitF<-c(as.dist(aitF))

key <- as.factor(sample_data(PMS_met)$ID_age)
mydf <- sample_data(PMS_met)

AgeD <- c(dist(mydf[,c("age_sampling")]))

data.dyad<-data.frame(Microbiome_A=aitMic,
                      Metabolome_A=aitMet,
                      AgeDiff=AgeD,
                      Bacteria_A=aitB,
                      Fungi_A=aitF,
                      Parasite_A=aitP)

list<-expand.grid(key, key)

# This created individual-to-same-individual pairs as well. Get rid of these:
list<-list[which(list$Var1!=list$Var2),]

# this still has both quantiles in--> add 'unique' key
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))
# sanity check that the Individual name combinations are in the same exact order as the lower quantile value vector of the matrices

# add the names of both individuals participating in each dyad into the data frame
data.dyad$IDA_s<-list$Var2
data.dyad$IDB_s<-list$Var1

data.dyad$IDA <- gsub("_.*", "", data.dyad$IDA_s)
data.dyad$IDB <- gsub("_.*", "", data.dyad$IDB_s)

# Make sure you have got rid of all self comparisons
data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]

######################### Now we model the data ####################
#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x, min.use, max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("Microbiome_A",
             "AgeDiff",
             "Metabolome_A",
             "Bacteria_A",
             "Parasite_A",
             "Fungi_A")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-
        range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }

library(brms)

#### model
modelMic_Met<-brm(Metabolome_A~ 1+ Microbiome_A+ AgeDiff+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 5000,
            cores = 20, chains = 4,
            inits=0)

model_Met<-brm(Metabolome_A~ 1+ Bacteria_A+ Fungi_A + Parasite_A+ AgeDiff+
                (1|mm(IDA,IDB)),
            data = data.dyad,
            family= "gaussian",
            warmup = 1000, iter = 5000,
            cores = 20, chains = 4,
            inits=0)

saveRDS(model_Met, "tmp/model_Met.rds")

model_Met <- readRDS("tmp/model_Met.rds")

print(summary(model_Met), digits=3)



###################Faeces color and consistency #########################################
#####################################################################
metab$faeces_consistency

metab$fcon_ordinal <- 20
metab$fcon_ordinal[metab$faeces_consistency=="dry"] <- 1
metab$fcon_ordinal[metab$faeces_consistency=="homogeneous"] <- 2
metab$fcon_ordinal[metab$faeces_consistency=="moderate-soft"] <- 3
metab$fcon_ordinal[metab$faeces_consistency=="very_soft"] <- 4
metab$fcon_ordinal[metab$faeces_consistency=="liquid"] <- 4

metab$fcon_ordinal

unique(dat$faeces_color)
metab$fcolor_ordinal <- 20
metab$fcolor_ordinal[metab$faeces_color=="light_brown"] <- 1
metab$fcolor_ordinal[metab$faeces_color=="brown"] <- 2
metab$fcolor_ordinal[metab$faeces_color=="dark_brown"] <- 3

unique(metab$fcolor_ordinal)

library(ggpmisc)
library(pdp)
library(patchwork)
library(cowplot)
library(corrplot)
library(caret)
library(microbiome)
library(ranger)

PMS.t <- transform(PMS_met, "compositional")
micro_ML<-as.data.frame(PMS.t@otu_table)
colnames(micro_ML) <- paste("ASV", seq(1, length(colnames(PMS.t@otu_table))), PMS.t@tax_table[,6], sep="_")

micro_ML$faecaes_consistency<-metab$fcon_ordinal

#corrplot(cor(micro_ML), type="upper", order="hclust")

set.seed(123)
train_index<-createDataPartition(micro_ML$faecaes_consistency, times=1, p=0.8, list=FALSE)
train_data <- micro_ML[train_index,]
test_data <- micro_ML[-train_index,]

tgrid <- expand.grid(
    mtry = 1:ncol(micro_ML-1),
    splitrule = "variance",
    min.node.size = c(5, 10))

#10-fold cross validation with 3 repeats
trainControl <- trainControl(method="repeatedcv", number=10, repeats=3, verboseIter = TRUE)
                                        #metric <- "Accuracy"

train_data$faecaes_consistency

set.seed(123)
fit.rf <- caret::train(faecaes_consistency~.,
                       data = train_data,
                       method="ranger",
                       tuneGrid=tgrid,
                       trControl = trainControl,
                       importance="permutation")



fit.rf

predict_treebag <- predict(fit.treebag, test <- data)
predict_rf <- predict(fit.rf, test_data)
print(postResample(predict_rf, test_data$faecaes_consistency))
pred_obs <- data.frame(obs=test_data$faecaes_consistency,
                           pred_rf= predict_rf)
cor.test(predict_rf, test_data$faecaes_consistency, method="spearman")

corrCon <-
    ggplot(pred_obs, aes(x=pred_rf, y=obs))+
    geom_point(size=3, color="orange")+
    geom_abline(linetype=5, color="blue", size=1)+
    stat_poly_line() +
    stat_poly_eq() +
    theme_classic()

## now for color
micro_ML$faecaes_consistency<-NULL
micro_ML$faecaes_color<-metab$fcolor_ordinal

set.seed(123)
train_index<-createDataPartition(micro_ML$faecaes_color, times=1, p=0.8, list=FALSE)
train_data <- micro_ML[train_index,]
test_data <- micro_ML[-train_index,]

set.seed(123)
fit.rf <- caret::train(faecaes_color~.,
                       data = train_data,
                       method="ranger",
                       tuneGrid=tgrid,
                       trControl = trainControl,
                       importance="permutation")



fit.rf

predict_rf <- predict(fit.rf, test_data)
print(postResample(predict_rf, test_data$faecaes_color))
pred_obs <- data.frame(obs=test_data$faecaes_color,
                           pred_rf= predict_rf)
cor.test(predict_rf, test_data$faecaes_color, method="spearman")

corrCol <-
    ggplot(pred_obs, aes(x=pred_rf, y=obs))+
    geom_point(size=3, color="orange")+
    geom_abline(linetype=5, color="blue", size=1)+
    stat_poly_line() +
    stat_poly_eq() +
    theme_classic()

corrCol



######################### NETWORK#############################################
##########################################################################
PMS_met

library(microbiome)

library(phyloseq)

## Filtering dataset
# prevalnce filter at 10%
KeepTaxap <- microbiome::prevalence(PMS_met)>0.10
PS.f <- phyloseq::prune_taxa(KeepTaxap, PMS_met)
#subset samples based on total read count (100 reads)
#PS.f <- phyloseq::prune_samples(sample_sums(PS.f)>100, PS.f)

met <- metabolite
asv <- (PS.f@otu_table)
tax <- data.frame(PS.f@tax_table)


colnames(asv) <- paste("cASV", seq(1:length(colnames(asv))), tax[,6], sep="_")
asv_met <- cbind(asv, met)


#spearman correlation and adjusting for multiple testing
otu.cor <- rcorr(as.matrix(asv_met), type="spearman")
otu.pval <- forceSymmetric(otu.cor$P) # Self-correlation as NA
cor.p <- p.adjust(otu.pval, method="BH")



# consider only significant and strong correlations
#cor.r1[which(!cor.r1 > 0.8 | !cor.r1 < -0.8)]=NA
#cor.p[which(cor.p>0.01)]=NA
otu.pval@x <- cor.p

    sel.tax <- tax[rownames(otu.pval),,drop=FALSE]
#sanity check
    all.equal(rownames(sel.tax), rownames(otu.pval))

p.yes <- otu.pval<0.01
r.val = otu.cor$r # select all the correlation values
p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion

## select asv based on rho
#p.yes.r <- abs(p.yes.r)>0.9 # output is logical vector
#p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.
adjm <- as.matrix(p.yes.r)

#colnames(adjm) <- as.vector(sel.tax$family)
#rownames(adjm) <- as.vector(sel.tax$family)

net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)



E(net.grph)$direction <- ifelse(edgew < 0,"red","blue")

E(net.grph)$weight <- abs(E(net.grph)$weight)

bad.vs<-V(net.grph)[degree(net.grph) == 0]
deg <- igraph::degree(net.grph, mode="all")

summary(deg)

net.grph <-delete.vertices(net.grph, bad.vs)

names2 <- gsub("cASV_.*_", "", names(V(net.grph)))

plot(net.grph,
#     vertex.size=deg/2,
     vertex.size=1,
     vertex.frame.color="black",
     edge.curved=F,
     edge.width=1.5,
     layout=layout.fruchterman.reingold,
     edge.color=E(net.grph)$direction,
#     vertex.label="",
     vertex.label.color="black",
     vertex.label.family="Times New Roman",
     vertex.label.font=2)
