library(ggplot2)
library(phyloseq)
library(vegan)
#library(tidyverse)



PMS <- readRDS("/SAN/Susanas_den/gitProj/AA_Hyenas_Pakt/tmp/fPMS.rds")
PMS <- prune_taxa(taxa_sums(PMS)>0, PMS)

sample_data(PMS)$AFR[sample_data(PMS)$Survival_Ad=="No"] # this makes sense

#sample_data(PMS)$LRS[sample_data(PMS)$Survival_Ad=="No"] <- 0 # gonna set this to zero, since dead animals don't reproduce
#sample_data(PMS)$LRS # still quite some NAs, why?

#let's input what we can
library(mice)
pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(PMS@sam_data,2,pMiss)

# to input: CSocialRank;
df <- data.frame(Sample=PMS@sam_data$Sample, SocialRank=PMS@sam_data$CSocialRank)
socialR <- mice(df, m=5, maxit=50, meth="pmm", seed=500)
socialR <- complete(socialR, 1)
PMS@sam_data$CSocialRank_inputed <- socialR$SocialRank

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

#defining Clan
PMS@sam_data$Clan <- substr(PMS@sam_data$hyena_ID, 1,1)

jac <- vegdist(PMS@otu_table, method="jaccard")
bray <- vegdist(PMS@otu_table, method="bray")

permajac <- adonis2(jac~
                    PMS@sam_data$CSocialRank_inputed+
                    PMS@sam_data$age_sampling+
                    PMS@sam_data$sex+
                    PMS@sam_data$IgA_inputed+
                    PMS@sam_data$mucin_inputed+
                    PMS@sam_data$Clan+
                    PMS@sam_data$Survival_Ad+
                    PMS@sam_data$season,
                    strata=PMS@sam_data$hyena_ID,
                    by="margin")

permabray <- adonis2(bray~
                    PMS@sam_data$CSocialRank_inputed+
                    PMS@sam_data$Clan+
                    PMS@sam_data$IgA_inputed+
                    PMS@sam_data$mucin_inputed+
                    PMS@sam_data$Survival_Ad+
                    PMS@sam_data$age_sampling+
                    PMS@sam_data$sex+
                    PMS@sam_data$season,
                    strata=PMS@sam_data$hyena_ID,
                    by="margin")


permajac
permabray

otu <- otu_table(PMS)
colnames(otu) <- as.vector(tax_table(PMS)[,"Genus"])
mds <- vegan::metaMDS(otu, try=350, k=3)
#plot(mds)

df <- PMS@sam_data
class(df) <- "data.frame"
cplot(mds)
names(df)

ordiplot(mds, type="n")
orditorp(mds, display="species", col="red", air=0.01)
orditorp(mds, display="sites")

class(mds)


mds.fit <- envfit(mds~fGCM+IgA+IgG+mucin+age_sampling+ CSocialRank_inputed+ season+ sex + Survival_Ad + hyena_ID + age_death + clan + LRS, data=df, na.rm=TRUE)

data.scores <- as.data.frame(vegan::scores(mds)$sites)

all(rownames(data.scores)==rownames(df))
data.scores <- cbind(data.scores, df[rownames(data.scores),c("fGCM", "IgA", "mucin", "age_sampling", "CSocialRank_inputed", "season", "sex", "Survival_Ad", "hyena_ID", "age_death", "clan", "LRS")])

mds.fit$vectors



EnCoord <- as.data.frame(rbind(vegan::scores(mds.fit, "vectors")*
      ordiArrowMul(mds.fit),
      vegan::scores(mds.fit, "factors")*
      ordiArrowMul(mds.fit)))


EnCoord$Cat <- ifelse(rownames(EnCoord)%in%c("fGCM", "IgA", "mucin", "age_death", "LRS"), "plot", "notplot")

NMDS_env <-ggplot(data.scores, aes(x=NMDS1, y=NMDS2))+
        geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                    data = subset(EnCoord, Cat=="plot"), size =1, alpha = 0.5) +
    geom_point(data=data.scores, shape=21, aes(fill=age_sampling), size=3, alpha=0.5)+
    geom_point(data = subset(EnCoord, Cat=="plot"), aes(x = NMDS1, y = NMDS2),
                  shape = 4, size = 4, alpha = 0.6)+
        geom_text(data = subset(EnCoord, Cat=="plot"), aes(x = NMDS1, y = NMDS2+0.04),
             label = rownames(subset(EnCoord, Cat=="plot")), colour = "red", fontface = "bold")+
    scale_fill_gradient2(midpoint=365, low="blue", mid="lightblue", high="red", space ="Lab" )+
    labs(x = "NMDS axis 1", y = "NMDS axis 2", fill = "Age at sampling") +
    theme_bw(base_size=12)+
    theme(legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          legend.key.size = unit(1.2, "lines"),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"))

NMDS_env

ggsave("Figures/NMDS_env.pdf", NMDS_env, width=150, height=150, units="mm", dpi=300)

## permanova on juveniles
## partial least squares regression
library(caret)
library(pls)

Juv <-df[df$age_sampling_cat=="sampled_as_juvenile",]
JPS <- subset_samples(PMS, age_sampling_cat=="sampled_as_juvenile")
Juv <- Juv[, colnames(Juv)%in%c("age_sampling", "CSocialRank_inputed", "season", "sex", "Survival_Ad", "hyena_ID", "clan")]
JuvM <- data.frame(otu_table(JPS), Juv)
JuvM$season <- as.factor(Juv$season)
JuvM$sex <- as.factor(Juv$sex)
JuvM$Survival_Ad <- as.factor(Juv$Survival_Ad)
JuvM$clan <- as.factor(Juv$clan)
JuvM$hyena_ID <- as.factor(Juv$hyena_ID)

PMS@sam_data$season <- as.factor(PMS@sam_data$season)
PMS@sam_data$sex <- as.factor(PMS@sam_data$sex)
PMS@sam_data$Survival_Ad <- as.factor(PMS@sam_data$Survival_Ad)
PMS@sam_data$clan <- as.factor(PMS@sam_data$clan)
PMS@sam_data$hyena_ID <- as.factor(PMS@sam_data$hyena_ID)


str(sample_data(PMS)$hyena_ID)

set.seed(123)
training.samples <- sample(sample_names(PMS), size=0.8*nrow(PMS@sam_data))

train.dataPS <- subset_samples(PMS, sample_names(PMS)%in%training.samples)
test.dataPS <- subset_samples(PMS, !sample_names(PMS)%in%training.samples)

levels(PMS@sam_data$hyena_ID)

train.df <- data.frame(otu_table(train.dataPS), train.dataPS@sam_data[, c("age_sampling", "clan", "sex", "season", "CSocialRank_inputed", "Survival_Ad")])

test.df <- data.frame(otu_table(test.dataPS), test.dataPS@sam_data[, c("age_sampling", "clan", "sex", "season", "CSocialRank_inputed", "Survival_Ad")])

set.seed(123)
model <- train(
    age_sampling~., data=train.df, method="pls")
#    scale=TRUE,
#    trControl=trainControl("cv", number=10),
#   tuneLength=10)

model$bestTune

plot(model)

model

tax <- as.data.frame(tax_table(PMS), stringsAsFactors=FALSE)
model$coefnames <- c(tax$Phylum, model$coefnames[1082:1087])
imp <- (varImp(model))
imp <- imp$importance
imp$names <- c(tax$Phylum, model$coefnames[1082:1087])
ggplot(imp, aes(Overall, names))+
    geom_point()

summary(model$finalModel)

# Make predictions
predictions <- model %>% predict(test.df)
predictions
sqrt(mean((test.df$age_sampling - predictions)^2)) # RMSE
cor(test.df$age_sampling, predictions) ^ 2 # R2
names(test.df)


### Individual repeatability

samplecount <- data.frame(sample_data(PMS)) %>%
    group_by(hyena_ID) %>%
    summarise(NoSamples=n())
sample_data(PMS)$Scount <- "single"
sample_data(PMS)$Scount[sample_data(PMS)$hyena_ID%in%samplecount$hyena_ID[samplecount$NoSamples>1]] <- "rep"

PMS.r <- subset_samples(PMS, Scount=="rep")
PMS.r@sam_data$age_sampling
PMS.r <- prune_taxa(taxa_sums(PMS.r)>0, PMS.r)
PMS.r



library(dplyr)

PMS.r

test <- PMS.r@sam_data %>%
    group_by(hyena_ID) %>%
    mutate(min_age=min(age_sampling),
           n= n(),
           max_age=max(age_sampling))

all(test$Sample==PMS.r@sam_data$Sample)

PMS.r@sam_data$TimeP <- 0
PMS.r@sam_data$TimeP[test$age_sampling==test$min_age] <- 1
PMS.r@sam_data$TimeP[test$age_sampling>test$min_age] <- 3
PMS.r@sam_data$TimeP[test$age_sampling==test$max_age] <- 2

### core analysis
library(microbiome)
taxa_names(PMS.r) <- paste("ASV", seq(1:length(taxa_names(PMS.r))), PMS.r@tax_table[,6], sep="_")

Core<-microbiome::core(PMS.r, detection = 0, prevalence = 0.30)

plot_core(Core,
             plot.type = "lineplot") +
    xlab("Relative Abundance (%)") +
      theme_bw()

library(viridis)

plot_core(Core_asv,
              plot.type = "heatmap")+
        scale_fill_viridis()

c
# shannon diversity repeatability
library(rptR)
shannon <- vegan::diversity(PMS.r@otu_table, index = "shannon")  
PMS.r@sam_data$shannon <- shannon
metadf <- as.data.frame(PMS.r@sam_data)
class(metadf) <- "data.frame"
rpt(shannon~(1|hyena_ID), data=metadf, datatype = "Gaussian", nboot = 100, npermut = 100, grname="hyena_ID")
# for core only
shannon <- vegan::diversity(Core@otu_table, index = "shannon")  
Core@sam_data$shannon <- shannon
metadf <- as.data.frame(Core@sam_data)
class(metadf) <- "data.frame"
rpt(shannon~(1|hyena_ID), data=metadf, datatype = "Gaussian", nboot = 100, npermut = 100, grname="hyena_ID")


PMS.r@sam_data$age_sampling

install.packages("see")
library(see)

ggplot(data=Core@sam_data, aes(x=as.factor(TimeP), y=shannon))+
    geom_point(size=2, alpha=0.5, fill="gray")+
     geom_line(aes(group=hyena_ID), colour="gray")+
    theme_bw(base_size=14)+
    ylab("Shannon diversity")+
    xlab("Sampling timepoint")+
    theme(panel.grid.minor=element_blank(),
          panel.background=element_blank(),axis.line=element_line(colour="black"))
 


ggplot(data=sample_data(PMS.r),aes(x=age_sampling,y=reorder(hyena_ID, age_sampling),group=hyena_ID))+
    geom_line(size=0.5,alpha=0.5)+
    geom_point(colour="black",pch=21,size=3, fill = "grey")+
    geom_vline(xintercept=730, colour="firebrick", size=1, alpha=0.2)+
    ylab("Hyena ID")+
    xlab("Age at sampling")+
    theme_bw(base_size=14)+
    theme(panel.grid.minor=element_blank(),
          panel.background=element_blank(),axis.line=element_line(colour="black"))+
    theme(axis.text.y=element_blank(),
                      axis.ticks.y=element_blank())


Juv <- subset_samples(PMS.r, age_sampling<730)
metadf <- as.data.frame(Juv@sam_data)
class(metadf) <- "data.frame"
rpt(shannon~(1|hyena_ID), data=metadf, datatype = "Gaussian", nboot = 100, npermut = 100, grname="hyena_ID")

JCore<-microbiome::core(Juv, detection = 0, prevalence = 0.30)
metadf <- as.data.frame(JCore@sam_data)
class(metadf) <- "data.frame"
rpt(shannon~(1|hyena_ID), data=metadf, datatype = "Gaussian", nboot = 100, npermut = 100, grname="hyena_ID")

Juv

asv_names<-taxa_names(scaled_core_asv)

bray <- PMS.r

jac <- vegdist(PMS.r@otu_table, method="jaccard")
bray <- vegdist(PMS.r@otu_table, method="bray")

permajac <- adonis2(jac~
                    PMS.r@sam_data$CSocialRank_inputed+
                    PMS.r@sam_data$age_sampling+
                    PMS.r@sam_data$hyena_ID,
                    by="margin")
permajac
permabray <- adonis2(bray~
                    PMS.r@sam_data$CSocialRank_inputed+
                    PMS.r@sam_data$age_sampling+
                    PMS.r@sam_data$hyena_ID,
                    by="margin")
permabray

# core
jac <- vegdist(Core@otu_table, method="jaccard")
bray <- vegdist(Core@otu_table, method="bray")

permajac <- adonis2(jac~
                    PMS.r@sam_data$CSocialRank_inputed+
                    PMS.r@sam_data$age_sampling+
                    PMS.r@sam_data$hyena_ID,
                    by="margin")

permabray <- adonis2(bray~
                    PMS.r@sam_data$CSocialRank_inputed+
                    PMS.r@sam_data$age_sampling+
                    PMS.r@sam_data$hyena_ID,
                    by="margin")

permajac
permabray

otu <- otu_table(PMS.r)
colnames(otu) <- as.vector(tax_table(PMS.r)[,"Genus"])
mds <- vegan::metaMDS(otu)
plot(mds)

df <- PMS.r@sam_data
class(df) <- "data.frame"
cplot(mds)
names(df)

data.scores <- vegan::scores(mds)

data.scores <- as.data.frame(data.scores$sites)

all(rownames(data.scores)==PMS.r@sam_data$Sample)

data.scores$hyena_ID <- PMS.r@sam_data$hyena_ID
data.scores$age_sampling <- PMS.r@sam_data$age_sampling
data.scores$age_sampling_cat <- PMS.r@sam_data$age_sampling_cat
head(data.scores)

ggplot(data.scores, aes(x=NMDS1, y=NMDS2, fill=age_sampling_cat))+
    geom_point(size=3, shape=21)+
    geom_line(aes(group=hyena_ID), colour="gray", alpha=0.8)+
    scale_fill_manual(values=c("#a5c3ab", "#eac161"))+
     theme_bw(base_size=14)+
    theme(panel.grid.minor=element_blank(),
          panel.background=element_blank(),axis.line=element_line(colour="black"))+
    theme(axis.text.y=element_blank(),
          axis.text.x=element_blank())


(scores(mds))

############# let's do dyad comparisons now
sample_data(PMS)$key <- paste(sample_data(PMS)$hyena_ID, sample_data(PMS)$Sample, sep="_")
key <- data.frame(ID=sample_data(PMS)$key)
metadt <- sample_data(PMS)
# renaming this now
sample_names(PMS) <- key$ID

####################
## 1) Jaccard distance
JACM <- as.matrix(phyloseq::distance(PMS, method="jaccard", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
dimnames(JACM)<- c(key, key)

## 2) Chisq distance
CHIM <- as.matrix(vegan::vegdist(PMS@otu_table, method="chisq"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
CHIM <- 1-CHIM
# sanity check
all(rownames(CHIM)==key)
dimnames(CHIM)<- c(key, key)

## 3) Bray distance
BRAYM <- as.matrix(phyloseq::distance(PMS, method="bray", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
BRAYM <- 1-BRAYM
# sanity check
all(rownames(BRAYM)==key)
dimnames(BRAYM)<- c(key, key)

## 4) Sex pairs
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
dimnames(SEXM)<-c(key, key)

# 5) Making age distances
#Create data frame with each sample name (character) and sampling time (numeric)
AGE_frame<-metadt[,c("key", "age_sampling")]
#Create an empty matrix to fill with distances
AGEM<-array(0,c(nrow(AGE_frame),nrow(AGE_frame)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(AGE_frame)){
    for (j in 1:nrow(AGE_frame))
    {AGEM[i,j]=abs(AGE_frame$age_sampling[i] -AGE_frame$age_sampling[j])
    }
}
dimnames(AGEM) <- c(key, key)


# 6) Making rank distances
RANK_frame<-metadt[,c("key", "CSocialRank_inputed")]
RANKM<-array(0,c(nrow(RANK_frame),nrow(RANK_frame)))
for (i in 1:nrow(RANK_frame)){
    for (j in 1:nrow(RANK_frame))
    {RANKM[i,j]=abs(RANK_frame$CSocialRank_inputed[i] -RANK_frame$CSocialRank_inputed[j])
    }
}
dimnames(RANKM) <- c(key, key)

# 7) Making IgA distances
IgA_frame<-metadt[,c("key", "IgA_inputed")]
IGAM<-array(0,c(nrow(IgA_frame),nrow(IgA_frame)))
for (i in 1:nrow(IgA_frame)){
    for (j in 1:nrow(IgA_frame))
    {IGAM[i,j]=abs(IgA_frame$IgA_inputed[i] -IgA_frame$IgA_inputed[j])
    }
}
dimnames(IGAM) <- c(key, key)

# 8) Making IgA distances
Mucin_frame<-metadt[,c("key", "mucin_inputed")]
MUCM<-array(0,c(nrow(Mucin_frame),nrow(Mucin_frame)))
for (i in 1:nrow(Mucin_frame)){
    for (j in 1:nrow(Mucin_frame))
    {MUCM[i,j]=abs(Mucin_frame$mucin_inputed[i] -Mucin_frame$mucin_inputed[j])
    }
}
dimnames(MUCM) <- c(key, key)

chi<-c(as.dist(CHIM))
bray <- c(as.dist(BRAYM))
jac<-c(as.dist(JACM))
age<-c(as.dist(AGEM))
sex<-c(SEXM[lower.tri(SEXM)])
iga<-c(as.dist(IGAM))
muc<-c(as.dist(MUCM))
rank <- c(as.dist(RANKM))

#Combine these vectors into a data frame
data.dyad<-data.frame(Microbiomejac=jac,Microbiomechi=chi, Microbiomebray=bray, Age=age, Sex=sex, IgA=iga, Mucin=muc, Rank=rank)

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
data.dyad$IDA<-list$Var2
data.dyad$IDB<-list$Var1
# Make sure you have got rid of all self comparisons
data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]

######################### Now we model the data ####################
#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("Rank","Age", "IgA", "Mucine")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }


library(brms)

#### model
modelJ<-brm(Microbiomejac~1+ Age+Sex+IgA+Mucin+Rank+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)

modelJ
