## plotting
## unique animals

ID_df <- data.frame(ID=PMS@sam_data$hyena_ID, motherID=PMS@sam_data$genetic_mum, Clan=PMS@sam_data$Clan)
ID_df <- (unique(ID_df))

mother.df <- as.data.frame(table(mother.df$motherID))

mother_hist <-
    ggplot(mother.df, aes(x=Freq-1))+
    geom_histogram()+
    labs(x="Number of siblings", y="Individual count")+
    theme_classic()

Clan_df <- as.data.frame(table(ID_df$Clan))

clan_hist <- ggplot(Clan_df, aes(y=Freq, x=Var1, color=Var1))+
    geom_segment(data=Clan_df, aes(x=Var1, xend=Var1, y=0, yend=Freq),
                 color="lightgray", size=3)+
    geom_point(size=10)+
    scale_color_manual(values = c("#D16103", "#52854C", "#293352"))+
    geom_text(aes(label = Freq), color = "white", size = 3) +
    labs(x="Clan", y="Individual count")+
    guides(color="none")+
    theme(legend.position="none")+
    theme_classic()


Rank.df <- data.frame(Rank=PMS@sam_data$CSocialRank)

Rank_hist <-  
    ggplot(Rank.df, aes(x=Rank))+
    geom_histogram()+
    labs(x="Social rank", y="Sample count")+
    theme_classic()
#hist(PMS@sam_data$CSocialRank)

data.df <- data.frame(date_sampling=as.Date(PMS@sam_data$date_sampling,
                                format='%m/%d/%Y'))


Sampling_hist <- ggplot(data.df, aes(x=date_sampling))+
    geom_histogram()+
    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="Year of sampling", y="Sample count")+
    theme_classic()
Sampling_hist


FigureS1 <- plot_grid(mother_hist, clan_hist, Rank_hist, Sampling_hist, labels="auto")
ggplot2::ggsave(file="Figures/FigureS1.pdf", FigureS1, width = 190, height = 150, dpi = 300, units="mm")


test <- PMS@sam_data %>%
    group_by(hyena_ID) %>%
    mutate(min_age=min(age_sampling),
           n= n(),
           max_age=max(age_sampling))

PMS@sam_data$TimeP <- 0
PMS@sam_data$TimeP[test$age_sampling==test$min_age] <- 1
PMS@sam_data$TimeP[test$age_sampling>test$min_age] <- 3
PMS@sam_data$TimeP[test$age_sampling==test$max_age] <- 2

repeatedS <- ggplot(data=sample_data(PMS),
                    aes(x=age_sampling/365,
                        y=reorder(hyena_ID, age_sampling),
                        group=hyena_ID,
                        fill=age_sampling_cat))+
    geom_point(colour="white",pch=21,size=3, alpha=0.7)+
    geom_line(size=0.5,alpha=0.5)+
    geom_vline(xintercept=2, colour="firebrick", size=2, alpha=0.5)+
    ylab("Hyena ID")+
    xlab("Age at sampling (years)")+
    scale_fill_manual(values=c("#a5c3ab", "#eac161"))+
    scale_x_continuous(breaks=0:16)+
    theme_classic()+
    theme(legend.position="none",
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          axis.line=element_line(colour="black"))+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())


IgA_hist <- ggplot(sample_data(PMS), aes(x=log(IgA)))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="Faecal IgA (RU, log scale)", y="Sample count")+
    theme_classic()

mucin_hist <- ggplot(sample_data(PMS), aes(x=log(mucin)))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="Faecal mucin (OE, log scale)", y="Sample count")+
    theme_classic()



Age_diff <- ggplot(data.dyad, aes(x=Age))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="Age distance", y="Pair count")+
    theme_classic()

IgA_diff <- ggplot(data.dyad, aes(x=IgAP))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="f-IgA distance", y="Pair count")+
    theme_classic()

Muc_diff <- ggplot(data.dyad, aes(x=MucinP))+
    geom_histogram()+
#    scale_x_date(date_breaks="year", date_labels="%Y")+
    labs(x="f-mucin distance", y="Pair count")+
    theme_classic()

AgeSD <- plot_grid(repeatedS, Age_diff, nrow=1, rel_widths=c(1, 0.5), labels=c("a", "b"))

ImmSD <- plot_grid(IgA_hist, mucin_hist, IgA_diff, Muc_diff, labels=c("c", "d", "e", "f"))

FigureSD <- plot_grid(AgeSD, ImmSD, nrow=2)

ggplot2::ggsave(file="Figures/FigureSD.pdf", FigureSD, width = 190, height = 190, dpi = 300, units="mm")

plot_grid(Age_diff, IgA_diff, Muc_diff, nrow=1, labels=c("d", "e", "f"))

hist(data.dyad$AgeDiff)


##############
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

mds.fit <- envfit(mds~fGCM+IgA+IgG+mucin+age_sampling+ CSocialRank+ season+ sex + Survival_Ad + hyena_ID + age_death + clan + LRS, data=df, na.rm=TRUE)

data.scores <- as.data.frame(vegan::scores(mds)$sites)

all(rownames(data.scores)==rownames(df))
data.scores <- cbind(data.scores, df[rownames(data.scores),c("fGCM", "IgA", "mucin", "age_sampling", "CSocialRank", "season", "sex", "Survival_Ad", "hyena_ID", "age_death", "clan", "LRS")])

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


######################################################
### Individual repeatability
library(dplyr)
samplecount <- data.frame(sample_data(PMS)) %>%
    group_by(hyena_ID) %>%
    summarise(NoSamples=n())

sample_data(PMS)$Scount <- "single"
sample_data(PMS)$Scount[sample_data(PMS)$hyena_ID%in%samplecount$hyena_ID[samplecount$NoSamples>1]] <- "rep"


PMS.r <- subset_samples(PMS, Scount=="rep")
PMS.r@sam_data$age_sampling
PMS.r <- prune_taxa(taxa_sums(PMS.r)>0, PMS.r)
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


#temp <- as.data.frame(sample_data(PMS) %>% group_by(hyena_ID) %>%
                              sample_n(1)) #subset phyloseq object to one randomly selected sample 
#one.ps <- subset_samples(PMS, Sample%in%temp$Sample)
#ps.core <- core(one.ps, detection = 0, prevalence = .3)
#Core<- phyloseq::prune_taxa(taxa_names(PMS)%in%taxa_names(ps.core), PMS)
Core <- core(PMS.r, detection = 0, prevalence = .3)


# shannon diversity repeatability
library(rptR)
shannon <- vegan::diversity(Core@otu_table, index = "shannon")  
Core@sam_data$shannon <- shannon
metadf <- as.data.frame(Core@sam_data)
class(metadf) <- "data.frame"
rpt(shannon~age_sampling+(1|hyena_ID), data=metadf, datatype = "Gaussian", nboot = 1000, npermut = 1000, grname=c("hyena_ID"))

repeatedS <- ggplot(data=sample_data(PMS.r),
                    aes(x=age_sampling,
                        y=reorder(hyena_ID, age_sampling),
                        group=hyena_ID,
                        fill=age_sampling_cat))+
    geom_line(size=0.5,alpha=0.5)+
    geom_point(colour="black",pch=21,size=3)+
    geom_vline(xintercept=730, colour="firebrick", size=2, alpha=0.5)+
    ylab("Hyena ID")+
    xlab("Age at sampling")+
    scale_fill_manual(values=c("#a5c3ab", "#eac161"))+
    theme_classic()+
    theme(legend.position="none",
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          axis.line=element_line(colour="black"))+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

library("GUniFrac", lib="/usr/local/lib/R/site-library")


jac <- vegdist(Core@otu_table, method="jaccard")
ait <- vegdist(Core@otu_table, method="aitchison",
               pseudocount=1)


dICC(jac, strata=Core@sam_data$hyena_ID)
dICC(ait, strata=Core@sam_data$hyena_ID)

permajac <- adonis2(jac~
                    Core@sam_data$age_sampling +
                    Core@sam_data$hyena_ID,
                    by="margin")
permajac

permaait <- adonis2(ait~
                    Core@sam_data$age_sampling+
                    Core@sam_data$hyena_ID,
                    by="margin")
permaait

otu <- otu_table(Core)
colnames(otu) <- as.vector(tax_table(Core)[,"Genus"])
mds <- vegan::metaMDS(otu, distance="jaccard", k=2)
#plot(mds)

stressplot(mds)

data.scores <- vegan::scores(mds)
data.scores <- as.data.frame(data.scores$sites)
data.scores$hyena_ID <- Core@sam_data$hyena_ID
data.scores$age_sampling <- Core@sam_data$age_sampling
data.scores$age_sampling_cat <- Core@sam_data$age_sampling_cat
head(data.scores)

Jac_R <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2, fill=age_sampling_cat))+
    geom_point(size=3, shape=21)+
    geom_line(aes(group=hyena_ID), colour="gray", alpha=0.8)+
    scale_fill_manual(values=c("#a5c3ab", "#eac161"), labels=c("Adult", "Juvenile"))+
    labs(fill="Age at sampling")+
    theme_classic()+
    theme(legend.position="none")

### aitchison
mds2 <- vegan::metaMDS(otu, distance="aitchison", pseudocount=1)
#plot(mds2)

data.scores2 <- vegan::scores(mds2)
data.scores2 <- as.data.frame(data.scores2$sites)
data.scores2$hyena_ID <- Core@sam_data$hyena_ID
data.scores2$age_sampling <- Core@sam_data$age_sampling
data.scores2$age_sampling_cat <- Core@sam_data$age_sampling_cat
head(data.scores)

Ait_R <- ggplot(data.scores2, aes(x=NMDS1, y=NMDS2, fill=age_sampling_cat))+
    geom_point(size=3, shape=21)+
    geom_line(aes(group=hyena_ID), colour="gray", alpha=0.8)+
    scale_fill_manual(values=c("#a5c3ab", "#eac161"), labels=c("Adult", "Juvenile"))+
    labs(fill="Age at sampling")+
    theme_classic()+
    theme(legend.position="none")


#plot_grid(Jac_R, Ait_R)

############# let's do dyad comparisons now
sample_data(Core)$key <- paste(sample_data(Core)$hyena_ID,
                              sample_data(Core)$Sample, sep="_")
key <- data.frame(ID=sample_data(Core)$key)
mt <- sample_data(Core)
# renaming this now
sample_names(Core) <- key$ID

####################
## 1) Jaccard distance
JM <- as.matrix(phyloseq::distance(Core,
                                     method="jaccard",
                                     type="samples",
                                     binary=T))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JM <- 1-JM
ja<-c(as.dist(JM))
## 2) Aitchison distance
AM <- as.matrix(vegan::vegdist(Core@otu_table,
                                 method="aitchison",
                                 pseudocount=1))
AM <- 1-AM
at<-c(as.dist(AM))

#Combine these vectors into a data frame
data.d<-data.frame(MS_J=ja,
                      MS_A=at)


list<-expand.grid(key$ID, key$ID)
list<-list[which(list$Var1!=list$Var2),]
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))
i=nrow(key)
JM[which(rownames(JM)==list$Var1[i]),which(colnames(JM)==list$Var2[i])]==ja[i]

# add the names of both individuals participating in each dyad into the data frame
data.d$IDA_s<-list$Var2
data.d$IDB_s<-list$Var1
data.d$IDA <- gsub("_.*", "", data.d$IDA_s)
data.d$IDB <- gsub("_.*", "", data.d$IDB_s)

data.d$Var <- "ABC"
data.d$Var[which(data.d$IDA!=data.d$IDB)] <- "Inter"
data.d$Var[which(data.d$IDA==data.d$IDB)] <- "Intra"

res <- (aov(data=data.d, MS_J~Var))
summary(res)

resA <- (aov(data=data.d, MS_A~Var))
summary(resA)

library(viridis)

B <- ggplot(data.dyad, aes(x=Var, y=MS_J, fill=Var))+
    geom_violin(alpha=0.5, width=1)+
     geom_boxplot(outlier.shape=NA, alpha=0.5, width=0.1)+
#    stat_summary(fun.data=mean_sdl, mult=1,
#                             geom="pointrange", color="red")+
#    geom_jitter(position=position_jitter(seed=1, width=0.2), alpha=0.1)+
    labs(x="Dyadic comparisons", y="Jaccard similarities")+
    scale_fill_viridis(discrete = TRUE) +
    scale_x_discrete(labels=c("Inter"="Inter-individual", "Intra"="Intra-individual"))+
    theme_classic()+
    theme(legend.position="none")

C <- ggplot(data.dyad, aes(x=Var, y=MS_A, fill=Var))+
        geom_violin(alpha=0.5, width=1)+
     geom_boxplot(outlier.shape=NA, alpha=0.5, width=0.1)+
#    stat_summary(fun.data=mean_sdl, mult=1,
#                             geom="pointrange", color="red")+
#    geom_jitter(position=position_jitter(seed=1, width=0.2), alpha=0.1)+
    labs(x="Dyadic comparisons", y="Aitchison similarities")+
    scale_fill_viridis(discrete = TRUE) +
    scale_x_discrete(labels=c("Inter"="Inter-individual", "Intra"="Intra-individual"))+
    theme_classic()+
    theme(legend.position="none")

dis_R <- plot_grid(B, C, labels=c("b", "c"))
Fig1 <- plot_grid(repeatedS, dis_R, ncol=1, labels=c("a", ""))

ggsave("Figures/Figure1.pdf", Fig1, width = 170, height = 170, dpi = 300, units="mm") 

