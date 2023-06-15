library(ggplot2)
library(phyloseq)
library(vegan)

PMS <- readRDS("/SAN/Susanas_den/gitProj/AA_Hyenas_Pakt/tmp/fPMS.rds")

head(sample_data(PMS))

sample_data(PMS)$Survival_Ad

sample_data(PMS)$AFR[sample_data(PMS)$Survival_Ad=="No"] # this makes sense

sample_data(PMS)$LRS[sample_data(PMS)$Survival_Ad=="No"] <- 0 # gonna set this to zero, since dead animals don't reproduce

sample_data(PMS)$LRS # still quite some NAs, why?

#let's input what we can
library(mice)

PMS@sam_data$Sample

pMiss <- function(x){sum(is.na(x))/length(x)*100}

apply(PMS@sam_data,2,pMiss)

                                        # to input: CSocialRank;
df <- data.frame(Sample=PMS@sam_data$Sample, SocialRank=PMS@sam_data$CSocialRank)
socialR <- mice(df, m=5, maxit=50, meth="pmm", seed=500)
socialR <- complete(socialR, 1)
PMS@sam_data$CSocialRank_inputed <- socialR$SocialRank

jac <- vegdist(PMS@otu_table, method="jaccard")

bray <- vegdist(PMS@otu_table, method="bray")

head(PMS@sam_data)

PMS@sam_data$Clan <- substr(PMS@sam_data$hyena_ID, 1,1)

## I think we need clan

summary(as.factor(PMS@sam_data$hyena_ID))


permajac <- adonis2(jac~
                        PMS@sam_data$CSocialRank_inputed+
                    PMS@sam_data$age_sampling+
                    PMS@sam_data$sex+
                    PMS@sam_data$Clan+
                    PMS@sam_data$Survival_Ad+
                    PMS@sam_data$season,
                    strata=PMS@sam_data$hyena_ID,
                    by="margin")

permabray <- adonis2(bray~
                    PMS@sam_data$CSocialRank_inputed+
                    PMS@sam_data$Clan+
                    PMS@sam_data$Survival_Ad+
                    PMS@sam_data$age_sampling+
                    PMS@sam_data$sex+
                    PMS@sam_data$season,
                    strata=PMS@sam_data$hyena_ID,
                    by="margin")



permajac
permabray

PMS

PMS <- prune_taxa(taxa_sums(PMS)>0, PMS)

otu <- otu_table(PMS)

colnames(otu) <- as.vector(tax_table(PMS)[,"Genus"])

mds <- vegan::metaMDS(otu, try=350, k=3)

plot(mds)


df <- PMS@sam_data

class(df) <- "data.frame"

cplot(mds)

names(df)

ordiplot(mds, type="n")
orditorp(mds, display="species", col="red", air=0.01)
orditorp(mds, display="sites")




mds.fit <- envfit(mds~fGCM+IgA+IgG+mucin+age_sampling+ CSocialRank_inputed+ season+ sex + Survival_Ad + hyena_ID + age_death + clan + LRS, data=df, na.rm=TRUE)

df$LRS

data.scores <- as.data.frame(scores(mds)$sites)

scores(mds)$species

all(rownames(data.scores)==rownames(df))

data.scores <- cbind(data.scores, df[rownames(data.scores),c("fGCM", "IgA", "mucin", "age_sampling", "CSocialRank_inputed", "season", "sex", "Survival_Ad", "hyena_ID", "age_death", "clan", "LRS")])

head(data.scores)

EnCoord <- as.data.frame(rbind(scores(mds.fit, "vectors")*
      ordiArrowMul(mds.fit),
      scores(mds.fit, "factors")*
      ordiArrowMul(mds.fit)))


EnCoord$Cat <- ifelse(rownames(EnCoord)%in%c("fGCM", "IgA", "mucin", "age_death", "LRS"), "plot", "notplot")

ggplot(data.scores, aes(x=NMDS1, y=NMDS2))+
    geom_point(data=data.scores, aes(colour=age_sampling), size=3, alpha=0.5)+
    geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                    data = subset(EnCoord, Cat=="plot"), size =1, alpha = 0.5) +
    geom_point(data = subset(EnCoord, Cat=="plot"), aes(x = NMDS1, y = NMDS2),
                  shape = "diamond", size = 4, alpha = 0.6)+
        geom_text(data = subset(EnCoord, Cat=="plot"), aes(x = NMDS1, y = NMDS2+0.04),
                                  label = rownames(subset(EnCoord, Cat=="plot")), colour = "red", fontface = "bold")


df$age_sampling_cat

EnCoord

subset(EnCoord$Cat=="plot")


## permanova on juveniles



head(Juv)

apply(Juv,2,pMiss)

Juv$Survival_Ad



str(data.scores)


## partial least squares regression

library(caret)
library(pls)

#Juv <-df[df$age_sampling_cat=="sampled_as_juvenile",]
#JPS <- subset_samples(PMS, age_sampling_cat=="sampled_as_juvenile")
#Juv <- Juv[, colnames(Juv)%in%c("age_sampling", "CSocialRank_inputed", "season", "sex", "Survival_Ad", "hyena_ID", "clan")]
#JuvM <- data.frame(otu_table(JPS), Juv)
#JuvM$season <- as.factor(Juv$season)
#JuvM$sex <- as.factor(Juv$sex)
#JuvM$Survival_Ad <- as.factor(Juv$Survival_Ad)
#JuvM$clan <- as.factor(Juv$clan)
#JuvM$hyena_ID <- as.factor(Juv$hyena_ID)

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

#  data partitioning strategy like k-fold cross-validation that resamples and splits our data many times. We then train the model on these samples and pick the best model. Caret makes this easy with the trainControl method.

set.seed(1)
ctrl <- trainControl(
    method = "cv",
    number = 10,
)


set.seed(123)
model1 <- train(
    age_sampling~., data=train.df, method="pls",
#    scale=TRUE,
    trControl=ctrl,
   tuneLength=10)

model1


predictions <- model1 %>% predict(test.df)
sqrt(mean((test.df$age_sampling - predictions)^2)) # RMSE
cor(test.df$age_sampling, predictions) ^ 2 # R2

plot(model1, plottype="coef")

model2 <- plsr(
    age_sampling~., data=pms.df,
#    scale=TRUE,
    vaidation="LOO",
    ncomp=10)

summary(model2)

plot(model2, ncomp = 2, asp = 1, line = TRUE)

scoreplot(model2, comps=1:2)

# Model performance metrics
data.frame(
    RMSE = caret::RMSE(predictions, test.df$age_sampling),
    Rsquare = caret::R2(predictions, test.df$age_sampling)
    )

pls_biplot <- list("loadings"=loadings(model$finalModel),
                   "scores"=scores(model$finalModel))

class(pls_biplot$scores)  <- "matrix"

sample_names(PMS)%in%rownames(na.omit())


pls_biplot$scores <- data.frame(sample_data(JPS.s), pls_biplot$scores)

tax <- as.data.frame(tax_table(JPS.s), stringsAsFactors=FALSE)

class(pls_biplot$loadings) <- "matrix"

pls_biplot$loadings <- data.frame(c(tax$Phylum, rownames(pls_biplot$loadings)[1082:1212])
, pls_biplot$loadings)

names(pls_biplot$loadings) <- c("Phylum", "Comp.1", "Comp.2")

pls_biplot$loadings$Phylum <- gsub("hyena.*", "", pls_biplot$loadings$Phylum)

names(pls_biplot$scores)

class(pls_biplot$loadings)

summary(pls_biplot$loadings$Comp.1)

summary(pls_biplot$scores$Comp.1)

ggplot()+
    geom_point(data=pls_biplot$loadings,
               aes(x=Comp.1*4000, y=Comp.2*4000, col=Phylum),
               size=2)+
    geom_point(data=pls_biplot$scores,
               aes(x=Comp.1, y=Comp.2, shape=Survival_Ad), size=3)+
#    scale_color_brewer(palette = "Set1")+
        scale_shape_discrete(solid=FALSE) +
            labs(x = "Axis1", y = "Axis2", col = "Phylum") +
            guides(col = guide_legend(override.aes = list(size = 3))) +
                theme_bw()

plot(model)

summary(model$finalModel)


predictions <- model %>% predict(test.data)

head(test.data)

(predictions)

test.data

data.frame(
    RMSE = caret::RMSE(predictions, test.data$Survival_Ad),
    Rsquare = caret::R2(predictions, test.data$Survival_Ad)
)

ggplot(data.frame(predictions, test.data$LRS.y),
       aes(predictions, test.data.LRS.y)) +
        geom <- point()


############# let's do dyad comparisons now

key <- data.frame(ID=sample_data(PMS)$

