library(GGally)
library(parallel)

if(!exists("PH", mode="S4")){
    source("1_HyenaPakt_MA.R")
}



## ## an analysis of (technical) replication of exact amplicon sequences

SPEZI <- psmelt(PH)

## OTUs were made unique but can acutally be repeated between 
as_tibble(SPEZI) %>%
    transform(OTU=sub("\\.\\d", "", OTU)) ->
    SPEZI

## ## I dont get this:
## ## see whether OTUs are now duplicated in each sample
## SPEZI %>%
##     group_by(sample_Sample) %>%
##     summarize(dupOTU=any(duplicated(OTU[Abundance>0]))) %>% print(n=500)

## Negative control samples are a concern in our dataset... this
## figure demonstates that in the Single-Amplicon sequencing runs they
## were indistinguishable from proper samples
pdf("./Figures/NegativeWeirdness_run.pdf")
SPEZI %>%
    transform(NC= grepl("Negative", sample_Sample))  %>%
    group_by(sampleID, run)  %>%
    summarize(totalCount=sum(Abundance),
              nOTUs=uniqueN(OTU[Abundance>0]),
              nGenus=uniqueN(genus[Abundance>0]),
              nPhylum=uniqueN(phylum[Abundance>0]),
              NC=unique(NC),
              ampMethod=unique(ampMethod),
              run=unique(run)) %>%
    ggplot(aes(NC, totalCount, color=ampMethod)) +
    geom_jitter() +
    facet_wrap(~run) +
    scale_x_discrete("Is it a negative control?") +
    scale_y_continuous("Number of reads")
dev.off()




#### repeatability for repeated samples: 1. ASVs ##########
getSampleCor <- function (tibble, sample, cols2group = c("ampMethod","OTU")){
    filter(tibble, sample_Sample%in%sample) %>%
        group_by(across(all_of(cols2group))) %>%
        summarize(abu=sum(Abundance, na.rm=TRUE), .groups="drop") %>%
        spread(key=cols2group[1], value=abu) %>%
        transform(sample=sample) %>% as_tibble()
}

getSampleCor(SPEZI, "B841") %>%
    select_if(is.numeric) %>%
    cor()

getSampleCor(SPEZI, "B841", c("run", "OTU")) %>%
    select_if(is.numeric) %>%
    cor()


    
pdf("./Figures/PairsOneSample.pdf")
getSampleCor(SPEZI, "B841", c("run", "OTU")) %>%
    select_if(is.numeric) %>%
ggpairs() +
    scale_y_log10() +
    scale_x_log10()
dev.off()


SAHY <- unique(SPEZI$sample_Sample)

corDataAmpASV <- mclapply(SAHY, function (x){
    getSampleCor(SPEZI, x)
}, mc.cores=48)

names(corDataAmpASV) <- SAHY

### does not work as left out when no ASVs found for one ampMethod and
### sample!! (how this could happen is another weird story
## Reduce(rbind, corDataAmpASV)
