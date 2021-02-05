## Please uncomment the first time you run this and re-install packages

# require(devtools)
## devtools::install_github("derele/MultiAmplicon", force= T)

library(MultiAmplicon)

library(ggplot2)
library(dada2)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)
library(pheatmap)
library(tidyr)
library(dplyr)


## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doQualEval <- FALSE

doFilter <- FALSE

doMultiAmpSort <- FALSE

doMultiAmpError <- FALSE

doMultiAmpPipe <- FALSE
    
doTax <- TRUE

###################Full run Microbiome#######################
## Preparation of files. These are the same steps that are followed by
## the DADA2 pipeline change according to where you downloaded


path <- c(
    ## Hyena Pool 2 (Single amplicon run) 2nd Full sequencing Run (Good run)
    "2018_22_hyena", 
    ## Hyena Pool 1 (Multiamplicon run) Preliminary test sequencing Run
    "2018_22_Hyena",
    ## Hyena Pool 1 (Multiamplicon run) Extra sequencing Run for more reads 1
    "2018_22_hyena1_extra3_part1",
    ## Hyena Pool 1 (Multiamplicon run) Extra sequencing Run for more reads 2
    "2018_22_hyena1_extra3_part2",
    ## Hyena Pool 1 (Multiamplicon run) Extra sequencing Run for more reads 3
    "2018_22_hyena1_extra3_part3",
    ## Hyena Pool 2 (Single amplicon run) Preliminary test sequencing Run
    "2018_22_hyena2_main_run_1",
    ## Hyena Pool 2 (Single amplicon run) 1st Full sequencing Run (Susan report some problems and few reads)
    "2018_22_hyena2_run2",
    ## Hyena Pool 1 (Multiamplicon run) Full sequencing Run
    "2018_22_hyena_main_run",
    ## Hyena Pool 2 (Single amplicon run) Extra sequencing Run for more reads 
    "2018_22_P2_extra") 


fullpath <- paste0("/SAN/Victors_playground/Metabarcoding/AA_Hyena/", path)

names(fullpath) <- path

fastqList <- lapply(fullpath, function (path) { 
    fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) 
    fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE)
    fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE)
    list(fastqF=fastqF, fastqR=fastqR)
})

if(doQualEval){
    readlenght <- lapply(fastqList, function (x) {
        con <- file(x[["fastqF"]][1],"r")
        ## first line
        secondLine <- readLines(con, n=2)[[2]]
### simple check whether it's V2 or V3 data
        nchar(secondLine)
    })

    allFastqF <- lapply(fastqList, function (x) {
        readFastq(x[["fastqF"]])
    })

    allFastqR <- lapply(fastqList, function (x) {
        readFastq(x[["fastqR"]])
    })


    sampleQual <- function (x) {
        ## sample quality scores of 100,000 sequences 
        qmat <- as(quality(x)[sample(100000)], "matrix")
        cols <- seq(1, ncol(qmat), by=10)
        sapply(cols, function (i) {
            mean(qmat[, i], na.rm=TRUE)
        })
    }

    qualityF <- lapply(allFastqF, sampleQual)
    qualityR <- lapply(allFastqR, sampleQual)

    shouldL <- max(unlist(lapply(qualityF, length)))

    qualityFilledF <- lapply(qualityF, function (x) {
        c(x, rep(NA, times=shouldL - length(x)))
    })

    qualityFilledR <- lapply(qualityR, function (x) {
        c(x, rep(NA, times=(shouldL - length(x))))
    })


    qualityDFF <- Reduce("cbind",  qualityFilledF)
    qualityDFR <- Reduce("cbind",  qualityFilledR)

    colnames(qualityDFF) <- path
    colnames(qualityDFR) <- path

    qualityDFFL <- reshape2::melt(qualityDFF)
    qualityDFFL$direction <- "forward"

    qualityDFRL <- reshape2::melt(qualityDFR)
    qualityDFRL$direction <- "reverse"

    qualityDFL <- rbind(qualityDFFL, qualityDFRL)

    qualityDFL$position <- qualityDFL$Var1*10 -10

    ggplot(qualityDFL, aes(position, value, color=Var2)) +
        geom_line() +
        facet_wrap(~direction)
}

## concluding from this that we can truncate at 220 and 200 for
## reverse and forward respectively

## concluding that we have to remove runs
exclude_runs <-  c("2018_22_Hyena",
                   "2018_22_hyena_main_run",
                   "2018_22_hyena2_main_run_1")

fastqList <- fastqList[!names(fastqList)%in%exclude_runs]

samplesList <- lapply (fastqList, function (x){
    samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(x[["fastqF"]]))
    paste(basename(dirname(x[["fastqF"]])), samples, sep="_-")
})

fastqFall <- unlist(lapply(fastqList, "[[", "fastqF"))
fastqRall <- unlist(lapply(fastqList, "[[", "fastqR"))

samplesAll <- unlist(samplesList)

#Creation of a folder for filtrated reads 
filt_path <- "/SAN/Victors_playground/Metabarcoding/AA_Hyena/filtered_Hyena_all"

if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(samplesAll, "_F_filt.fastq.gz"))
names(filtFs) <- samplesAll
filtRs <- file.path(filt_path, paste0(samplesAll, "_R_filt.fastq.gz"))
names(filtRs) <- samplesAll

if(doFilter){
  filter.track <- lapply(seq_along(fastqFall),  function (i) {
      filterAndTrim(fastqFall[i], filtFs[i], fastqRall[i], filtRs[i],
                    truncLen=c(220,200), minLen=c(220,200), 
                    maxN=0, maxEE=2, truncQ=2, 
                    compress=TRUE, verbose=TRUE,
                    matchIDs=TRUE) ## forward and reverse not matching otherwise 
  })
  saveRDS(filter.track, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/filter.Rds")
} else {
  filter.track <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/filter.Rds")
}

##Check the proportion of reads that passed the filtering 
filter <- as.data.frame(do.call(rbind, filter.track))
sum(filter[,"reads.out"])/sum(filter[,"reads.in"])

### Over 80% passed for all runs...
filter$run <- unlist(lapply(strsplit(samplesAll, "_-"), "[", 1))
##filter$run <- gsub("\\d_part\\d", "", filter$run)

## but less in some runs...
by(filter, filter$run, function (x) sum(x[,"reads.out"]/sum(x[,"reads.in"])))


files <- PairedReadFileSet(filtFs, filtRs)

### SAMPLES
sampleIDs <- read.csv("Data/Index_Pool_1.csv")
sampleIDs <- merge(sampleIDs, read.csv("Data/Index_Pool_2.csv"), by="Sample", all=TRUE)

sampleIDs <- pivot_longer(sampleIDs, cols=c(BeGenDiv_Pool_1, BeGenDiv_Pool_2)) 

sampleIDs <- as.data.frame(sampleIDs)

colnames(sampleIDs)[colnames(sampleIDs)%in%"value"] <- "sampleID"

filter$sampleID <- gsub(".*?(P\\d)\\.(FLD\\d{4}).*",
                        "\\1_\\2", rownames(filter))

filter$SnumIDs <- gsub("(S\\d{3,4})\\.(P\\d)\\.(FLD\\d{4}).*", "\\1_\\2_\\3",
                       rownames(filter))

sampleIDs <- merge(sampleIDs, filter, by="sampleID", all=TRUE)

sampleIDs$ampMethod <- ifelse(grepl("P1", sampleIDs$sampleID),
       "MultiAmp", ifelse(grepl("P2", sampleIDs$sampleID), "SingleAmp", NA))

#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel

#Primers used in the arrays, primer pairs in single processin are part of this
ptable <- read.csv(file = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/primer_list.csv",
                   sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)

M1 <- MultiAmplicon(primer, files)

rownames(sampleIDs) <- make.unique(paste(sampleIDs$run,
                                         gsub("_", "-", sampleIDs$SnumIDs),
                                         sep="_-"))

MA <- addSampleData(M1, sampleIDs)

sumSample <- tibble(sampleIDs) %>%
    group_by(Sample, ampMethod) %>% drop_na() %>%
    summarise(FilteredReads = sum(reads.out), DNA_conc=unique(DNA_conc),
              P260_280 = unique(P260_280), P260_230 = unique(P260_230)) %>%
    transform(fewReads = FilteredReads < quantile(FilteredReads, 0.1))

fewData <- subset(sumSample, fewReads)
goodData <- subset(sumSample, !fewReads)

write.csv(sumSample, file="sequencingOutputBySample.csv")

## devtools::install_github("slowkow/ggrepel")
library(ggrepel)


pos <- position_jitter(width = 0.3, seed = 1)
pdf("Figures/sequencingOutputBySample.pdf", width=7, height=7)
ggplot(goodData, aes(ampMethod, FilteredReads, label=Sample)) +
    geom_point(position=pos, aes(color=fewReads)) +
    geom_label_repel(data = fewData,
                     position = pos) +
    scale_y_log10() +
    labs(color = "#Reads < \n0.1 quantile\n(6072)") +
    theme_bw()
dev.off()


##Multi amplicon pipeline
if(doMultiAmpSort){
  filedir <- "/SAN/Victors_playground/Metabarcoding/AA_Hyena/stratified_All"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  ## This step sort the reads into amplicons based on the number of primer pairs
  MA <- sortAmplicons(MA, n=1e+07, filedir=filedir) 
  pdf("Figures/overview_all_heat.pdf", width=16, height=61)
  pheatmap(log10(getRawCounts(MA)+1)) #, 
  ##          annotation_col=MA@sampleData[, c("run", "reads.in")])
  dev.off()
  saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_sorted.Rds")
} else {
    MA <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_sorted.Rds")
}

## seperate sample data for each run

if(doMultiAmpError){
  ## doing things seperately per run from here to allow sperate error
  ## profiles per run
  errorList <- lapply(unique(MA@sampleData$run), function (run) { 
      i <- which(MA@sampleData$run %in% run)
      errF <-  learnErrors(unlist(getStratifiedFilesF(MA[, i])), nbase=1e8,
                           verbose=0, multithread = 12)
      errR <- learnErrors(unlist(getStratifiedFilesR(MA[, i])), nbase=1e8,
                          verbose=0, multithread = 12)
      list(errF, errR)
  })

    MAList <- lapply(seq_along(unique(MA@sampleData$run)), function (j) {
        run <- unique(MA@sampleData$run)[j]
        i <- which(MA@sampleData$run %in% run)
        dadaMulti(MA[,i], Ferr=errorList[[j]][[1]],
                  Rerr=errorList[[j]][[2]],  pool=FALSE,
                  verbose=0, mc.cores = 1)
    })

  ## combining into one MA object again
  saveRDS(MAList, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MAList_error.Rds")
} else {
    MAList <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MAList_error.Rds")
}

if(doMultiAmpPipe){
    MAListMerged <- lapply(MAList, mergeMulti, mc.cores=12)
    ###  so there is a strange error if I concatenate the list before
    ###  the merging: the dada and derep objects get out of sync. It
    ###  migh be worth to revisit this.
    MAMerged <- Reduce("concatenateMultiAmplicon", MAListMerged)
    
    propMerged <- MultiAmplicon::calcPropMerged(MAMerged)

    MAListMerged <- lapply(MAList, mergeMulti, justConcatenate=propMerged<0.7)

    MA <- Reduce("concatenateMultiAmplicon", MAListMerged)
    
    ## consider removing the clutter (lot of space in RAM)
    ## rm(MAMerged, MAListMerged, MAList)

        
    MA <- makeSequenceTableMulti(MA, mc.cores=12)
    MA <- removeChimeraMulti(MA, mc.cores=12) 
  
    saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_piped.Rds")
} else {
    MA <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_piped.Rds")
}

trackingF <- getPipelineSummary(MA) 
PipSum <- plotPipelineSummary(trackingF) + scale_y_log10()
ggsave("Figures/Pipeline_track.pdf", PipSum,height = 15, width = 15)


## ## an analysis of (technica) replication WORK ON ME!!!
## foo <- getSequenceTable(MA)[[8]]

## bar <- MA@sampleData

## foobar <- merge(foo, bar, by=0)

## ## the others are the 
## non.seq.cols <- colnames(foobar)[nchar(colnames(foobar)) < 50]
## seq.cols <- colnames(foobar)[nchar(colnames(foobar)) >= 50]

## table(foobar$ampMethod)
## foobarBAZ <- pivot_longer(foobar, cols=all_of(seq.cols), names_repair="unique")
## colnames(foobarBAZ)[c(17, 25)] <- c("pool", "ASV")


## drop_na(foobarBAZ) %>% select(ampMethod, ASV, value, Sample) %>%
##     ## some are in two runs --- so this doesn't fully address things
##     ## happening in sequencing (differently between runs)
##     pivot_wider(names_from=c(ampMethod), values_from=value, values_fn = sum) ->
##     BAT 

## BAT
## ## a value for each sample for each ASV! That's what we need to look at!
## rowwise(BAT) %>% mutate(bothSum=sum(MultiAmp, SingleAmp),
##                         bothPos=MultiAmp>0&&SingleAmp>0) -> BAT

## cor(BAT$MultiAmp, BAT$SingleAmp, use="pairwise.complete.obs")

## cor(filter(BAT, bothPos)$MultiAmp, filter(BAT, bothPos)$SingleAmp,
##     use="pairwise.complete.obs")

## ## well, damn!

## corRaw <- ggplot(filter(BAT, bothPos), aes(x=MultiAmp, y=SingleAmp)) +
##     geom_density2d_filled()   +
##     scale_y_log10() +
##     scale_x_log10()

## ## wow this is sooo BAD. Can't believe it
## ggsave("Figures/CorRAW.pdf", corRaw, height = 15, width = 15)


###New taxonomic assignment
## should be set sysetem-wide, here just for clarity
## Sys.setenv("BLASTDB" = "/SAN/db/blastdb/") 



MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/Hyena_in.fasta",
                    outblast = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/Hyena_out.blt",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 20)

PH <- toPhyloseq(MA, samples=colnames(MA))

## this STILL! buggs
## PH.list <- toPhyloseq(MA, samples=colnames(MA), multi2Single=FALSE)


## Some first quick view at bacterial taxa...
TT <- tax_table(PH)

BacTT <- TT[TT[, "superkingdom"]%in%"Bacteria", ]

BacTT[BacTT[, "genus"]%in%"Vibrio", ]


genusTab <- table(BacTT[, "genus"])

write.csv(genusTab[order(genusTab, decreasing=TRUE)],
          file="genus_table.csv", row.names=FALSE)


## First collapsing technical replcates!!

## This should acutally work with phyloseq's merge_samples function
## but doesn't as this messes up sample_data.
## CANDIDATE FOR INCLUSION IN PACKAGE...!
sumTecRep <- function (PS, by.sample, fun=sum){
    otab <- setDT(apply(otu_table(PS), 2, as.list))
    ## the columns giving numbers for sequences
    numcols <- colnames(otab)[nchar(colnames(otab))>10]
    sdat <- sample_data(PS, errorIfNULL = FALSE)
    otab[, (numcols):=lapply(.SD, as.numeric), .SDcols=numcols]
    otab[, sfac := as.factor(sdat[[by.sample]])]
    setkey(otab, sfac)
    otabN <- otab[, lapply(.SD, fun), by=sfac]
    setkey(otabN, sfac)
    OTN <- as.matrix(otabN, rownames=TRUE)
    ## now select the entries from colums that have the same values in
    ## the sample table...
    sdatN <- by(sdat, sdat[[by.sample]], function(x){
        sapply(x, function (y){
            uy <- unique(y)
            if(length(uy)==1) uy else paste(uy, collapse=";")
        })
    })
    sdatN <- as.data.frame(do.call(rbind, sdatN))
    phyloseq(otu_table(OTN, taxa_are_rows=FALSE),
             sample_data(sdatN),
             tax_table(PS))
}


PM <- sumTecRep(PH, by.sample="Sample")

## Now adding the annotation realy
SDat <- read.csv("Data/Covariates_int_biomes.csv")

newSdat <- merge(sample_data(PM), SDat, by.x=0, by.y="sample_ID", all.x=TRUE)
rownames(newSdat) <- newSdat$Row.names
newSdat$Row.names <- NULL
sample_data(PM) <- newSdat


PG <- tax_glom(PM, "genus")
PG <- subset_taxa(PG, superkingdom%in%c("Bacteria", "Eukaryota"))

PBak <- subset_taxa(PG, superkingdom%in%c("Bacteria"))

PEuk <- subset_taxa(PG, superkingdom%in%c("Eukaryota"))

pdf("Figures/Genera_per_phylum_Bak.pdf")
ggplot(data.frame(tax_table(PBak)), 
       aes(x=phylum, y=..count..)) +
    geom_bar() +
    coord_polar() +
    scale_y_log10("Number of genera per phylum") +
    theme_bw()
dev.off()

pdf("Figures/Genera_per_phylum_Euk.pdf")
ggplot(data.frame(tax_table(PEuk)), 
       aes(x=phylum, y=..count..)) +
    geom_bar() +
    coord_polar() +
    scale_y_log10("Number of genera per phylum") +
    theme_bw()
dev.off()


pdf("Figures/Reads_per_phylum_Bak.pdf")
ggplot(data.frame(psmelt(PBak)), 
       aes(x=phylum, y=..count.., weight=Abundance)) +
    geom_bar() +
    coord_polar() +
    scale_y_log10("Seqeuncing read counts per phylum") +
    theme_bw()
dev.off()

pdf("Figures/Reads_per_phylum_Euk.pdf")
ggplot(data.frame(psmelt(PEuk)), 
       aes(x=phylum, y=..count.., weight=Abundance)) +
    geom_bar() +
    coord_polar() +
    scale_y_log10("Seqeuncing read counts per phylum") +
    theme_bw()
dev.off()


pdf("Figures/Bacterial_richness_categorical.pdf")
plot_richness(PBak, measures=c("Observed", "Chao1", "Shannon"), "age_sampling_cat") +
    geom_boxplot()
dev.off()

pdf("Figures/Bacterial_richness_continuous.pdf")
plot_richness(PBak, measures=c("Observed", "Chao1", "Shannon"), "age_sampling")
dev.off()

pdf("Figures/Euk_richness_RankMom.pdf")
plot_richness(prune_samples(!is.na(sample_data(PEuk)$social_rank_genetic_mum),
                            PEuk), 
                            measures=c("Observed", "Chao1", "Shannon"),
              "social_rank_genetic_mum") + geom_smooth()
dev.off()


pdf("Figures/Euk_richness_RankID.pdf")
plot_richness(prune_samples(!is.na(sample_data(PEuk)$social_rank_hyena_ID),
                            PEuk), 
                            measures=c("Observed", "Chao1", "Shannon"),
              "social_rank_hyena_ID") + geom_smooth()
dev.off()


prune_prune <- function (ps) {
    s <- prune_samples(sample_sums(ps) > 100, ps)
    S <- prune_samples(!grepl("Negative", sample_names(s)), s)
    prune_taxa(taxa_sums(S) > 1, S)
    
}

logBAC <- transform_sample_counts(
    prune_prune(PBak),
    function(x) log10(1+x))

out.bc.log <- ordinate(logBAC, method = "NMDS", distance = "bray",
                       maxit=100)

pdf("Figures/age_bac_ordi.pdf")
plot_ordination(logBAC, out.bc.log, color="age_sampling", shape="age_sampling_cat") +
    labs(col = "Binned Age") +
    scale_shape_discrete(solid=FALSE) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    ggtitle("Ordination on log 1+x transformed abundance of bacterial genera")
dev.off()


logEuk <- transform_sample_counts(
    prune_prune(PEuk),
    function(x) log10(1+x))

out.eu.log <- ordinate(logEuk, method = "NMDS", distance = "bray")

pdf("Figures/age_Euk_ordi.pdf")
plot_ordination(logEuk, out.eu.log, color="age_sampling", shape="age_sampling_cat") +
    labs(col = "Binned Age") +
    scale_shape_discrete(solid=FALSE) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    ggtitle("Ordination on log 1+x transformed abundance of Eukaryote genera")
dev.off()



pdf("Figures/RankMom_Euk_ordi.pdf")
plot_ordination(logEuk, out.eu.log, color="social_rank_genetic_mum")+
    labs(col = "Binned Age") +
    scale_shape_discrete(solid=FALSE) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    ggtitle("Ordination on log 1+x transformed abundance of Eukaryote genera")
dev.off()
