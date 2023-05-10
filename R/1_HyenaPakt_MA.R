## Please uncomment the first time you run this and re-install packages

require(devtools)
## load the devel version
devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")
## devtools::install_github("derele/MultiAmplicon", force= T)

library(ggplot2)
library(dada2)
## library(reshape)
library(phyloseq)
library(data.table)
## library(taxonomizr)
## library(taxize)
library(parallel)
library(pheatmap)
library(tidyr)
library(dplyr)


## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doQualEval <- TRUE

doFilter <- TRUE

doMultiAmpSort <- TRUE

doMultiAmpError <- TRUE ## c("errEst", "direct")

doMultiAmpPipe <- TRUE
    
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
    mutate(fewReads = FilteredReads < quantile(FilteredReads, 0.1))

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

if("errEst"%in%doMultiAmpError){
  ## doing things seperately per run from here to allow sperate error
  ## profiles per run
  errorList <- lapply(unique(MA@sampleData$run), function (run) { 
      i <- which(MA@sampleData$run %in% run)
      errF <-  learnErrors(unlist(getStratifiedFilesF(MA[, i])), nbase=1e8,
                           verbose=0, multithread = 84)
      errR <- learnErrors(unlist(getStratifiedFilesR(MA[, i])), nbase=1e8,
                          verbose=0, multithread = 84)
      list(errF, errR)
  })
    MAList <- lapply(seq_along(unique(MA@sampleData$run)), function (j) {
        run <- unique(MA@sampleData$run)[j]
        i <- which(MA@sampleData$run %in% run)
        dadaMulti(MA[,i], Ferr=errorList[[j]][[1]],
                  Rerr=errorList[[j]][[2]],  pool=FALSE,
                  verbose=0, mc.cores = 84)
    })
  ## combining into one MA object again
  saveRDS(MAList, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MAList_error.Rds")
## } else if ("direct"%in%doMultiAmpError){
##     MA <- derepMulti(MA, mc.cores=24)
##     MA <- dadaMulti(MA, dadaMulti(MA, Ferr=NULL, selfConsist=TRUE,
##                                   Rerr=NULL,  pool=TRUE,
##                                   verbose=0, mc.cores = 64))
##     saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MAdirectErr.Rds")
}else {
    MAList <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MAList_error.Rds")
##    MA <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MAdirectErr.Rds")
}

if(doMultiAmpPipe){
    MAMerged <- Reduce("concatenateMultiAmplicon", MAList)
    MAMerged <- mergeMulti(MAMerged, mc.cores=84)

    propMerged <- calcPropMerged(MAMerged)

    MA.P <- mergeMulti(MAMerged, justConcatenate=propMerged<0.7, mc.cores=84)
    
    ## consider removing the clutter (lot of space in RAM)
    ## rm(MAMerged, MAListMerged, MAList)

    MA.P <- makeSequenceTableMulti(MA.P, mc.cores=84)

    ## get the sequence table fill it, bind it, coerce it to integer
    STF <- getSequenceTable(MA.P, dropEmpty=FALSE)
    STFU <- do.call(cbind, STF)
    mode(STFU) <- "integer"

    ## collapse same identical ASVs transpose data: calculate rowsum
    ## per group (colnames of the original data), then transpose the
    ## result back to original structure.
    STFU <- t(rowsum(t(STFU), group = colnames(STFU), na.rm = TRUE))

    isCruelBimera <- dada2::isBimeraDenovoTable(STFU,
                                                multithread=TRUE, minSampleFraction=0.5,
                                                allowOneOff=TRUE, maxShift = 32,
                                                ignoreNNegatives=4)

    isPooledBimera <- dada2::isBimeraDenovo(STFU, multithread=86, 
                                            allowOneOff=TRUE, maxShift = 32)

    superCruel <- isCruelBimera  | isPooledBimera
    NoBimeras <- colnames(STFU[, !superCruel])

    MA.final <- MA.P
    
    MA.final@sequenceTableNoChime <- sapply(getSequenceTable(MA.final), function (x) {
        noBim <- intersect(NoBimeras, colnames(x))
        x[, noBim, drop=FALSE] ## DROP FALSE TO KEEP MATRIX, ARRAY structure!!
    })
  
    saveRDS(MA.final, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_piped.Rds")
} else {
    MA.final <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_piped.Rds")
}

## ## FIXME in package!
## trackingF <- getPipelineSummaryX(MA.final) 
## PipSum <- plotPipelineSummary(trackingF) + scale_y_log10()
## ggsave("Figures/Pipeline_track.pdf", PipSum, height = 15, width = 15)

if(doTax){
    unlink("/SAN/Victors_playground/Metabarcoding/AA_Hyena/Hyena_in.fasta")
    unlink("/SAN/Victors_playground/Metabarcoding/AA_Hyena/Hyena_out.blt")
MA.A <- blastTaxAnnot(MA.final,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/Hyena_in.fasta",
                    outblast = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/Hyena_out.blt",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 64)

### more sample data
saveRDS(MA.A, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_blast.Rds")
} else {
MA.A <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_blast.Rds")
    }

#### alternative taxonomic annotation
# bit messy here now.
# we first need to know the target of each amplicon
primerL <- read.csv("/SAN/Susanas_den/gitProj/Eimeria_AmpSeq/data/primerInputUnique.csv")
ptable$Primer_name <- paste(ptable$Name_F, ptable$Name_R, sep=".")
primerL$Primer_name[6]<- "18S_0067a_deg_3Mod_53_F.NSR399_3Mod_53_R"
primerL$Primer_name[120]<- "Bgf_132_F.Bgr_132_R"
p.df <- primerL[which(primerL$Primer_name%in%ptable$Primer_name),]
Unk <- ptable[which(!ptable$Primer_name %in% primerL$Primer_name),]
Unk$Gen <- c("COI", "na", "tRNA", "COI")
Unk <- Unk[,c("Primer_name", "Gen")]
P.df <- rbind(p.df[,c("Primer_name", "Gen")], Unk)
# now reorder
idx <- match(ptable$Primer_name, P.df$Primer_name)
P.df <- P.df[idx,]
all(P.df$Primer_name==ptable$Primer_name)# sanity check

# now the taxonomic annotation per target
taxT1 <- list()
seqs <- getSequencesFromTable(MA.final)
seqs <- lapply(seqs, DNAStringSet)
for (i in 1:48){
    if (P.df$Gen[i]=="16S"){
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
               "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa",
                                         multithread=20,
                                         tryRC = TRUE,
                                         verbose=TRUE))
    }
    else if (P.df$Gen[i]=="18S"){
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
               "/SAN/Susanas_den/AmpMarkers/RESCRIPt/SSURef_NR99/Fastas/Slv138.dada2.fa",
                                         multithread=20,
                                         tryRC = TRUE,
                                         verbose=TRUE))
    }
    else if (P.df$Gen[i]=="28S"){
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
               "/SAN/Susanas_den/AmpMarkers/RESCRIPt/LSURef_NR99/Fastas/Slv138LSU.dada2.fa",
                                         multithread=20,
                                         tryRC = TRUE,
                                         verbose=TRUE))
    }
    else if (P.df$Gen[i]=="ITS"){
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
               "/SAN/Susanas_den/AmpMarkers/UNITE/sh_general_release_s_all_10.05.2021/sh_general_release_dynamic_s_all_10.05.2021.fasta",
                                         multithread=20,
                                         tryRC = TRUE,
                                         verbose=TRUE))
    }
    else {
        try(taxT1[[i]] <- assignTaxonomy(seqs[[i]],
               "/SAN/Susanas_den/AmpMarkers/RESCRIPt/other/Fastas/other.dada2.fa",
                                         multithread=90,
                                         tryRC = TRUE,
                                         verbose=TRUE))
    }
}
MA.final@taxonTable <- taxT1

saveRDS(MA.final, file="/SAN/Susanas_den/gitProj/AA_Hyenas_Pakt/tmp/MA_Tax.Rds")



source("/SAN/Susanas_den/gitProj/Eimeria_AmpSeq/R/toPhyloseq.R") # function is broken in the package

PH <- TMPtoPhyloseq(MA.final, samples=colnames(MA.final))

sample_names(PH)

## Does single Amplicon-Data differ from multiAmplicon
pdf("./Figures/Bacterial_ampMethod_Richness.pdf")
plot_richness(subset_taxa(PH, superkingdom%in%"Bacteria"),
              measures=c("Observed", "Chao1", "Shannon"),
              "ampMethod")
dev.off()


## this STILL! buggs FIXME in package!!
PH.list <- TMPtoPhyloseq(MA.final, samples=colnames(MA.final), multi2Single=FALSE)

## Some first quick view at bacterial taxa...
TT <- tax_table(PH)

BacTT <- TT[TT[, "Kingdom"]%in%"k__Bacteria", ]

genusTab <- table(BacTT[, "Genus"])

#write.csv(genusTab[order(genusTab, decreasing=TRUE)],
#          file="genus_table.csv", row.names=FALSE)


## First collapsing technical replcates!!
sample_data(PH)$Sample <- gsub("\\.2", "", sample_data(PH)$Sample)

##  testing the oucome if only using singleAmpData
PSS <- subset_samples(PH, ampMethod%in%"MultiAmp")

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
    if(!is.null(access(PS, "tax_table", errorIfNULL=FALSE))){
        phyloseq(otu_table(OTN, taxa_are_rows=FALSE),
                 sample_data(sdatN),
                 tax_table(PS))
    } else {
        phyloseq(otu_table(OTN, taxa_are_rows=FALSE),
                 sample_data(sdatN))
    }
}

PM <- sumTecRep(PH, by.sample="Sample")
PMS <- sumTecRep(PSS, by.sample="Sample")

## Now adding the annotation realy
## SDat <- read.csv("Data/Covariates_int_biomes.csv")

SDat <- read.csv("Data/microbiome_tagged_fitness_2021-06-30.csv")


SDat <- read.csv("/SAN/Susanas_den/gitProj/AA_Hyenas_Pakt/Data/Covariates_int_biomes_21-05-04.csv")

SDat$sample_ID.x[SDat$sample_ID.x=="C47"]  <- "C0047"
SDat$sample_ID.x[SDat$sample_ID.x=="C52"]  <- "C0052"
SDat$sample_ID.x[SDat$sample_ID.x=="C129"]  <- "C0129"

SDat$CSocialRank <- rowMeans(SDat[, c("social_rank_hyena_ID",
                                      "social_rank_genetic_mum")],
                             na.rm=TRUE)

newSdat <- merge(sample_data(PM), SDat, by.x=0,
                 by.y="sample_ID.x", all.x=TRUE)
rownames(newSdat) <- newSdat$Row.names
newSdat$Row.names <- NULL

sample_data(PM) <- newSdat
sample_data(PMS) <- newSdat

setdiff(SDat$sample_ID.x, sample_data(PM)$Sample)
## Sequence covers them all!

non.immuno <- setdiff(sample_data(PM)$Sample, SDat$sample_ID.x)
grep("Negative", non.immuno, value=TRUE, invert=TRUE)
##   "B3456" "B6423" "X6674"

saveRDS(PMS, "/SAN/Susanas_den/gitProj/AA_Hyenas_Pakt/tmp/PMS.rds")
saveRDS(PM, "/SAN/Susanas_den/gitProj/AA_Hyenas_Pakt/tmp/PM.rds")

### ONE DATASET to start with P --- ONLY MULTIAMPLICON FOR NOW

### exclude super weird taxa
P <- subset_taxa(PMS, superkingdom%in%c("Bacteria", "Eukaryota"))

## Exclude off-target Eukaryote taxa 
P <- subset_taxa(P, !phylum%in%c("Chordata",
                                 "Chlorophyta",
                                 "Streptophyta"
                                 ))

PS <- tax_glom(PMS, "species")
PG <- tax_glom(PMS, "genus")

### SIX datasets to rule them all!!!
PBac <- subset_taxa(P, superkingdom%in%c("Bacteria"))
PEuk <- subset_taxa(P, superkingdom%in%c("Eukaryota"))

PSBac <- subset_taxa(PS, superkingdom%in%c("Bacteria"))
PSEuk <- subset_taxa(PS, superkingdom%in%c("Eukaryota"))

PGBac <- subset_taxa(PG, superkingdom%in%c("Bacteria"))
PGEuk <- subset_taxa(PG, superkingdom%in%c("Eukaryota"))



library(phyloseq)
EukTib <- as_tibble(psmelt(PSEuk))


EukTib %>% filter(genus%in%"Ancylostoma") %>%
    select(Abundance, Sample, OTU, Ancylostoma_egg_load) %>%
    group_by(Sample) %>%
    summarize(sumAbu=sum(Abundance),
              Sample=unique(Sample),
              AncyMic=unique(Ancylostoma_egg_load)
              ) %>% na.omit() -> AncyloCorDat


AncyloCorDat %>% select_if(is.numeric) %>% cor()

nrow(AncyloCorDat)

cor.test(AncyloCorDat$sumAbu, AncyloCorDat$AncyMic)

summary(lm(log10(sumAbu+1)~log10(AncyMic+1), AncyloCorDat))

pdf("Figures/Ancylostoma_Cor.pdf")
ggplot(AncyloCorDat, aes(AncyMic+1, sumAbu+1)) +
    geom_point() +
    stat_smooth(method="lm") +
    scale_y_log10("Ancylostoma sequence abundance") +
    scale_x_log10("Ancylostoma egg load")
dev.off()


EukTib %>% filter(order%in%c("Eucoccidiorida")) %>% 
    select(Abundance, Sample, OTU, Cystoisospora_oocyst_load) %>%
    group_by(Sample) %>%
    summarize(sumAbu=sum(Abundance),
              Sample=unique(Sample),
              CystoMic=unique(Cystoisospora_oocyst_load)
              ) %>% na.omit() -> CoccidiaCorDat

CoccidiaCorDat %>% select_if(is.numeric) %>% cor()

pdf("Figures/Coccidia_Cor.pdf")
CoccidiaCorDat %>%
    ggplot(aes(CystoMic+1, sumAbu+1)) +
    geom_point() +
    stat_smooth(method="lm") +
    scale_y_log10("Coccidia sequence abundance") +
    scale_x_log10("Cystoisospora oocyst load")
dev.off()

cor.test(CoccidiaCorDat$CystoMic, CoccidiaCorDat$sumAbu, method="spearman")

EukTib %>% filter(genus%in%"Cystoisospora") %>% 
    select(Abundance, Sample, OTU, Cystoisospora_oocyst_load) %>%
    group_by(Sample) %>%
    summarize(sumAbu=sum(Abundance),
              Sample=unique(Sample),
              CytoMic=unique(Cystoisospora_oocyst_load)
              ) %>% na.omit() -> CystoisoCorDat

pdf("Figures/Cystoisospora_Cor.pdf")
CystoisoCorDat %>%
    ggplot(aes(CytoMic+1, sumAbu+1)) +
    geom_point() +
    stat_smooth(method="lm") +
    scale_y_log10("Cytoisospora sequence abundance") +
    scale_x_log10("Cystoisospora oocyst load")
dev.off()

cor.test(CystoisoCorDat$CytoMic, CystoisoCorDat$sumAbu, method="spearman")

