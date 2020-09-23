## Please uncomment the first time you run this and re-install packages

#require(devtools)
#devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)

library(ggplot2)
library(MultiAmplicon)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- FALSE

doMultiAmpSort <- FALSE

doMultiAmpError <- FALSE

doMultiAmpPipe <- TRUE

doTax <- FALSE
## But remember: if you change the MultiAmplicon Analysis, the
## taxonomic annotation might be out of sync...

###################Full run Microbiome#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

## change according to where you downloaded


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
    fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
    fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
    fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads
    list(fastqF=fastqF, fastqR=fastqR)
})

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
                    compress=TRUE, verbose=TRUE)
  })
  saveRDS(filter.track, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/filter.Rds")
} else {
  filter.track <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/filter.Rds")
}

## Error in filterAndTrim(fastqFall[i], filtFs[i], fastqRall[i], filtRs[i],  : 
##   All output files must be distinct.
## In addition: There were 50 or more warnings (use warnings() to see the first 50)


##Check the proportion of reads that passed the filtering 
filter <- as.data.frame(do.call(rbind, filter.track))
sum(filter[,"reads.out"])/sum(filter[,"reads.in"])

### Over 80% passed for all runs...
filter$run <- unlist(lapply(strsplit(samplesAll, "_-"), "[", 1))

## but less in some runs...
by(filter, filter$run, function (x) sum(x[,"reads.out"]/sum(x[,"reads.in"])))


files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel

#Primers used in the arrays, primer pairs in single processin are part of this
ptable <- read.csv(file = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/primer_list.csv",
                   sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)


##Multi amplicon pipeline
if(doMultiAmpSort){
  MA <- MultiAmplicon(primer, files)
  filedir <- "/SAN/Victors_playground/Metabarcoding/AA_Hyena/stratified_All"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  ## This step sort the reads into amplicons based on the number of primer pairs
  MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) 
  saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_sorted.Rds")
} else {
    MA <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_sorted.Rds")
}

  ## seperate sample data for each run
  rownames(filter) <- rownames(MA@sampleData)
  MA <- addSampleData(MA, filter)

  library(pheatmap)
  pdf("Figures/overview_all_heat.pdf", width=16, height=61)
  pheatmap(log10(getRawCounts(MA)+1)) #, 
  ##          annotation_col=MA@sampleData[, c("run", "reads.in")])
  dev.off()

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

  
  MAList <- lapply(unique(MA@sampleData$run), function (run) { 
      i <- which(MA@sampleData$run %in% run)
      derepMulti(MA[, i], mc.cores=12)
  })
      
  MAList <- lapply(seq_along(errorList), function (i) { 
     dadaMulti(MAList[[i]], Ferr=errorList[[i]][[1]],
                      Rerr=errorList[[i]][[2]],  pool=FALSE,
                      verbose=0, mc.cores = 12)
  })

  ## combining into one MA object again
  MA <- Reduce(concatenateMultiAmplicon, MAList)
  saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_error.Rds")
} else {
    MA <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_error.Rds")
}


if(doMultiAmpPipe){
    MA <- mergeMulti(MA, mc.cores=12) 
  
    propMerged <- MultiAmplicon::calcPropMerged(MA)
  
    MA <- mergeMulti(MA, justConcatenate=propMerged<0.8, mc.cores=12) 

    MA <- makeSequenceTableMulti(MA, mc.cores=12) 
  
    MA <- removeChimeraMulti(MA, mc.cores=12)
    saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_piped.Rds")
} else {
    MA <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_piped.Rds")
}


trackingF <- getPipelineSummary(MA) 
plotPipelineSummary(trackingF) 
PipSum <- plotPipelineSummary(trackingF) + scale_y_log10()
#ggsave("Sequencing_summary_Hyena_1.pdf", PipSum, path = "~/AA_HMHZ/", height = 15, width = 15)


Heatmap <- plotAmpliconNumbers(MA)
#ggsave("Sequencing_reads_Hyena_2.pdf", Heatmap, path = "~/AA_Hyena/", height = 15, width = 15)

#pdf("~/AA_Hyena/Sequencing_reads.pdf", plotAmpliconNumbers(MA), height = 15, width = 15)
#dev.off()

###New taxonomic assignment

Sys.setenv("BLASTDB" = "/SAN/db/blastdb/") #To make the annotation work, boss will fix this in the package
#library("vctrs", lib.loc="/usr/local/lib/R/site-library")
#MA <- blastTaxAnnot(MA,  dataBaseDir = Sys.getenv("BLASTDB"), negative_gilist = "/SAN/db/blastdb/uncultured.gi", num_threads = 20)


MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Metabarcoding/AA_Hyenas_EpiRank/Hyena_1_in.fasta",
                    outblast = "/SAN/Metabarcoding/AA_Hyenas_EpiRank/blast2_1_out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 20)


### TODO: fix blastTaxAnnot... somehow had to blast manually and read then

#saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA1_1Tax.Rds") ##Just Test run 

##Start from here after the taxonomic annotation
#MA<- readRDS(file= "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA1_1Tax.Rds") ###Test run


###Load sample information

sample.data <- read.csv("~/AA_Hyena_Pakt/Index_Pool_1.csv",
                        dec=",", stringsAsFactors=FALSE)
##Adding sample data 

MAsample <- addSampleData(MA, sample.data)
