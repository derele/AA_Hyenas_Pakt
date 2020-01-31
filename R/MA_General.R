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
doFilter <- TRUE

doMultiAmp <- TRUE

doTax <- TRUE
## But remember: if you change the MultiAmplicon Analysis, the
## taxonomic annotation might be out of sync...

###################Full run Microbiome#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

## change according to where you downloaded
path <- "/SAN/Victors_playground/Metabarcoding/AA_Hyena/2018_22_Hyena_1/"

fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 


samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("S\\d+_", "\\1", basename(samples))

#Extra step in the pipeline: quality plots of the reads 
## plotQualityProfile(fastqF[[1]]) ### Really low quality after 150bp :/ for test run data 
## plotQualityProfile(fastqF[[200]])
## plotQualityProfile(fastqR[[1]])
## plotQualityProfile(fastqR[[200]])

#Creation of a folder for filtrated reads 
filt_path <- "/SAN/Victors_playground/Metabarcoding/AA_Hyena/filtered_Hyena_1"
#filt_path <- "/SAN/Victors_playground/Metabarcoding/AA_Hyena/filtered_Hyena_2"

#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering
if(doFilter){
  filter.track <- lapply(seq_along(fastqF),  function (i) {
    filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(200,200), minLen=c(200,200), 
                  maxN=0, maxEE=2, truncQ=2, 
                  compress=TRUE, verbose=TRUE)
  })
  saveRDS(filter.track, file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/filter.Rds")
} else {
  filter.track <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Hyena/filter.Rds")
}

##Check the proportion of reads that passed the filtering 
filter <- do.call(rbind, filter.track)
colSums(filter)[2]/colSums(filter)[1]

###Uuuuu :S just 40% passed for the test run. 

names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel

#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/primer_list.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)


##Multi amplicon pipeline
if(doMultiAmp){
  MA <- MultiAmplicon(primer, files)
  #filedir <- "/SAN/Victors_playground/Metabarcoding/AA_Hyena/stratified_Hyena_2"
  filedir <- "/SAN/Victors_playground/Metabarcoding/AA_Hyena/stratified_Hyena_1"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
  
  errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                       verbose=0, multithread = 12)
  errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                      verbose=0, multithread = 12)
  
  MA <- derepMulti(MA, mc.cores=12) 
  
  MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                  verbose=0, mc.cores=12)
  
  MA <- mergeMulti(MA, mc.cores=12) 
  
  propMerged <- MultiAmplicon::calcPropMerged(MA)
  
  MA <- mergeMulti(MA, justConcatenate=propMerged<0.8, mc.cores=12) 
  
  MA <- makeSequenceTableMulti(MA, mc.cores=12) 
  
  MA <- removeChimeraMulti(MA, mc.cores=12)
  
  saveRDS(MA, "/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_1.RDS") ##Pool Hyena 1 preliminary run
  #saveRDS(MA, "/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_2.RDS") ##Pool Hyena 2 full run (2nd batch of data)
} else{
  MA <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_1.RDS") ###START from here now!
  #MA <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Hyena/MA_2.RDS")
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
                    db = "/SAN/db/blastdb/",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/Hyena_1_in.fasta",
                    outblast = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/blast2_1_out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 20)

###ERROR in taxassignment
#BLAST Database error: No alias or index file found for nucleotide database [BLASTDB] in search path [/localstorage/victor:/SAN/db/blastdb:]

#FINISHED running blast
#Error in read.table(file = file, header = header, sep = sep, quote = quote,  : 
#                     no lines available in input

#saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA1_1Tax.Rds") ##Just Test run 

##Start from here after the taxonomic annotation
#MA<- readRDS(file= "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA1_1Tax.Rds") ###Test run


###Load sample information

sample.data <- read.csv("~/AA_Hyena_Pakt/Index_Pool_1.csv",
                        dec=",", stringsAsFactors=FALSE)
##Adding sample data 

MAsample <- addSampleData(MA, sample.data)
