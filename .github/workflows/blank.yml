library(VariantAnnotation) ## load package for reading vcf file
library(TxDb.Hsapiens.UCSC.hg19.knownGene) ## load package for gene model
library(org.Hs.eg.db) ## load package for gene name
library(rtracklayer)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
## set the track range
gr <- GRanges("22", IRanges(50968014, 50970514, names="TYMP"))
## read in vcf file
tab <- TabixFile(fl)
vcf <- readVcf(fl, "hg19", param=gr)
## get GRanges from VCF object 
mutation.frequency <- rowRanges(vcf)
## keep the metadata
mcols(mutation.frequency) <-
    cbind(mcols(mutation.frequency),
          VariantAnnotation::info(vcf))
## set colors
mutation.frequency$border <- "gray30"
mutation.frequency$color <-
    ifelse(grepl("^rs", names(mutation.frequency)),
           "lightcyan", "lavender")
## plot Global Allele Frequency based on AC/AN
mutation.frequency$score <- round(mutation.frequency$AF*100)
## change the SNPs label rotation angle
mutation.frequency$label.parameter.rot <- 45
## keep sequence level style same
seqlevelsStyle(gr) <- seqlevelsStyle(mutation.frequency) <- "UCSC"
seqlevels(mutation.frequency) <- "chr22"
## extract transcripts in the range
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                         org.Hs.eg.db, gr=gr)
## subset the features to show the interested transcripts only
features <- c(range(trs[[1]]$dat), range(trs[[5]]$dat))
## define the feature labels
names(features) <- c(trs[[1]]$name, trs[[5]]$name)
## define the feature colors
features$fill <- c("lightblue", "mistyrose")
## define the feature heights
features$height <- c(.02, .04)
## set the legends
legends <- list(labels=c("known", "unkown"),
                fill=c("lightcyan", "lavender"),
                color=c("gray80", "gray80"))
## lollipop plot
lolliplot(mutation.frequency,
          features, ranges=gr, type="circle",
          legend=legends)
