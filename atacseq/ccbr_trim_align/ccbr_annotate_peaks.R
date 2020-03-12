
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ChIPseeker"))

suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg38.knownGene"))
suppressPackageStartupMessages(library("TxDb.Mmusculus.UCSC.mm9.knownGene"))
suppressPackageStartupMessages(library("TxDb.Mmusculus.UCSC.mm10.knownGene"))

suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("org.Mm.eg.db"))

parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument("-n", "--narrowpeak", required=TRUE,
                    dest="narrowpeak", help="narrowpeak file")
parser$add_argument("-a", "--annotated", required=TRUE, dest="annotated",
                    help="annotated output file")
parser$add_argument("-u", "--uptss", required=FALSE, type="integer", default=2000, help="upstream bases from TSS")
parser$add_argument("-d", "--downtss", required=FALSE, type="integer", default=2000, help="upstream bases from TSS")
parser$add_argument("-t", "--toppromoterpeaks", required=FALSE, type="integer", default=1000, help="filter top N peaks in promoters for genelist output")
parser$add_argument("-l", "--genelist", required=TRUE, help="list of genes with peaks in promoter regions")
parser$add_argument("-f", "--atypefreq", required=TRUE, help="frequency of different annotation types")
parser$add_argument("-g", "--genome", required=TRUE, dest="genome",
                    help="hg38/hg19/mm10/mm9")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

if (args$genome=="mm9" | args$genome=="mm10"){
  adb="org.Mm.eg.db"
}
if (args$genome=="hg19" | args$genome=="hg38"){
  adb="org.Hs.eg.db"
}
if (args$genome=="hg19") {tdb=TxDb.Hsapiens.UCSC.hg19.knownGene}
if (args$genome=="hg38") {tdb=TxDb.Hsapiens.UCSC.hg38.knownGene}
if (args$genome=="mm9") {tdb=TxDb.Mmusculus.UCSC.mm9.knownGene}
if (args$genome=="mm10") {tdb=TxDb.Mmusculus.UCSC.mm10.knownGene}


np=read.table(args$narrowpeak,sep="\t")
np=np[,seq(1,10)]
colnames(np)=c("chrom",
               "chromStart",
               "chromEnd",
               "name",
               "score",
               "strand",
               "signalValue",
               "pValue",
               "qValue",
               "peak")
np=np[order(-np$qValue),]
np$peakID=paste(np$chrom,":",np$chromStart,"-",np$chromEnd,sep="")

peaks=GRanges(seqnames=np$chrom,ranges=IRanges(np$chromStart,np$chromEnd))

# using annotatePeak from ChIPseeker
pa=annotatePeak(peak = peaks,
                tssRegion = c(-2000,2000),
                TxDb = tdb,
                level = "transcript",
                genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                annoDb = adb, 
                sameStrand = FALSE, 
                ignoreOverlap = FALSE,
                ignoreUpstream = FALSE, 
                ignoreDownstream = FALSE, 
                overlap = "TSS")

padf=as.data.frame(pa)
padf$peakID=paste(padf$seqnames,":",padf$start,"-",padf$end,sep="")
merged=merge(padf,np,by="peakID")
merged=merged[,
              c("peakID",
                "chrom",
                "chromStart",
                "chromEnd",
                "width",
                "annotation",
                "geneChr",
                "geneStart",
                "geneEnd",
                "geneLength",
                "geneStrand",
                "geneId",
                "transcriptId",
                "distanceToTSS",
                "ENSEMBL",
                "SYMBOL",
                "GENENAME",
                "score",
                "signalValue",
                "pValue",
                "qValue",
                "peak"
              )]
# Adding the hash to the first colname
colnames(merged)=c("#peakID",
                   "chrom",
                   "chromStart",
                   "chromEnd",
                   "width",
                   "annotation",
                   "geneChr",
                   "geneStart",
                   "geneEnd",
                   "geneLength",
                   "geneStrand",
                   "geneId",
                   "transcriptId",
                   "distanceToTSS",
                   "ENSEMBL",
                   "SYMBOL",
                   "GENENAME",
                   "score",
                   "signalValue",
                   "pValue",
                   "qValue",
                   "peak")


# merge annotation with narrowPeak file 
merged=merged[order(-merged$qValue),]
write.table(merged,args$annotated,sep = "\t",quote = FALSE, row.names = FALSE)
l=paste("# Median peak width : ",median(merged$width),sep="")
write(l,args$annotated,append=TRUE)
l=paste("# Median pValue : ",median(merged$pValue),sep="")
write(l,args$annotated,append=TRUE)
l=paste("# Median qValue : ",median(merged$qValue),sep="")
write(l,args$annotated,append=TRUE)


# get promoter genes 
# ... all lines with annotation starting with "Promoter"
promoters1=dplyr::filter(merged,grepl("Promoter",annotation))
# ... all lines with annotation is "5' UTR"
promoters2=merged[merged$annotation=="5' UTR",]
promoters=rbind(promoters1,promoters2)
promoters=promoters[order(-promoters$qValue),]
promoters=head(promoters,n=args$toppromoterpeaks)
promoter_genes=unique(promoters[,c("ENSEMBL","SYMBOL")])
colnames(promoter_genes)=c("#ENSEMBL","SYMBOL")
write.table(promoter_genes,args$genelist,sep = "\t",quote = FALSE, row.names = FALSE)
l=paste("# Median peak width : ",median(promoters$width),sep="")
write(l,args$genelist,append=TRUE)
l=paste("# Median pValue : ",median(promoters$pValue),sep="")
write(l,args$genelist,append=TRUE)
l=paste("# Median qValue : ",median(promoters$qValue),sep="")
write(l,args$genelist,append=TRUE)

# annotation type frequence table

l=paste("#annotationType","frequency","medianWidth","medianpValue","medianqValue",sep="\t")
write(l,args$atypefreq)
atypes=c("3' UTR",
         "5' UTR",
         "Distal Intergenic",
         "Downstream (<1kb)",
         "Downstream (1-2kb)",
         "Downstream (2-3kb)",
         "Promoter (<=1kb)",
         "Promoter (1-2kb)")
for (ann in atypes) {
  x=merged[merged$annotation==ann,]
  w=median(x$width)
  p=median(x$pValue)
  q=median(x$qValue)
  l=paste(gsub(" ","",ann),nrow(x),w,p,q,sep="\t")
  write(l,args$atypefreq,append=TRUE)
}
for (ann in c("Exon","Intron")){
  x=dplyr::filter(merged,grepl(ann,annotation))
  w=median(x$width)
  p=median(x$pValue)
  q=median(x$qValue)
  l=paste(gsub(" ","",ann),nrow(x),w,p,q,sep="\t")
  write(l,args$atypefreq,append=TRUE)  
}

