#optimized for R/3.6.1
#for details: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html

suppressPackageStartupMessages(library("S4Vectors"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("GenomeInfoDb"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("DelayedArray"))
suppressPackageStartupMessages(library("matrixStats"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SingleR"))
suppressPackageStartupMessages(library("Seurat"))

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument("-i", "--inputrds", required=TRUE,
                    dest="inputrds", help="inputRDS file containing Seurat object")
parser$add_argument("-s", "--species", required=TRUE, dest="species",
                    help="human/mouse")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

inputRds = as.character(args$inputrds)
species = as.character(args$species)

so=readRDS(inputRds)

runSingleR = function(obj,refFile,fineORmain){
  #convert Seurat object to SingleCellExperiment object, as required by SingleR
  #use SCTransform assay if available. Otherwise use base RNA counts

  if(length(grep("SCT",names(obj@assays))!=0)){
    sce = as.SingleCellExperiment(obj,assay = "SCT")
  }else{
    sce = as.SingleCellExperiment(obj,assay = "RNA")
  }

  ref = refFile
  s = SingleR(test = sce, ref = ref,labels = ref[[fineORmain]]) #SingleR call
  # return(s$pruned.labels)
  return(s$labels)
  # return(s)
}

if(species == "human"){
  #Human Primary Cell Atlas
  ref.se <- readRDS("/opt2/db/HumanPrimaryCellAtlasData.rds")
  so$HPCA_main <- runSingleR(so,ref.se,"label.main")
  so$HPCA <-  runSingleR(so,ref.se,"label.fine")

  #Blueprint/ENCODE database
  ref.se <- readRDS("/opt2/db/BlueprintEncodeData.rds")
  so$BP_encode_main <-  runSingleR(so,ref.se,"label.main")
  so$BP_encode <-  runSingleR(so,ref.se,"label.fine")

  #Monaco Immune Cell data
  ref.se <- readRDS("/opt2/db/MonacoImmuneData.rds")
  so$monaco_main <-  runSingleR(so,ref.se,"label.main")
  so$monaco <-	 runSingleR(so,ref.se,"label.fine")

  #Database for Immune Cell Expression
  ref.se <- readRDS("/opt2/db/DatabaseImmuneCellExpressionData.rds")
  so$dice_main <-  runSingleR(so,ref.se,"label.main")
  so$dice <- runSingleR(so,ref.se,"label.fine")

  #Novershtern Hematopoietic Cell Data
  ref.se <- readRDS("/opt2/db/NovershternHematopoieticData.rds")
  so$Novershtern_main <-  runSingleR(so,ref.se,"label.main")
  so$Novershtern <- runSingleR(so,ref.se,"label.fine")
}

if(species == "mouse"){
  #ImmGen
  ref.se <- readRDS("/opt2/db/ImmGenData.rds")
  so$immgen_main <-  runSingleR(so,ref.se,"label.main")
  so$immgen <- runSingleR(so,ref.se,"label.fine")

  #Curated collection of mouse RNASeq data from GEO
  ref.se <- readRDS("/opt2/db/MouseRNAseqData.rds")
  so$mouseRNAseq_main <-  runSingleR(so,ref.se,"label.main")
  so$mouseRNAseq <- runSingleR(so,ref.se,"label.fine")
}

outputRds=gsub(".rds","_annotated.rds",inputRds)
saveRDS(so,outputRds)

