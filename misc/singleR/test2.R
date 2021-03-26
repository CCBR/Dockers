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
  #Novershtern Hematopoietic Cell Data
  ref.se <- readRDS("/opt2/db/NovershternHematopoieticData.rds")
  print("#Novershtern Hematopoietic Cell Data READ DONE!")
  so$Novershtern <- runSingleR(so,ref.se,"label.fine")
  print("#Novershtern Hematopoietic Cell Data FINE DONE!")
  so$Novershtern_main <-  runSingleR(so,ref.se,"label.main")
  print("#Novershtern Hematopoietic Cell Data MAIN DONE!")
  
  #Human Primary Cell Atlas
  ref.se <- readRDS("/opt2/db/HumanPrimaryCellAtlasData.rds")
  print("#Human Primary Cell Atlas DONE!")
  so$HPCA_main <- runSingleR(so,ref.se,"label.main")
  print("#Human Primary Cell Atlas DONE!")
  so$HPCA <-  runSingleR(so,ref.se,"label.fine")
  print("#Human Primary Cell Atlas DONE!")

  #Blueprint/ENCODE database
  ref.se <- readRDS("/opt2/db/BlueprintEncodeData.rds")
  print("#Blueprint/ENCODE database DONE!")
  so$BP_encode_main <-  runSingleR(so,ref.se,"label.main")
  print("#Blueprint/ENCODE database DONE!")
  so$BP_encode <-  runSingleR(so,ref.se,"label.fine")
  print("#Blueprint/ENCODE database DONE!")

  #Monaco Immune Cell data
  ref.se <- readRDS("/opt2/db/MonacoImmuneData.rds")
  print("#Monaco Immune Cell data DONE!")
  so$monaco_main <-  runSingleR(so,ref.se,"label.main")
  print("#Monaco Immune Cell data DONE!")
  so$monaco <-	 runSingleR(so,ref.se,"label.fine")
  print("#Monaco Immune Cell data DONE!")

  #Database for Immune Cell Expression
  ref.se <- readRDS("/opt2/db/DatabaseImmuneCellExpressionData.rds")
  print("#Database for Immune Cell Expression DONE!")
  so$dice_main <-  runSingleR(so,ref.se,"label.main")
  print("#Database for Immune Cell Expression DONE!")
  so$dice <- runSingleR(so,ref.se,"label.fine")
  print("#Database for Immune Cell Expression DONE!")

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
