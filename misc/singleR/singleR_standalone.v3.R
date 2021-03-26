#optimized for R/3.6.1
#for details: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html

library(S4Vectors)
library(GenomicRanges)
library(SummarizedExperiment)
library(IRanges)
library(GenomeInfoDb)
library(Biobase)
library(DelayedArray)
library(matrixStats)
library(SingleCellExperiment)
library(SingleR)
library(Seurat)

args<-commandArgs(trailingOnly=T)

inputRds = as.character(args[1])
species = as.character(args[2])

so=readRDS(inputRds)

runSingleR = function(obj,refFile,fineORmain){
  #convert Seurat object to SingleCellExperiment object, as required by SingleR
  #use SCTransform assay if available. Otherwise use base RNA counts
  
  if(length(grep("SCT",Assays(obj))!=0)){
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
  ref.se <- readRDS("HumanPrimaryCellAtlasData.rds")
  so$HPCA_main <- runSingleR(so,ref.se,"label.main")
  so$HPCA <-  runSingleR(so,ref.se,"label.fine")
  
  #Blueprint/ENCODE database
  ref.se <- readRDS("BlueprintEncodeData.rds")
  so$BP_encode_main <-  runSingleR(so,ref.se,"label.main")
  so$BP_encode <-  runSingleR(so,ref.se,"label.fine")
  
  #Monaco Immune Cell data
  ref.se <- readRDS("MonacoImmuneData.rds")
  so$monaco_main <-  runSingleR(so,ref.se,"label.main")
  so$monaco <-	 runSingleR(so,ref.se,"label.fine")
  
  #Database for Immune Cell Expression
  ref.se <- readRDS("DatabaseImmuneCellExpressionData.rds")
  so$dice_main <-  runSingleR(so,ref.se,"label.main")
  so$dice <- runSingleR(so,ref.se,"label.fine")
  
  #Novershtern Hematopoietic Cell Data
  ref.se <- readRDS("NovershternHematopoieticData.rds")
  so$Novershtern_main <-  runSingleR(so,ref.se,"label.main")
  so$Novershtern <- runSingleR(so,ref.se,"label.fine")
}

if(species == "mouse"){
  #ImmGen
  ref.se <- readRDS("ImmGenData.rds")
  so$immgen_main <-  runSingleR(so,ref.se,"label.main")
  so$immgen <- runSingleR(so,ref.se,"label.fine")
  
  #Curated collection of mouse RNASeq data from GEO
  ref.se <- readRDS("MouseRNAseqData.rds")
  so$mouseRNAseq_main <-  runSingleR(so,ref.se,"label.main")
  so$mouseRNAseq <- runSingleR(so,ref.se,"label.fine")
}

outputRds=gsub(".rds","_annotated.rds",inputRds)
saveRDS(so,outputRds)
