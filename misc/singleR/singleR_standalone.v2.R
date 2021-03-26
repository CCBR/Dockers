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
  so$HPCA_main <- runSingleR(so,HumanPrimaryCellAtlasData(),"label.main")
  so$HPCA <-  runSingleR(so,HumanPrimaryCellAtlasData(),"label.fine")
  
  #Blueprint/ENCODE database
  so$BP_encode_main <-  runSingleR(so,BlueprintEncodeData(),"label.main")
  so$BP_encode <-  runSingleR(so,BlueprintEncodeData(),"label.fine")
  
  #Monaco Immune Cell data
  so$monaco_main <-  runSingleR(so,MonacoImmuneData(),"label.main")
  so$monaco <-	 runSingleR(so,MonacoImmuneData(),"label.fine")
  
  #Database for Immune Cell Expression
  so$dice_main <-  runSingleR(so,DatabaseImmuneCellExpressionData(),"label.main")
  so$dice <- runSingleR(so,DatabaseImmuneCellExpressionData(),"label.fine")
  
  #Novershtern Hematopoietic Cell Data
  so$Novershtern_main <-  runSingleR(so,NovershternHematopoieticData(),"label.main")
  so$Novershtern <- runSingleR(so,NovershternHematopoieticData(),"label.fine")
}

if(species == "mouse"){
  #ImmGen
  so$immgen_main <-  runSingleR(so,ImmGenData(),"label.main")
  so$immgen <- runSingleR(so,ImmGenData(),"label.fine")
  
  #Curated collection of mouse RNASeq data from GEO
  so$mouseRNAseq_main <-  runSingleR(so,MouseRNAseqData(),"label.main")
  so$mouseRNAseq <- runSingleR(so,MouseRNAseqData(),"label.fine")
}

outputRds=gsub(".rds","_annotated.rds",inputRds)
saveRDS(so,outputRds)
