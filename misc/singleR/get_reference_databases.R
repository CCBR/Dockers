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

// HumanPrimaryCellAtlasData.ref.se = HumanPrimaryCellAtlasData()
// saveRDS(HumanPrimaryCellAtlasData.ref.se,"HumanPrimaryCellAtlasData.rds")
// 
// BlueprintEncodeData.ref.se = BlueprintEncodeData()
// saveRDS(BlueprintEncodeData.ref.se, "BlueprintEncodeData.rds")
// 
// MonacoImmuneData.ref.se = MonacoImmuneData()
// saveRDS(MonacoImmuneData.ref.se, "MonacoImmuneData.rds")
// 
// DatabaseImmuneCellExpressionData.ref.se = DatabaseImmuneCellExpressionData()
// saveRDS(DatabaseImmuneCellExpressionData.ref.se, "DatabaseImmuneCellExpressionData.rds")
// 
// NovershternHematopoieticData.ref.se = NovershternHematopoieticData()
// saveRDS(NovershternHematopoieticData.ref.se,"NovershternHematopoieticData.rds")
// 
// ImmGenData.ref.se = ImmGenData()
// saveRDS(ImmGenData.ref.se,"ImmGenData.rds")
// 
// MouseRNAseqData.ref.se = MouseRNAseqData()
// saveRDS(MouseRNAseqData.ref.se,"MouseRNAseqData.rds")

readRDS("HumanPrimaryCellAtlasData.rds")
readRDS("BlueprintEncodeData.rds")
readRDS("MonacoImmuneData.rds")
readRDS("DatabaseImmuneCellExpressionData.rds")
readRDS("NovershternHematopoieticData.rds")
readRDS("ImmGenData.rds")
readRDS("MouseRNAseqData.rds")

save(HumanPrimaryCellAtlasData.ref.se,BlueprintEncodeData.ref.se,MonacoImmuneData.ref.se,DatabaseImmuneCellExpressionData.ref.se,NovershternHematopoieticData.ref.se,ImmGenData.ref.se,MouseRNAseqData.ref.se, file="singleR_databases_workspace.RData")