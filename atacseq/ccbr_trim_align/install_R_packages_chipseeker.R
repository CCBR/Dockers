install.packages("caTools", repos="http://cran.us.r-project.org")
install.packages("XML", repos="http://cran.us.r-project.org")
install.packages("RCurl", repos="http://cran.us.r-project.org")

BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")

