# Promoter regions
# biowulf location: /data/CCBR_Pipeliner/db/PipeDB/db/Promoters
Definition:
	* GTF file with "gene" in column 3 and "protein_coding" in description
	* 2kb upstream from TSS and 200 bp downstream
GTF versions:
	* Gencode version 30 for hg38
	* Gencode version 19 for hg19
	* Gencode version M1 for mm9
	* Gencode version M21 for mm10

# DNAase Hypersensitivity Regions (DHS)
# biowulf location: /data/CCBR_Pipeliner/db/PipeDB/db/DHS/ENCODE
* mm9 and hg19 versions are downloaded from ENCODE (source: https://www.encodeproject.org/search/?type=Annotation&encyclopedia_version=4&annotation_type=representative+DNase+hypersensitivity+sites)
* overlapping regions are collaped using "bedtools merge"
* mm10 and hg38 versions are then crossmapped from these collapsed versions from the mm10 and hg19 versions, respectively.

# Enhancers 
# biowulf location: /data/CCBR_Pipeliner/db/PipeDB/db/Enhancers/merged
Sources:
	* HACER
		* hg19 : http://bioinfo.vanderbilt.edu/AE/HACER/download/T1.txt
	* FANTOM
		* hg38 : http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz
		* hg19 : http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz
		* mm9 : http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/mouse_permissive_enhancers_phase_1_and_2.bed.gz
		* mm10 : http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.bed.gz
	* DBSuper
               	* mm9 : http://asntech.org/dbsuper/data/bed/mm9/all_mm9_bed.bed
               	* hg19 : http://asntech.org/dbsuper/data/bed/hg19/all_hg19_bed.bed
	* EnhancerAtlas: 
		* fasta files were downloaded using the following commands for hg19 and mm9 data:
"""
for i in 3T3-L1	416B	AtT-20	BAT	Bone_marrow	Brain_E14.5	Brown_preadipocyte_E18.5	C3H10Thalf	CD19+	CD4+CD8+	CD4+Treg	CD4+	CD8+	Cerebellum	Cerebellum_neonate	CH12	CMP	Cortex	Dendritic_cell	EpiLC	EpiSC	Erythroid_fetal_liver	Erythroid_spleen	ESC_Bruce4	ESC_J1	ESC_KH2	ESC_NPC	Forebrain_E11.5	Forebrain_E12.5	Forelimb_bud_embryo	Forelimb_E11	Forelimb_E13	G1E-ER4	G1E	GMP	Heart_E11.5	Heart_E12.5	Heart_E14.5	Heart	Hepatocyte	HFSC	Hindbrain_E11.5	Intestine	IPSC	Kidney	Large_intestine_epithelial	Lens_P1	Limb_E11.5	Limb_E14.5	Liver_E14.5	Liver	Lung_E14.5	Lung	Lung_neonate	MC3T3-E1	Megakaryocyte	MEL	Microglia	Midbrain_E11	Neuron_cortical	NIH-3T3	NKC_spleen	NKT	NPC	Olfactory_bulb	Pancreas	Pancreatic_islet	PDC_BM	PDC	Peritoneal_macrophage	Placenta	Pre-B	Pre-pro-B	Pro-B_BM	Prostate	Rib_chondrocyte_P1	Spermatid	Spleen	Stomach_neonate	Striatum	Testis	Th1	Th2	Thymus	Treg_cell	Uterus	V6.5	WAT	ZHBTc4;do
wget http://www.enhanceratlas.org/data/AllEPs/mm/${i}_EP.txt
done
for f in A549	BJ	Caco-2	CD4+	CD8+	CD14+	CD19+	CD20+	CD34+	CD36+	CD133+	CMK	CUTLL1	ECC-1	GM10847	GM12878	GM12891	GM12892	GM18486	GM18505	GM18507	GM18508	GM18516	GM18522	GM18526	GM18951	GM19099	GM19141	GM19193	GM19238	GM19239	GM19240	H1	  H9	H54	H128	H2171	HCT116	HEK293	HEK293T	Hela	Hela-S3	HepG2	HL-60	  HMEC	HSMM	HUES64	HUVEC	IMR90	Jurkat	K562	Kasumi-1	LNCaP	LoVo	LS174T	  MCF-7	MCF10A	ME-1	MM1S	NB4	NH-A	NHDF	NHEK	NHLF	NT2-D1	OCI-Ly1	P493-6	PANC-1	Raji	SK-N-MC	SK-N-SH	T47D	th1	U2OS	U87	VCaP;do
echo $f;wget http://www.enhanceratlas.org/data/enhseq/${f}.fasta;done
"""
		* chromsome, start, end info is extracted from the fasta file sequence ids and saved as bed files for hg19 and mm9
Processing:
	* overlapping regions are collaped using "bedtools merge" for each source
	* if a version is not available from a particular source then it is created using the collapsed regions file from previous step
	* collapsed regions files from all four sources, HACER/FANTOM/DBSuper/EnhancedAtlas are concatenated using "cat" and collaped again using "bedtools merge"

