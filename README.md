# Roborovski dwarf hamster scRNA-seq (Dexamethasone/antibody treatments, lung samples)
Single-cell RNA-sequencing 

Raw and processed files are available at NCBI GEO, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191080

1. Starting with raw data
The twelve samples (three animals each for untreated, antibody, dexamethasone, and antibody+dexamethasone treatments) were sequenced in three pools, and every pool was sequenced twice. Therefore, via GEO there are 4 SRR entries, two each for the gene expression and the cellplex sequencing. The two sequencings can be placed in separate folders, and the cellranger multi command can by run using the multiconfig.csv file provided here.
For the Phodopus roborovskii reference folder, the fasta and annotation can be downloaded here: https://figshare.com/articles/dataset/Phodopus_roborovskii_assembly/16695457 – the main part was done with the gtf file phorob_curated_rev2_S2.gtf in the customized_annotations.zip folder. The only difference of phorob_curated_rev3_S2.gtf compared to phorob_curated_rev2_S2.gtf is an annotation of the Nr3c1 gene (glucocorticoid receptor). The reference using the fasta and gtf file can then be done using the cellranger mkref command. And overview of the data can be found in the single-cell_metadata.txt.gz file at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191080.

2. Initial processing of single-cell RNA-seq data
The data is processed using the the merge_integrate.R script, which yields the DexAb_combined_integrated.rds file that is then the input for the R scripts in PhoRob_DexAb_scRNAseq.R used for data analysis and figure imaging. At the GEO entry, the file DexAb_combined_integrated_b.rds can be found which also contains the Nr3c1 gene.

3. Overview of processed files in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191080
The h5 files of the individual samples are in the RAW.tar file. Note that the data from the three pools mentioned above are not merged here, i.e. for every sample there are three h5 files (e.g. Hamster1_AB1_sample_feature_bc_matrix.h5, Hamster2_AB1_sample_feature_bc_matrix.h5 and Hamster3_AB1_sample_feature_bc_matrix.h5 for antibody sample 1). DexAb.int.rds.gz contains the fully annotated dataset used/generated in PhoRob_DexAb_scRNAseq. DexAb.neu.rds.gz contains Neutrophil dataset.
