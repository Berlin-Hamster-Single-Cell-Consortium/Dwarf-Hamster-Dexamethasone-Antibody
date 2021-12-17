# Info

This is the python code for the RNA velocity analysis and diffusion analysis of the
dwarf hamster scRNA-seq data (PhoRob_DexAb). The loom file is available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191080.

# Install & Snakemake

We use snakemake to make this code reproducible. To get it to work do the following:
- install python and the required packages with conda/pip (see "requirements.txt")
- download the required loom data file from GEO
- go to "Snakefile" and change the paths where I left a NOTE to your paths
- run the command "snakemake" to produce all data files

# Notebooks
In the folder "notebooks/" you can find the jupyter notebook "PhoRob_DexAb_Figures.ipynb", which we
used to produce the figures in the paper from the snakemake generated files.
