# BIOT 6900 Coursework

Coursework repository for BIOT 6900: AI-Powered Biotech Workflows.

Grader should theoretically be able to do:

- git clone MY_REPOSITORY_URL
- cd biot6900
- conda env create -f environment.yml
- conda activate biot6900
- jupyter notebook

## Modules

- Module 1: Computational Environment Setup \& Tool Ecosystem
- Module 2: Multi-Omics Target Identification & Validation
   - "biot6900_module2_pdac_data_preparation.ipynb" : Since I selected Pancreatic Ductal Adenocarcinoma data, pulled required genomic, proteomics, transcriptomics data from https://linkedomics.org/data_download/CPTAC-PDAC/
  - "biot6900_module2_pdac_target_ranking.ipynb" : Does the asked work for the homework - ranking the target genes - genes that have all 3 layers of data to do the computation
- Module 3: PDAC target triage with protein language models
   - "biot6900_pdac_module3.ipynb" : completed notebook
   - "targets_pdac.csv": input from module 2 (attached to canvas)
   - "targets_pdac_w3.csv": ranked 15-target output for Week 4 (attached to canvas)



