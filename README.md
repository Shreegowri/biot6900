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

########################### START OF PART I ###########################

- Module 2: Multi-Omics Target Identification & Validation
   - "biot6900_module2_pdac_data_preparation.ipynb" : Since I selected Pancreatic Ductal Adenocarcinoma data, pulled required genomic, proteomics, transcriptomics data from https://linkedomics.org/data_download/CPTAC-PDAC/
  - "biot6900_module2_pdac_target_ranking.ipynb" : Does the asked work for the homework - ranking the target genes - genes that have all 3 layers of data to do the computation
- Module 3: PDAC target triage with protein language models
   - "biot6900_pdac_module3.ipynb" : completed notebook
   - "targets_pdac.csv": input from module 2 (attached to canvas)
   - "targets_pdac_w3.csv": ranked 15-target output for Week 4 (attached to canvas)
- Module 4: PDAC structure-based target assessment
   - "biot6900_pdac_module4.ipynb" : completed notebook for the structure assessment
   - "targets_pdac.csv" & "targets_pdac_w3.csv": input from previous modules
   - "results": Results folder for structure assessment of teh selected proteins and fpocket - uploaded to canvas - this would be produced by the notebook
   - "portfolio/figures": Produced by "plot_portfolio_figures.py" - script generating images used in the final report

########################### END OF PART I ###########################

- Module 5: DPYD variant interpretation
   - "biot6900_module5_assignment.ipynb" : completed notebook
 - Module 6: PRS portability
   - "module6_part3_assignment.ipynb" : completed notebook
   - input used from professor's data.zip file at canvas

