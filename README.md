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

I selected pancreatic ductal adenocarcinoma (PDAC) because it is an aggressive cancer with substantial unmet need.
I used processed, publicly downloadable CPTAC-PDAC matrices from LinkedOmics (https://linkedomics.org/data_download/CPTAC-PDAC/). So, it is OPEN data.

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

########################### START OF PART II ###########################
- Module 5: DPYD variant interpretation
   - "biot6900_module5_assignment.ipynb" : completed notebook

 - Module 6: PRS portability
   - "module6_part3_assignment.ipynb" : completed notebook
   - input used from professor's data.zip file at canvas

- Module 7: PDAC knowledge graph
     - "biot6900_module7_assignment.ipynb" : completed notebook
          - Updated biot6900_env.yml with new packages needed to run this notebook
          - Need Aura credentials to run this notebook - which I have kept in my .env file - the username, password, database name of the NEO4J - this was in credential file I recieved from Aura
     - "data.zip" - uploaded to canvas - used the files kept inside data folder as input for the module 7 excercise
     - "results.zip" - exported csv files with nodes, edges, relationship, centrality score etc., are kept in this folder and uplaoded to canvas

- Module 8: Single-Cell & Spatial Omics with AI
    - Updated python version to 3.12 as squidpy requires later version of python, 3.11 does not work
    - Target gene choosen : HER2 (ERBB2)
    - The cache contains the original pre-QC count matrix used for the report. It was exported from the CZ CELLxGENE Census release 2025-11-08 using the primary breast cancer filter in the notebook. The sample contains 2,000 cells from each of 10 donors. The notebook verifies the release, filter and fingerprint of the selected cell IDs before analysis. On a supported system, the notebook can repeat the Census query when the cache is absent.
    - "biot6900_module8_assignment.ipynb" : completed notebook with results

- Module 9: Clinical Data Analysis & Natural Language Processing
  - "module9_part3_assignment.ipynb" : completed notebook with results and my comments
       - The notebook uses numpy, pandas, scikit-learn, and matplotlib, plus the standard library. It needs no downloaded data or clinical NLP package.

- Module 10: Medical Imaging AI
   - "biot6900_module10_assignment.ipynb" : completed notebook with results
- Module 11:
########################### END OF PART II ###########################
