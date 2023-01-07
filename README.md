# Replication Materials for "Uses and Limitations of Dichotomous Aggregate Relational Data"

## Software Requirement

To undertake replication, you need R and Stan installed on your computer. Required packages are listed at the top of the R scripts.  

## Data

The folder "Dichotomous-ARD-Replication-Materials/Data" contains the data sets already transformed and recoded from the source versions. Please refer to the "Analytic Strategy" segment of the paper for directions to the raw data and steps for preparing them for analysis. 

## Replication

Below are steps for replicating the results presented in the main text. To reproduce the results in the appendix, follow the same steps using the files stored in the "Appendix" subfolders of the "Data," "Run Stan Models," and "Replicate Results" folders.

  1. Save the Stan files stored in "Dichotomous-ARD-Replication-Materials/Stan Files/". This folder contains five subfolders corresponding to the models presented in the paper. Each subfolder, in turn, includes versions of the models for different data treatments (count, dichotomous, or 0/1/2+). 
  2. Go to "Dichotomous-ARD-Replication-Materials/Run Stan Models/Main Text/". Follow the headers and subheaders of "run_stan_maintext.R" to run the Stan models of Step 1. 
  3. Save the objects storing the extracted samples from each model. Their names end with "Extract" in "run_stan_maintext.R" (e.g., mccartyCount_MaltielRDM_Extract).
  4. Using the extracted samples of Step 3 and the data files stored in "Dichotomous-ARD-Replication-Materials/Data/Main Text/", follow the headers and subheaders of "replicate_results_maintext.R" in "Dichotomous-ARD-Replication-Materials/Replicate Results/Main Text/". They provide directions for reproducing each table and figure presented in the paper.
