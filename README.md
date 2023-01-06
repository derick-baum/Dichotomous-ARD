# Replication Materials for "Uses and Limitations of Dichotomous Aggregate Relational Data"

## Data

The data folder contains the data sets used in the paper, already transformed and recoded from the source versions. Please refer to the "Analytic Strategy" segment for directions to the raw data and steps for preparing them for analysis. 

## Replication

Below I present steps for replicating the results presented in the main text. To reproduce the results in the appendix, follow the same steps using the files stored in the "Appendix" subfolders of the "Data," "Run Stan Models," and "Replicate Results" folders.

  1. Download the Stan files stored in "Dichotomous-ARD-Replication-Materials/Stan Files/". This folder contains five subfolders corresponding to models presented in the paper. Each subfolder, in turn, includes versions of the models for different data treatments (count, dichotomous, or 0/1/2+). 
  2. Go to "Dichotomous-ARD-Replication-Materials/Run Stan Models/Main Text/". Follow the headers and subheaders of "run_stan_maintext.R" to run the Stan models of Step 1. 
  3. Save the objects storing the extracted samples from each model. Their names end with "Extract" (e.g., mccartyCount_MaltielRDM_Extract) in "run_stan_maintext.R".
  4. Download the R script "replicate_results_maintext.R" from "Dichotomous-ARD-Replication-Materials/Replicate Results/Main Text/". Using the extracted samples of Step 3 and the data files stored in "Dichotomous-ARD-Replication-Materials/Data/Main Text/", follow the headers and subheaders of "run_stan_maintext.R". They provide guidance for reproducing each table and figure presented in the paper.
