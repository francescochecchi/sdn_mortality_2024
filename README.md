## Population mortality in Sudan, 2023-2024: A capture-recapture analysis
### Description of input dataset and R analysis scripts
October 2024

Funding support: United States Centers for Disease Control and Prevention; United Kingdom Foreign, Commonwealth and Development Office

## General description
This repository contains data and R scripts needed to replicate the above analysis. The two datasets required are found in the `\in` folder, and are read automatically when the code is run. Every row in dataset `sdn_cam_data_final_public.xlsx` contains an individual record of a deceased person, as contained in one of three independent lists generated by the study. The Excel file also contains a variable dictionary.
All the code is found in the `\code` folder. To replicate the analysis, follow these steps:
* Download and unzip the repository to any folder in your computer (other than the Downloads folder, which usually gets wiped automatically). The folder is identified automatically when the code is run.
* Download R and RStudio (see download links on [https://posit.co/download/rstudio-desktop/]). While R is sufficient to run the analysis, it is recommended to instead run the scripts from the RStudio interface.
* Open and run the entire `00_master_code.R` script (just press Alt+Ctrl+R). This will create an `\out` folder with further sub-folders, to which output tables and graphs will be saved automatically. As this scripts calls all the others, it alone is sufficient to replicate the analysis. On a cheap laptop, it should take about 5-10 hours to run the analysis (95% of this time is spent on the simulation in code 05, but the time can be cut down to minutes if the user specifies a lower number of runs (search for the 'n_runs' parameter in script 05 and modify it from 1000 to, say, 10).

## Description of each R script
* `00_master_code.R` installs or loads R packages needed for this analysis, sets general parameters and sources the other scripts in logical order.
* `01_functions.R` contains functions needed for different steps in the analysis.
* `02_read_prepare_data.R` carries out several data cleaning and management tasks, de-duplicates records within each list based on a given duplication score (confidence) threshold, and lastly matches records across lists based on a given match (overlap) confidence score threshold. A clean, analysis-ready dataset is generated for each combination of duplication threshold ('d') and match (overlap) threshold ('o'), and saved in the corresponding `\d[i]o[j]` folder, where i = 1 to 5 and j = 1 to 5, namely the possible levels of the two thresholds.
* `03_describe_data.R` generates various descriptive analyses for each d[i]o[j] dataset, including attrition, characteristics of each list, numbers listed by date, matching deaths by cause and the overlap of deaths. All outputs are saved to each respective `\d[i]o[j]` folder.
* `04_estimate_score_probs.R` performs a simulation to estimate the probability distributions of being a within-list duplicate or a cross-list match, by duplication/match score attributed manually. The output is a set of probability distributions that are then sampled from in the subsequent script to estimate mortality. NOTE: this script reads a private version of the dataset, which contains unique identifiers and is therefore not uploaded to this repository. The script is provided so as to show what was done. The output of the script is also added to the repository so that users can replicate the analysis in full.
* `05_estimate_mortality.R` performs capture-recapture estimation for deaths within Khartoum State only (other states feature too much data sparsity), by cause (all causes and intentional injuries only). All possible models are fitted to the three-list matching structure of the dataset, confounders are added to the model and lastly the alternative models are averaged using their goodness-of-fit as a weight, yielding estimates of the unlisted and total number of deaths, with confidence intervals, by model and overall. Results are saved in different .csv files, where the 'raw' version contains more information, and the 'pretty' version is meant for publication. The sensitivity of each list is also estimated. The main analysis relies on estimated probabilities of duplication and matching, and its results are stored in the `\out` folder. Outputs of each sensitivity analysis are stored in the corresponding `\d[i]o[j]` folder. Tables and graphs summarising all the different sensitivity analyses ('d', 'o' combinations) are also generated: these are stored in the `\out` folder.

