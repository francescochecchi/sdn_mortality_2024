#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - SUDAN (2024) ++++++ ###
#...............................................................................

#...............................................................................
## ------ R SCRIPT TO LOAD PACKAGES AND SOURCE OTHER ANALYSIS SCRIPTS  ------ ##
#...............................................................................


#...............................................................................
### Preparatory steps
#...............................................................................

  #...................................      
  ## Install or load required R packages
    
    # Install or load packages from CRAN
    pacman::p_load(
      ggplot2,       # Data visualization
      ggpubr,        # Arranging multiple plots into a single plot
      gridExtra,     # Add tables to plot areas
      gtools,        # Assist various programming tasks
      lubridate,     # Makes it easier to work with dates and times
      MASS,          # For various statistical functions
      readxl,        # Read Excel files
      remotes,       # Install packages from github
      scales,        # Scaling and formatting data for visualizations
      sf,            # Work with spatial data and produce maps
      tidyverse,     # Tidyverse suite of packages
      viridis       # Colour palettes
    )

    # Install or load packages not on CRAN
    if (! "rgeoboundaries" %in% rownames(installed.packages())) 
    {remotes::install_github("wmgeolab/rgeoboundaries")}
    library("rgeoboundaries")
  
  #...................................      
  ## Starting setup

    # Clean up from previous code / runs
    rm(list=ls(all=T) )
  
    # Set font for Windows or Mac
    suppressWarnings(windowsFonts(Arial = windowsFont("Arial")))
    suppressWarnings(par(family = "Arial"))

    # Set working directory to where this file is stored
    dir_path <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/")
    setwd(dir_path)
    dir_path <- gsub("/code", "", dir_path)
    
    # Initialise random numbers
    set.seed(123)
    
    # Colour-blind palette for graphing
    palette_gen <- viridis(16)
    show_col(palette_gen)
          

#...............................................................................
### Implement each R script in order
#...............................................................................
  
  #...................................      
  ## Source functions needed for the analysis
  source(paste0(dir_path, "code/01_functions.R"))
    
  #...................................      
  ## Read and prepare datasets for all sensitivity analyses
  source(paste0(dir_path, "code/02_read_prepare_data.R"))
    
  #...................................      
  ## Carry out descriptive analyses
  source(paste0(dir_path, "code/03_describe_data.R"))
    
  #...................................      
  ## Estimate mortality using capture-recapture analysis
  source(paste0(dir_path, "code/04_estimate_mortality.R"))
    
    

#...............................................................................
### ENDS
#...............................................................................
