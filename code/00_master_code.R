#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - GENERIC CODE ++++++ ###
#...............................................................................

#...............................................................................
## ------ R SCRIPT TO LOAD PACKAGES AND SOURCE OTHER ANALYSIS SCRIPTS  ------ ##
#...............................................................................


#...............................................................................
### Preparatory steps
#...............................................................................

  #...................................      
  ## Install or load required R packages
  pacman::p_load(
    ggplot2,       # Data visualization
    ggpubr,        # Arranging multiple plots into a single plot
    gtools,        # Assist various programming tasks
    lubridate,     # Makes it easier to work with dates and times
    MASS,          # For various statistical functions
    readxl,        # Read Excel files
    scales,        # Scaling and formatting data for visualizations
    tidyverse,     # Tidyverse suite of packages
    viridis        # Colour palettes
  )

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
    print(getwd() )
    
    # Initialise random numbers
    set.seed(123)
    
    # Colour-blind palette for graphing
    palette_gen <- viridis(16)
    show_col(palette_gen)
          

#...............................................................................
### Reading in required inputs
#...............................................................................

  #...................................      
  ## Identify data file and read in parameters
    
    # Search for data file in the directory based on string pattern
    filename <- list.files(paste0(dir_path, "in"), pattern = "list_data")[1]

    # Read in parameters from Excel file
    pars <- read_excel(paste0(dir_path, "in/", filename), sheet = "parameters")
    pars <- as.data.frame(pars)

  #...................................      
  ## Read in mortality datasets
  
    # Create names for datasets of each location (also used for analysis output)
    x <- gsub("[^[:alnum:] ]", "", pars$location) # remove non-alphanumeric char
    x <- trimws(gsub("\\s+", " ", x)) # trim whitespace
    x <- tolower(gsub(" ", "_", x)) # names to lower case, spaces to underscores
    pars$loc_df <- x

    # Read individual datasets for each location and assign objects
        # only read in datasets that should be analysed
    for (i in 1:nrow(pars)) {
      if (pars[i, "analyse"] == "Y") {
        
        # read the sheet corresponding to the location name
        x <- read_excel(paste0(dir_path, "in/", filename), 
          sheet = pars[i, "location"])
        x <- as.data.frame(x)
        
        # name the dataframe after the location
        assign(pars[i, "loc_df"], x)
      }
    }
    
    
  #...................................      
  ## Source functions
  source(paste(dir_path, "code/01_functions.R"), echo = T)

#...............................................................................
### Preparing data for analysis
#...............................................................................

  #...................................      
  ## Recognise characteristics of each dataset
  
    # Calculate the number of lists for each dataset
    pars$n_lists <- apply(pars[, paste0("list", 1:4)], 1, 
      function(x) {length(which(!is.na(x)))})
  
    # Subset only parameters for datasets that should be analysed
    pars_in <- subset(pars, analyse == "Y")

    
  #...................................      
  ## Manage each dataset
  for (i in pars_in$loc_df) { 
    
    # Clean dataset
    assign(i, f_clean(get(i), i))

    # Prepare datasets for each analysis stratum as specified by the user
    assign(i, f_stratify(get(i), i))
  }

    
    
#...............................................................................
### Visualising data and preparing contingency tables for analysis
#...............................................................................

  #...................................      
  ## Describe patterns and overlap among lists
  for (i in pars_in$loc_df) { 
    
    # Describe patterns in period, age and gender, by list and overall
    f_describe(get(i), i)
    
    # Visualise and quantify overlap among lists
    assign(paste0(i, "_overlap"), f_overlap(get(i), i))
  }

    
    
#...............................................................................
### Analysing data for each location (site)
#...............................................................................

for (i in pars_in$loc_df) { 

  #...................................      
  ## Preparatory steps: assemble inputs, initialise output

    # Print a message indicating which site is being analysed
    print(paste0("now doing analysis for this site: ", 
      pars_in[pars_in$loc_df == i, "location"]))
    
    # Assemble required inputs
      # dataset
      df <- get(i) 
      
      # output of preparatory steps for this dataset
      overlap <- get(paste0(i, "_overlap"))
      
      # number of lists
      n_lists <- pars_in[pars_in$loc_df == i, "n_lists"]
    
    # Initialise output table
      
      # output file name  
      filename <- paste0(dir_path, "out/", i, "_analysis_output.csv")
    
      # write the location of analysis
      write.table(paste0("LOCATION: ", 
        pars_in[pars_in$loc_df == i, "location"]), filename, 
        row.names = F, col.names = F, sep = ",")
    
      # write the period of analysis
      write.table(paste0("  Period: ", min(subset(df, eligible == T)$date_clean, 
        na.rm = T), " to ", 
        max(subset(df, eligible == T)$date_clean, na.rm = T)), 
        filename, append = T, row.names = F, col.names = F, sep = ",")
      
      # write a separator
      write.table(rbind(rep("......................", 5)), filename, 
        append = T, row.names = F, col.names = F, sep = ",")

  #...................................      
  ## Fit candidate models for both all observations and each desired stratum
  for (j in 1:nrow(overlap)) {
    
    # Control message: which stratum is being analysed?
    print(paste0("  working on this stratum: ", 
      gsub("_", " ", overlap[j, "stratum"])))

    # Continue table
    write.table(paste0("Stratum: ", gsub("_", " ", overlap[j, "stratum"])), 
      filename, append = T, row.names = F, col.names = F, sep = ",")

    # If two lists only, simple Chapman estimator
    if (n_lists == 2) {out <- f_chapman(overlap[j, ], i, pars_in)}

    # If three or four lists, log-linear models and model averaging
    if (n_lists %in% c(3, 4)) {
      
      # if analysis is on all observations...
      if (j == 1) {df_j <- df; pars_j <- pars_in}
      
      # if analysis is on a single stratum...
      if (j > 1) {
        
        # subset only stratum of interest
        df_j <- df[which(df[, overlap[j, "variable"]] ==overlap[j, "stratum"]),] 
        
        # check that the stratum isn't also specified as a confounder...
        pars_j <- pars_in
        x1 <- which(pars_j$loc_df == i)
        x2 <- trimws(unlist(strsplit(pars_j[x1, "confounders"], ",")))
        x3 <- gsub("_stratum", "", overlap[j, "variable"])
        if (x3 %in% x2) {pars_j[x1, "confounders"] <- paste0(x2[x2 != x3])}

        #... or an exposure
        x2 <- pars_j[x1, "exposure"]
        x3 <- overlap[j, "variable"]
        if (x3 %in% x2) {pars_j[x1, "exposure"] <- NA}
      }
    
      # implement log-linear models
      x <- f_logl(data_f = df_j, data_name_f = i, pars_f = pars_j)
  
      # perform model averaging
      out <- f_model_average(f_out = x, data_name_f = i, pars_f = pars_j)
    }
    
    # Write results to table
      # 2-list scenario
      if (n_lists == 2) {
        write.table(out, filename, append = T, row.names = F, 
          col.names = T, sep = ",")
        write.table(rbind(rep("......................", 2)), filename, 
          append = T, row.names = F, col.names = F, sep = ",")
      }        
      
    # 3- or 4-list scenario
    if (n_lists %in% c(3, 4)) {
      
      # unformatted output
      write.table("  Unformatted output:", filename, append = T, 
        row.names = F, col.names = F, sep = ",")      
      write.table("    Candidate models:", filename, append = T, 
        row.names = F, col.names = F, sep = ",")
      write.table(out$out_raw, filename, append = T, row.names = F, 
        na = "", sep = ",")
      write.table("    Estimated deaths:", filename, append = T, 
        row.names = F, col.names = F, sep = ",")
      write.table(out$out_est_raw, filename, append = T, 
        row.names = F, col.names = T, na = "", sep = ",")
      write.table("    List sensitivity:", filename, append = T, 
        row.names = F, col.names = F, sep = ",")
      write.table(out$out_sens_raw, filename, append = T, 
        row.names = F, na = "", sep = ",")
      write.table(rbind(rep("-----", 2)), filename, append = T, 
        row.names = F, col.names = F, sep = ",")
      
      # formatted output
      write.table("  Formatted output:", filename, append = T, 
        row.names = F, col.names = F, sep = ",")      
      write.table("    Candidate models:", filename, append = T, 
        row.names = F, col.names = F, sep = ",")
      write.table(out$out_pretty, filename, append = T, 
        row.names = F, na = "", sep = ",")
      write.table("    Estimated deaths:", filename, append = T, 
        row.names = F, col.names = F, sep = ",")
      write.table(out$out_est_pretty, filename, append = T, 
        row.names = F, col.names = T, na = "", sep = ",")
      write.table("    List sensitivity:", filename, append = T, 
        row.names = F, col.names = F, sep = ",")
      write.table(out$out_sens_pretty, filename, append = T, 
        row.names = F, na = "", sep = ",")
      write.table(rbind(rep("......................", 5)), filename, 
        append = T, row.names = F, col.names = F, sep = ",")
    }
  }
}

#...............................................................................
### ENDS
#...............................................................................
