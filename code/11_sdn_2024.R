#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - SUDAN (2024) ++++++ ###
#...............................................................................

#...............................................................................
## ----- R SCRIPT TO PREPARE DATA FOR ANALYSIS THROUGH GENERIC SCRIPTS  ----- ##
#...............................................................................


#...............................................................................
### Reading inputs and managing the dataset
#...............................................................................

  #...................................      
  ## Load additional packages needed for the Sudan analysis
  pacman::p_load(
    remotes,       # Install packages from github
    sf,            # Work with spatial data and produce maps
    terra          # Work with geospatial raster datasets
  )
  remotes::install_github("wmgeolab/rgeoboundaries") 
  library("rgeoboundaries")

  
  #...................................      
  ## Read and prepare administrative boundary data
  
    # Download admin1 (state) boundaries
    sdn_boundaries <- geoboundaries("Sudan", "adm1")
  
    # Identify state names and codes
    sdn_adm1 <- sf::st_drop_geometry(sdn_boundaries)
    
  
  #...................................      
  ## Read main dataset
    
    # Read dataset
    df <- as.data.frame(read_excel(
      paste0(dir_path, "in/all_lists_CAM_Sudan_Data_2_Sep.xlsx"),
      sheet = "all_lists_v8"))

    # Fix date
    df$death_date <- as.Date(df$death_date)
    
    # Streamline columns
    x <- grep("PB|PV|SM", colnames(df))
    df <- df[, -x]
    x <- c("cleaned_name", "name", "nickname", "titles", "loc_death_rwd",
      "loc_death_rwd_1", "loc_death_rwd_2", "month_death_rwd", "cod_rwd", 
      "year_death_rwd", "job", "gender_rwd", "hcw", "cod_avoidable", "source", 
      "date_death_rwd", "comments", "comments_2", "publication_date","source_2")
    df <- df[! colnames(df) %in% x]

  #...................................      
  ## Check and manage some variables
    
    # Overlap scores
    for (i in grep("ovrlp_score", colnames(df))) {
      print(i)
      print(table(df[, i], useNA = "ifany"))
    }
    
    # Number repeated
    table(df$n_repeated_list, useNA = "ifany")  

    # Duplication score
    table(df$dup_score, useNA = "ifany")  

    # Harmonise location of death
    table(df$loc_death, useNA = "ifany")
    df[which(df$loc_death == "West Darfour"), "loc_death"] <- "West Darfur"
    df[which(df$loc_death %in% c("NA", "Egypt borders", "Libya borders", 
      "Sudan")), "loc_death"] <- "unknown / unclear"
    df[which(is.na(df$loc_death)), "loc_death"] <- "unknown / unclear"
    table(df$loc_death, useNA = "ifany")
        
    # Cause of death
    table(df$cod, useNA = "ifany")
    df[which(df$cod %in% c("natural death","other","starvation_disease")), 
      "cod"] <- "disease"
    df[which(df$cod == "NA"), "cod"] <- "unknown / unclear"
    df[which(is.na(df$cod)), "cod"] <- "unknown / unclear"
    table(df$cod, useNA = "ifany")

    # Membership of resistance committees   
    table(df$resistance_committees, useNA = "ifany")
    df$resistance_committees <- as.logical(df$resistance_committees)
    table(df$resistance_committees, useNA = "ifany")
        
    # Age of death
    table(df$age, useNA = "ifany") # cannot clean
    table(df$age_cat, useNA = "ifany")
    table(df$age_cat_imp, useNA = "ifany")
    
    # Time of death
      # month  
      table(df$month_death, useNA = "ifany")
      table(df$month_death_imp, useNA = "ifany")
      
      # year    
      table(df$year_death, useNA = "ifany")
      table(df$year_death_imp, useNA = "ifany")
  
      # date (reconstruct)
      x <- ifelse(day(df$death_date) == 15, 15, day(df$death_date))
      df$date_death <- as.Date(paste(df$year_death, df$month_death, x, sep="-"),
        "%Y-%m-%d")
      table(df$date_death, useNA = "ifany")
      df$date_death_imp <- as.Date(paste(df$year_death_imp, df$month_death_imp, 
        x, sep="-"), "%Y-%m-%d")
      table(df$date_death_imp, useNA = "ifany")
      

  #...................................      
  ## Apply some exclusion criteria
    
    # Remove decedents whose del_score is < 5 (insufficient info to analyse)
    table(df$del_score)
    df <- subset(df, del_score == 5)
    
    # Remove decedents who died outside of Sudan
      # variable to exclude because death outside Sudan
      df$excl_loc_death <- F
    
      # based on location of death  
      x <- which(df$loc_death %in% c("Australia", "Chad", "Egypt", "Eritrea",
        "Ethiopia", "Germany", "India", "Kuwait", "Libya", "Mali", "Qatar",
        "Saudi Arabia", "South Sudan", "Tanzania", "Türkiye", "UAE", "Uganda",
        "UK", "USA"))
      df[x, "excl_loc_death"] <- T
          
      # based on date of departure from Sudan
      df$date_dep <- as.Date(paste(df$sudan_leave_year, df$sudan_leave_month, 
        "15", sep = "-"), "%Y-%m-%d")
      x <- which(df$date_dep < df$death_date)
      df[x, "excl_loc_death"] <- T
      
      # tabulate exclusion
      table(df$excl_loc_death)
      df <- df[which(df$excl_loc_death == F), ]
      
##### NEED TO REMOVE DATES OF DEATH OUTSIDE POSSIBLE RANGE, BUT ONLY AFTER
        # GENERATING DATASET FOR c-rc FOR ANY SENSITIVITY SCENARIO
                    
      
#...............................................................................
### ENDS
#...............................................................................
