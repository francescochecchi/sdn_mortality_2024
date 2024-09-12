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
    sf            # Work with spatial data and produce maps
    # terra          # Work with geospatial raster datasets
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
      paste0(dir_path, "in/sdn_cam_data_12sep2024.xlsx"), sheet = "data"))

    # Fix date
    df$death_date <- as.Date(df$death_date)
    
    # Streamline columns
    x <- c("cleaned_name", "name", "nickname", "titles", "loc_death_rwd", "age",
      "loc_death_rwd_1", "loc_death_rwd_2", "month_death_rwd", "cod_rwd", 
      "year_death_rwd", "job", "gender_rwd", "hcw", "cod_avoidable", "source", 
      "date_death_rwd", "comments", "comments_2", "publication_date","source_2")
    df <- df[! colnames(df) %in% x]

    # Remove decedents whose del_score is < 5 (insufficient info to analyse)
    table(df$del_score)
    df <- subset(df, del_score == 5)

  #...................................      
  ## Fix overlap ID variables so they only contain IDs of possible matches

    # Check uniqueness of ID and duplicate ID
    table(table(df$id))
    table(table(df$dup_id))
    
    # Remove end suffixes ("_1", "_2") that denote mirror sets of overlap pairs
    x <- grep("ovrlp_id", colnames(df), value = T)
    for (i in x) {df[, i] <- substr(df[, i], 1, nchar(df[, i]) - 2) }

    # Remove "_" before/after each list acronym
    for (i in c("PB", "PV", "SM")) {
      df$id <- gsub(paste0(i, "_"), i, df$id)
      df$dup_id <- gsub(paste0(i, "_"), i, df$dup_id)
      for (j in x) {
        df[, j] <- gsub(paste0("_", i), i, df[, j])
        df[, j] <- gsub(paste0(i, "_"), i, df[, j])
      }
    }
    
    # Remove same-row duplicate ID from ovrlp IDs
    for (i in 1:nrow(df)) {
      for (j in x) {df[i, j] <- gsub(df[i, "dup_id"], "", df[i, j])}
    }

  #...................................      
  ## Create, check and manage other variables
    
    # 'Parent' ID that each duplicate belongs to
    df$parent_id <- sapply(strsplit(df$dup_id, "_"), function(x) {x[1]})

    # List that each observation belongs to
    df$list_name <- substr(df$id, 1, 2)
    table(df$list_name, useNA = "ifany")
    df$list <- NA
    df[which(df$list_name == "PB"), "list"] <- "list1"
    df[which(df$list_name == "PV"), "list"] <- "list2"
    df[which(df$list_name == "SM"), "list"] <- "list3"
    table(df$list, useNA = "ifany")
    
    # Overlap scores
    for (i in grep("ovrlp_score", colnames(df))) {
      print(i)
      print(table(df[, i], useNA = "ifany"))
    }
    
    # Number repeated
    table(df$n_repeated_list, useNA = "ifany")
    df$n_rep <- na.replace(df$n_repeated_list, 1)
    df[which(df$n_rep == 0), "n_rep"] <- 1
    table(df$n_rep, useNA = "ifany")
#### THINK ABOUT THIS MORE
        
    # Duplication score
    table(df$dup_score, useNA = "ifany")  

    # Harmonise location of death
    table(df$loc_death, useNA = "ifany")
    df[which(df$loc_death == "West Darfour"), "loc_death"] <- "West Darfur"
    df[which(df$loc_death %in% c("NA", "Egypt borders", "Libya borders", 
      "Sudan")), "loc_death"] <- "unknown / unclear"
    df[which(is.na(df$loc_death)), "loc_death"] <- "unknown / unclear"
    x <- c("Australia", "Chad", "Egypt", "Eritrea",
        "Ethiopia", "Germany", "India", "Kuwait", "Libya", "Mali", "Qatar",
        "Saudi Arabia", "South Sudan", "Tanzania", "Türkiye", "UAE", "Uganda",
        "UK", "USA")
    df[which(df$loc_death %in% x), "loc_death"] <- "outside Sudan"
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
    table(df$age_cat, useNA = "ifany")
    table(df$age_cat_imp, useNA = "ifany")
    
    # Date of death
      # Specify start and end of analysis / data collection period  
      date_start <- as.Date(paste(2023, 4, 15, sep = "-"), "%Y-%m-%d")
      date_end <- as.Date(paste(2024, 6, 4, sep = "-"), "%Y-%m-%d")
      
      # day
      df$day_death <- ifelse(is.na(df$death_date), 15, day(df$death_date))
      table(df$day_death, useNA = "ifany")
      df$day_death_imp <- df$day_death
      table(df$day_death_imp, useNA = "ifany")
    
      # month  
      table(df$month_death, useNA = "ifany")
      table(df$month_death_imp, useNA = "ifany")
      
      # year    
      table(df$year_death, useNA = "ifany")
      table(df$year_death_imp, useNA = "ifany")
   
    # Date of leaving Sudan
      # month
      table(df$sudan_leave_month, useNA = "ifany")
      df[which(df$sudan_leave_month == "unknown"), "sudan_leave_month"] <- NA
      df$sudan_leave_month <- as.integer(df$sudan_leave_month)
      
      # year
      table(df$sudan_leave_year, useNA = "ifany")
      df$sudan_leave_year <- as.integer(df$sudan_leave_year)
      
      
#...............................................................................
### Generating datasets for each sensitivity scenario
#...............................................................................

  #...................................      
  ## Define parameters
    
    # Which sensitivity analyses
    sens <- expand.grid(dup_score = 5:0, ovrlp_score = 1:5) 
    sens$sens <- paste0("d", sens$dup_score, "_o", sens$ovrlp_score)
    
    # Variables that need to be merged when there is overlap, and how  
    vars_merge <- data.frame(
      var = c("list", "list_name", "parent_id", "n_rep", "loc_death", "gender", 
        "age_cat", "age_cat_imp", "month_death", "day_death", "day_death_imp",
        "month_death_imp", "year_death", "year_death_imp",  
        "resistance_committees", "cod", "sudan_leave_month","sudan_leave_year"),
      how = c("mode", "mode", "mode", "mean", "mode", "mode", 
        "mode", "mode", "mean", "mean", "mean",
        "mean", "mean", "mean", 
        "mode", "mode", "mean", "mean")
    )
      # if how == mode and values are discordant, take the most common of 3, 
        # the sole value out of 2/3 or, if 2 values out of 2 are discordant, 
        # set to NA

    # List all possible overlap pairs
    ovrlp <- c()
    for (i in 1:5) {
      x <- na.omit(df[, c("id", "dup_id", "parent_id", 
        paste0(c("ovrlp_id_", "ovrlp_score_"), i))])
      colnames(x) <- c("id", "dup_id", "parent_id", "ovrlp_id", "ovrlp_score")
      ovrlp <- rbind(ovrlp, x)
    }
      
    # Identify all parent observations that have 1 or more possible duplicates
    x <-  by(df, df$parent_id, nrow)
    x <- names(x[which(x > 1)])
    df$parent <- ifelse(df$dup_id %in% x, T, F)    
    table(df$parent, useNA = "ifany")


  
#...............................................................................
### Function for first merging duplicates and returning data with 
      # dates reconstructed from merged values
#...............................................................................

f_dup <- function(df_f = df, vars_merge_f = vars_merge, threshold_dup = 3,
  date_start_f = date_start, date_end_f = date_end) {
    
  #...................................      
  ## Prepare for merging
    
    # Initialise a generic one-row dataframe to hold single merged row
    single <- as.data.frame.matrix(matrix(NA, ncol=nrow(vars_merge_f), nrow=1))
    colnames(single) <- vars_merge_f$var

    # Apply duplicate threshold and identify the observations that need to merge
    x <- which((df_f$dup_score >= threshold_dup) | df_f$parent)
    df_merge <- df_f[x, ]
    
    # Set aside observations that don't need to merge aside
    df_nomerge <- df_f[-x, ]

  #...................................      
  ## Merge duplicates
  out_merge <- by(df_merge, df_merge$parent_id, function(xx) {

    # Initialise merged row
    merged <- single
    
    # Merge variables by mode
    vars_mode <- vars_merge_f[which(vars_merge_f$how == "mode"), "var"]
    for (i in vars_mode) {
      
      # compute mode(s)
      tab <- table(xx[, i])
      mode_x <- suppressWarnings(names(tab[which(tab == max(tab))]))
      
      # if a single most common value, return it, else return NA
      merged[, i] <- ifelse(length(mode_x) == 1, mode_x, NA)
    }

    # Merge variables by mean
    vars_mean <- vars_merge_f[which(vars_merge_f$how == "mean"), "var"]
    merged[, vars_mean] <- colMeans(xx[, vars_mean], na.rm = T)
    
    # Return
    return(merged)
  })
  out_merge <- do.call(rbind, out_merge)

  #...................................      
  ## Return prepared dataset

    # Update dataset after merging duplicates
    df_out <- rbind(df_nomerge[,vars_merge_f$var], out_merge[,vars_merge_f$var])
    
    # Reconstruct date of death based on merged values
    x <- grep("day|month|year", colnames(df_out))
    for (i in x) {df_out[, x] <- round(df_out[, x], 0)}
    df_out$date_death <- as.Date(paste(df_out$year_death, df_out$month_death, 
      df_out$day_death, sep="-"), "%Y-%m-%d")
    df_out$date_death_imp <- as.Date(paste(df_out$year_death_imp, 
      df_out$month_death_imp, df_out$day_death_imp, sep="-"), "%Y-%m-%d")
      
      # set dates to NA if they are after end of data collection
      df_out[which(df_out$date_death > date_end_f), "date_death"] <- NA
      df_out[which(df_out$date_death_imp > date_end_f), "date_death_imp"] <- NA

    # Return
    return(df_out)
}

    
    
###### n_rep: sum or mean/mode? CHECK AGAIN

    

##### ANY POSSIBLE DUPLICATES WITH POSSIBLE OVERLAP WITH DIFFERENT PEOPLE?
      # UNLIKELY (MAYBE CHECK LATER)
    
        
##### FIRST DO OVERLAP, THEN SPLIT BY STATE    
    
      # generate exclusion criterion for dates before analysis period
      df_out$excl_date <- ifelse(df_out$date_death < date_start_f, T, F)
      df_out$excl_date_imp <- ifelse(df_out$date_death_imp < date_start_f, T, F)

    # Generate exclusion criterion for decedents who died outside of Sudan
      # variable to exclude because death outside Sudan
      df_out$excl_loc_death <- F
    
      # based on location of death  
      df_out[which(is.na(df_out$loc_death)), "loc_death"] <- "unknown / unclear"
      df_out[which(df_out$loc_death == "outside Sudan"), "excl_loc_death"] <- T
          
      # based on date of departure from Sudan
      df_out$date_dep <- as.Date(paste(df_out$sudan_leave_year, 
        df_out$sudan_leave_month, 15, sep = "-"), "%Y-%m-%d")
      df_out[which(df_out$date_dep < df_out$death_date), "excl_loc_death"] <- T
        
    
#...............................................................................
### ENDS
#...............................................................................
