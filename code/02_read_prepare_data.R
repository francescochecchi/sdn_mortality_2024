#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - SUDAN (2024) ++++++ ###
#...............................................................................

#...............................................................................
## ------- R SCRIPT TO PREPARE ALL SENSITIVITY DATASETS FOR ANALYSIS  ------- ##
#...............................................................................


#...............................................................................
### Reading inputs
#...............................................................................
  
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
      paste0(dir_path, "in/sdn_cam_data_final_public.xlsx"), sheet = "data"))

    # Fix date
    df$death_date <- as.Date(df$death_date)
    
    # Streamline columns
    x <- c("cleaned_name", "name", "nickname", "titles", "loc_death_rwd", "age",
      "loc_death_rwd_1", "loc_death_rwd_2", "month_death_rwd", "cod_rwd", 
      "year_death_rwd", "job", "gender_rwd", "hcw", "cod_avoidable", "source", 
      "date_death_rwd", "comments", "comments_2", "publication_date","source_2")
    df <- df[! colnames(df) %in% x]

#...............................................................................
### Managing different variables in the main dataset
#...............................................................................

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

    # Whether the observation is a parent
    x <- names(which(table(df$parent_id) > 1))
    df$parent <- ifelse(df$parent_id %in% x & df$dup_score == 0, T, F)
    table(df$parent, useNA = "ifany")
    
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
    
    # Recreate number of times the person was mentioned (everyone gets 1,
        # then this gets summed during merging)
    df$n_rep <- 1
    
    # Duplication score
    table(df$dup_score, useNA = "ifany")  

    # Harmonise location of death
    table(df$loc_death, useNA = "ifany")
    df[which(df$loc_death == "Al Qadarif"), "loc_death"] <- "Gedaref"
    df[which(df$loc_death == "West Darfour"), "loc_death"] <- "West Darfur"
    df[which(df$loc_death %in% c("NA", "Egypt borders", "Libya borders", 
      "Sudan")), "loc_death"] <- "unknown / unclear"
    df[which(is.na(df$loc_death)), "loc_death"] <- "unknown / unclear"
    x <- c("Australia", "Chad", "Egypt", "Eritrea",
      "Ethiopia", "France", "Germany", "India", "Kuwait", "Libya", "Mali", 
      "New Zealand", "Qatar", "Saudi Arabia", "South Sudan", "Tanzania", 
      "Türkiye", "UAE", "Uganda", "UK", "USA")
    df[which(df$loc_death %in% x), "loc_death"] <- "outside Sudan"
    table(df$loc_death, useNA = "ifany")
        
    # Cause of death
    table(df$cod, useNA = "ifany")
    df[which(df$cod %in% c("natural death","other","starvation_disease")), 
      "cod"] <- "disease"
    df[which(df$cod == "NA"), "cod"] <- "unknown / unclear"
    df[which(is.na(df$cod)), "cod"] <- "unknown / unclear"
    df$cod <- factor(df$cod, levels = c("disease", "injury_accident", 
      "injury_intentional", "unknown / unclear"), labels = c("disease", 
        "accidental injury", "intentional injury", "unknown / unclear"))
    table(df$cod, useNA = "ifany")

    # Membership of resistance committees   
    table(df$resistance_committees, useNA = "ifany")
    df$resistance_committees <- as.logical(df$resistance_committees)
    table(df$resistance_committees, useNA = "ifany")
        
    # Age of death
    table(df$age_cat, useNA = "ifany")
    
    # Date of death
      # Specify start and end of analysis / data collection period  
      date_start <- as.Date(paste(2023, 4, 15, sep = "-"), "%Y-%m-%d")
      date_end <- as.Date(paste(2024, 6, 4, sep = "-"), "%Y-%m-%d")
      
      # day
      df$day_death <- ifelse(is.na(df$death_date), 15, day(df$death_date))
      table(df$day_death, useNA = "ifany")
    
      # month  
      table(df$month_death, useNA = "ifany")
      
      # year    
      table(df$year_death, useNA = "ifany")
   
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
    sens <- expand.grid(dup_score = 1:5, ovrlp_score = 1:5) 
    sens$sens <- paste0("d", sens$dup_score, "_o", sens$ovrlp_score)
    
    # Variables that need to be merged when there is duplication, and how  
    vars_dup <- data.frame(
      var = c("list", "list_name", "parent_id", "n_rep", "loc_death", "gender", 
        "age_cat", "month_death", "day_death", "year_death", "del_score",
        "resistance_committees", "cod", "sudan_leave_month","sudan_leave_year"),
      how = c("mode", "mode", "mode", "sum", "mode", "mode", 
        "mode", "mean", "mean", "mean", "mean",
        "mode", "mode", "mean", "mean")
    )
      # if how == "mode" and values are discordant, take the most common of 3, 
        # the sole value out of 2/3 or, if 2 values out of 2 are discordant, 
        # set to NA

    # Variables that need to be merged when there is overlap, and how  
    vars_ovrlp <- data.frame(
      var = c(paste0("list",1:3),paste0("list_name",1:3),paste0("final_id",1:3), 
        "n_rep", "loc_death", "gender", "age_cat", "month_death", "day_death",
        "year_death", "del_score", 
        "resistance_committees", "cod", "sudan_leave_month","sudan_leave_year"),
      how = c(rep("which", 3), rep("which", 3), rep("which", 3),   
        "sum", "mode", "mode", "mode", "mean", "mean",
        "mean", "mean",
        "mode", "mode", "mean", "mean")
    )
      # if how == "which", reshape the 1, 2 or 3 values horizontally
    
    
  #...................................      
  ## Identify all possible overlap pairs

    # List all possible overlap pairs
    ovrlp <- c()
    for (i in 1:5) {
      x <- na.omit(df[, c("dup_id", "parent_id", 
        paste0(c("ovrlp_id_", "ovrlp_score_"), i))])
      colnames(x) <- c("dup_id", "parent_id", "ovrlp_id", "ovrlp_score")
      ovrlp <- rbind(ovrlp, x)
    }
      
    # Make sure all possible final IDs after duplicate merging feature in the
        # overlap dataset, for both the observation and the possible match;
        # at this stage we don't know which will be found in the de-duplicated
        # dataset
    colnames(ovrlp) <- c("dup_id", "parent_id", "match_dup_id", "ovrlp_score")
    ovrlp2 <- ovrlp
    ovrlp$final_id <- ovrlp$dup_id
    ovrlp$match_final_id <- ovrlp$match_dup_id
    ovrlp2$final_id <- ovrlp2$parent_id
    ovrlp2$match_final_id <- sapply(strsplit(ovrlp2$match_dup_id, "_"), 
      function(x) {x[1]})
    ovrlp <- rbind(ovrlp, ovrlp2)      
    ovrlp <- ovrlp[, c("final_id", "match_final_id", "ovrlp_score")]
    
    # Eliminate redundancies and mirror images of each pairing
    ovrlp <- unique(ovrlp)
    x <- c("final_id", "match_final_id")
    ovrlp[, x] <- t(apply(ovrlp[, x], 1, sort))
    ovrlp <- unique(ovrlp)
    
      # there are still redundancies due to discordant overlap scores
      # (59, 1 and 1 instances of 2, 3 and 4 redundant pairs, respectively)
      # resolve by adopting the highest overlap score (most conservative)
      ovrlp <- aggregate(list(ovrlp_score = ovrlp$ovrlp_score),
        by = ovrlp[, x], FUN = max)
      
    # Rename columns  
    colnames(ovrlp) <- c("match1_id", "match2_id", "ovrlp_score")
    
        
  #...................................      
  ## Generate all possible sensitivity analysis datasets
    
    # For all possible duplication score thresholds...
    for (i in 1:5) {
      print(paste0("now creating dataset for duplication threshold ", i))
      
      # ... and for all possible overlap score thresholds...
      for (j in 1:5) {
        print(paste0("  ...and overlap threshold ", j))
                
        # generate dataset        
        df_out <- f_dup(threshold_dup = i)
        df_out <- f_ovrlp(threshold_ovrlp = j)
        
        # save to its own directory
        x <- paste0(dir_path, "out/d", i, "o", j)
        suppressWarnings(dir.create(x))
        saveRDS(df_out, paste0(x, "/list_data.rds"))
      }
    }  


#...............................................................................
### ENDS
#...............................................................................
