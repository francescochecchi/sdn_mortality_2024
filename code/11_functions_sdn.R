#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - SUDAN (2024) ++++++ ###
#...............................................................................

#...............................................................................
## ---- R SCRIPT WITH FUNCTIONS USED SPECIFICALLY FOR THE SUDAN ANALYSIS  --- ##
#...............................................................................



#...............................................................................
### Function for first merging duplicates and returning data with 
      # dates reconstructed from merged values
#...............................................................................

f_dup <- function(df_f = df, vars_dup_f = vars_dup, threshold_dup = 3,
  date_start_f = date_start, date_end_f = date_end) {
    
  #...................................      
  ## Prepare for merging
    
    # Initialise a generic one-row dataframe to hold single merged row
    single <- as.data.frame.matrix(matrix(NA, ncol=nrow(vars_dup_f), nrow=1))
    colnames(single) <- vars_dup_f$var

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
    vars_mode <- vars_dup_f[which(vars_dup_f$how == "mode"), "var"]
    for (i in vars_mode) {
      
      # compute mode(s)
      tab <- table(xx[, i])
      mode_x <- suppressWarnings(names(tab[which(tab == max(tab))]))
      
      # if a single most common value, return it, else return NA
      merged[, i] <- ifelse(length(mode_x) == 1, mode_x, NA)
    }

    # Merge variables by mean
    vars_mean <- vars_dup_f[which(vars_dup_f$how == "mean"), "var"]
    merged[, vars_mean] <- colMeans(xx[, vars_mean], na.rm = T)
    
    # Return
    return(merged)
  })
  out_merge <- do.call(rbind, out_merge)

  #...................................      
  ## Return prepared dataset

    # Create final ID variable
    df_nomerge$final_id <- df_nomerge$dup_id
    out_merge$final_id <- out_merge$parent_id
  
    # Update dataset after merging duplicates
    x <- c("final_id", vars_dup_f$var)
    df_out <- rbind(df_nomerge[, x], out_merge[, x])
    
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

#...............................................................................
### Function for first merging overlapping observations and returning data with 
      # average co-variate values and exclusion criteria applied, ready for
      # analysis with generic code
#...............................................................................

f_ovrlp <- function(df_f = df_out, ovrlp_f = ovrlp, vars_ovrlp_f = vars_ovrlp, 
  threshold_ovrlp = 1, date_start_f = date_start, date_end_f = date_end) {   

  #...................................      
  ## Identify a set of unique matches that meet the overlap threshold
  
    # Apply overlap threshold to overlap list
    ovrlp_yes <- subset(ovrlp_f, ovrlp_score >= threshold_ovrlp)
    
    # Only retain pairs for which both people actually feature in the dataset
    ovrlp_yes <- ovrlp_yes[which(ovrlp_yes$match1_id %in% df_f$final_id & 
      ovrlp_yes$match2_id %in% df_f$final_id), ]
    
    # Resolve instances in which a person is paired with >1 others 
        # from the same list by choosing the pairing with highest overlap score;
        # if > 1 matches have the same overlap score, choose the one that isn't
        # a potential duplicate, or else choose the first (very, very few cases)
      # find any such instances and separate them out for fixing
      ovrlp_yes$match1_list <- paste0(ovrlp_yes$match1_id, "x",
        substr(ovrlp_yes$match2_id, 1, 2))
      ovrlp_yes$match2_list <- paste0(ovrlp_yes$match2_id, "x",
        substr(ovrlp_yes$match1_id, 1, 2))
      x <- unique(c(names(which(table(ovrlp_yes$match1_list) > 1)),
        names(which(table(ovrlp_yes$match2_list) > 1))))
      x <- which(ovrlp_yes$match1_list %in% x | ovrlp_yes$match2_list %in% x)
      ovrlp_tofix <- ovrlp_yes[x, ]
      ovrlp_ok <- ovrlp_yes[-x, ]
  
      # group instances
      x <- duplicated(ovrlp_tofix$match1_list)
      x <- unique(ovrlp_tofix[x, "match1_list"])
      x <- data.frame(match1_list = x, group1 = 10000 + 1:length(x))
      ovrlp_tofix <- merge(ovrlp_tofix, x, by = "match1_list", all.x = T)
      x <- duplicated(ovrlp_tofix$match2_list)
      x <- unique(ovrlp_tofix[x, "match2_list"])
      x <- data.frame(match2_list = x, group2 = 20000 + 1:length(x))
      ovrlp_tofix <- merge(ovrlp_tofix, x, by = "match2_list", all.x = T)
      ovrlp_tofix$group <- rowSums(ovrlp_tofix[, c("group1", "group2")],na.rm=T)
  
      # fix instances
      x <- by(ovrlp_tofix, ovrlp_tofix$group, function(xx) {
      
        # is there a single row with the highest overlap score?
        test <- length(which(xx$ovrlp_score == max(xx$ovrlp_score))) == 1
        
        # if yes, adopt that row
        if (test) {return(xx[which.max(xx$ovrlp_score), ])}
        
        # otherwise, adopt row with the shortest ids (= no possible duplicates)
        if (! test) {
          return(xx[which.min(nchar(xx$match1_id) + nchar(xx$match2_id)), ])}
      })
      ovrlp_fixed <- do.call(rbind, x)
      
      # rebind parts together
      ovrlp_yes <- rbind(ovrlp_ok, ovrlp_fixed[, colnames(ovrlp_ok)])
      ovrlp_yes <- ovrlp_yes[, c("match1_id", "match2_id")]

  #...................................      
  ## Identify all overlapping triplets and pairs

    # Identify all overlap triplets
### THIS IS FINE BUT PROBABLY OVERKILL - FEWER LINES OF CODE SHOULD DO IT
    x <- by(ovrlp_yes, ovrlp_yes$match1_id, function(xx) {
      if (nrow(xx) > 1) {return(c(unique(xx$match1_id), xx$match2_id))}
    })
    x1 <- as.data.frame(do.call(rbind, x))
    colnames(x1) <- c("match1_id", "match2_id", "match3_id")
    x <- by(ovrlp_yes, ovrlp_yes$match2_id, function(xx) {
      if (nrow(xx) > 1) {return(c(unique(xx$match2_id), xx$match1_id))}
    })
    x2 <- as.data.frame(do.call(rbind, x))
    colnames(x2) <- c("match1_id", "match2_id", "match3_id")
    x <- which(ovrlp_yes$match1_id %in% ovrlp_yes$match2_id)
    x3 <- ovrlp_yes[x, ] 
    x <- ovrlp_yes[which(ovrlp_yes$match2_id %in% x3$match1_id), ]
    colnames(x) <- c("match3_id", "match1_id")
    x3 <- merge(x3, x, by = "match1_id")
    x <- which(ovrlp_yes$match2_id %in% ovrlp_yes$match1_id)
    x4 <- ovrlp_yes[x, ] 
    x <- ovrlp_yes[which(ovrlp_yes$match1_id %in% x4$match2_id), ]
    colnames(x) <- c("match2_id", "match3_id")
    x4 <- merge(x4, x, by = "match2_id")
    triplets <- rbind(x1, x2, x3, x4)
    triplets <- t(apply(triplets, 1, sort))
    triplets <- as.data.frame(unique(triplets))
    colnames(triplets) <- c("match1_id", "match2_id", "match3_id")
    
    # Identify pairs (not triplets)
    x <- as.vector(unlist(triplets))
    x <- which(ovrlp_yes$match1_id %in% x | ovrlp_yes$match2_id %in% x)
    pairs <- ovrlp_yes[-x, ]
    pairs <- as.data.frame(t(apply(pairs, 1, sort)))
   
  #...................................      
  ## Merge dataset to compose a three-list system

    # Initialise a generic one-row dataframe to hold single merged row
    single <- as.data.frame.matrix(matrix(NA, ncol=nrow(vars_ovrlp_f), nrow=1))
    colnames(single) <- vars_ovrlp_f$var
        
    # Split dataset into observations that are part of triplets, pairs or single
    df_t <- df_f[which(df_f$final_id %in% as.vector(unlist(triplets))), ]
    df_d <- df_f[which(df_f$final_id %in% as.vector(unlist(pairs))), ]
    df_s <- df_f[which(! df_f$final_id %in% df_t$final_id & 
      ! df_f$final_id %in% df_d$final_id), ]    

  #...................................      
  ## Merge triplets and pairs
  for (i in c("_t", "_d")) {
    
    # Select observations
    if (i == "_t") {sets <- triplets} else {sets <- pairs}
    df_i <- get(paste0("df", i))
    
    # Group into sets
    df_i$set <- NA
    for (j in 1:nrow(sets)) {
      df_i[which(df_i$final_id %in% as.character(sets[j, ])), "set"] <- j
    }

    # Remove sets that have less that 2 (for pairs) or 3 (for triplets)
        # (happens in two instances when ovrlp threshold < 2)
    x <- ifelse(i == "_t", 3, 2)
    x <- names(table(df_i$set)[table(df_i$set) < x])
    df_i <- df_i[which(! df_i$set %in% as.integer(x)), ]
 
    # Merge
    out_merge <- by(df_i, df_i$set, function(xx) {
  
      # initialise merged row
      merged <- single
      
      # merge variables by which
      vars_which <- vars_ovrlp_f[which(vars_ovrlp_f$how == "which"), "var"]
      for (j in vars_which) {
        vars_which_j <- gsub("[[:digit:]]", "", j)
        vars_which_row <- as.integer(gsub("[^0-9]", "", j))
        merged[, j] <- xx[vars_which_row, vars_which_j]          
      }
      
      # merge variables by mode
      vars_mode <- vars_ovrlp_f[which(vars_ovrlp_f$how == "mode"), "var"]
      for (j in vars_mode) {
        
        # compute mode(s)
        tab <- table(xx[, j])
        mode_x <- suppressWarnings(names(tab[which(tab == max(tab))]))
        
        # if a single most common value, return it, else return NA
        merged[, j] <- ifelse(length(mode_x) == 1, mode_x, NA)
      }
  
      # merge variables by mean
      vars_mean <- vars_ovrlp_f[which(vars_ovrlp_f$how == "mean"), "var"]
      merged[, vars_mean] <- colMeans(xx[, vars_mean], na.rm = T)
      
      # return
      return(merged)
    })
    
    # Assemble merged dataset
    out_merge <- do.call(rbind, out_merge)
    for (j in c("list", "list_name", "final_id")) {
      x <- paste0(j, 1:3)
      if (i == "_t") {out_merge[,x[1:3]] <- t(apply(out_merge[,x[1:3]],1,sort))}
      if (i == "_d") {out_merge[,x[1:2]] <- t(apply(out_merge[,x[1:2]],1,sort))}
    }
    
    # Assign name
    if (i == "_t") {assign("out_t", out_merge)}
    if (i == "_d") {assign("out_d", out_merge)}
  }    
   
  #...................................      
  ## Assemble the merged dataset
    
    # Harmonise singles dataframe with the other pieces
    x <- c("list", "list_name", "final_id")
    for (i in x) {df_s[, paste0(i, 1)] <- df_s[, i]}
    for (i in x) {df_s[, paste0(i, 2:3)] <- NA}
    
    # Bind all pieces back together
    x <- colnames(out_t)
    df_out <- rbind(df_s[, x], out_t[, x], out_d[, x])
    
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
    
    # Fix list, list name and IDs of matches
    df_out$lists <- paste(df_out$list1, df_out$list2, df_out$list3, sep = "+")
    for (i in paste0("list", 1:3)) {
      df_out[, i] <- as.numeric(grepl(i, df_out$lists))
    }
    df_out$list_names <- paste(df_out$list_name1, df_out$list_name2, 
      df_out$list_name3, sep = "+")
    df_out$ids <- paste(df_out$final_id1, df_out$final_id2, 
      df_out$final_id3, sep = ",")
    x <- c("lists", "list_name1", "list_name2", "list_name3", "final_id1",
      "final_id2", "final_id3")
    df_out <- df_out[, ! colnames(df_out) %in% x]

  #...................................      
  ## Generate exclusion criteria and return prepared dataset
  
    # Exclusion criterion for dates before analysis period
    df_out$excl_date <- ifelse(df_out$date_death < date_start_f, T, F)
    df_out$excl_date_imp <- ifelse(df_out$date_death_imp < date_start_f, T, F)

    # Generate exclusion criterion for decedents who died outside of Sudan
      # variable to exclude because death outside Sudan
      df_out$excl_loc_death <- F
    
      # based on location of death  
      df_out[which(is.na(df_out$loc_death)), "loc_death"] <- "unknown / unclear"
      df_out[which(df_out$loc_death == "outside Sudan"), "excl_loc_death"] <- T
          
      # # based on date of departure from Sudan (not implemented)
      # df_out$date_dep <- as.Date(paste(df_out$sudan_leave_year, 
      #   df_out$sudan_leave_month, 15, sep = "-"), "%Y-%m-%d")
      # df_out[which(df_out$date_dep < df_out$death_date), "excl_loc_death"] <-T

    # Return prepared dataset
    return(df_out)
}        


#...............................................................................
### ENDS
#...............................................................................
