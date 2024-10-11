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

  #...................................
  ## Set probabilities of duplication and overlap, by score

    # Duplication (within each list)
    dup_probs <- data.frame(dup_score = 0:5,
      dup_prob = c(0, 0.2, 0.4, 0.6, 0.8, 1))

    # Overlap (match across lists)
    ovrlp_probs <- data.frame(ovrlp_score = 0:5,
      ovrlp_prob = c(0, 0.2, 0.4, 0.6, 0.8, 1))
    
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
### Reading in, preparing and visualising expert elicitation data
#...............................................................................

  #...................................      
  ## Read in and and manage data

    # Read data
    see <-data.frame(readxl::read_excel(paste0(dir_path,"in/sdn_see_data.xlsx"), 
      sheet = "see_data"))
    see <- see[order(see$expert, see$pair_id), ]
    
    # Remove missing observations
    x <- grep("value", colnames(see), value = T)
    see <- see[complete.cases(see[, x]), ]
    
    # Rescale the expert answers to proportions
    see[, x] <- see[, x] / 10
    
    # Compute and save experts' scores and weights, and add them to dataset
    experts <- f_see()
    write.csv(experts, paste0(dir_path, "out/see_experts_wts.csv"), row.names=F)
    see <- merge(see, experts[, c("expert", "wt")], by = "expert", all.x = T)
    
    # Streamline dataset
    see$par <- ifelse(is.na(see$dup_score), "o", "d")
    see$score <- ifelse(see$par == "d", see$dup_score, see$ovrlp_score)
    see <- see[, c("expert", "par", "score", x, "wt")]
   
    # If we do not want to take into account experts' weights, set equal weights
    expert_wt <- "yes"
    if (expert_wt == "no") {see$wt <- 1 / nrow(experts)}
    
  #...................................      
  ## Average and visualise expert-elicited probabilities
    
    # Take weighted averages by parameter and score
    x <- by(see, see[, c("par", "score")], function(x) {c(unique(x$par), 
      unique(x$score), weighted.mean(x$value, x$wt))})
    x <- as.data.frame(do.call(rbind, x))
    colnames(x) <- c("par", "score", "prob")
    x$score <- as.numeric(x$score)
    x$prob <- as.numeric(x$prob)
    dup_probs <- data.frame(dup_score = 0:5,
      dup_prob = c(0, x[which(x$par == "d"), "prob"], 1))
    ovrlp_probs <- data.frame(ovrlp_score = 0:5,
      ovrlp_prob = c(0, x[which(x$par == "o"), "prob"], 1))
          
    # Visualise expert answers
    see_agg <- aggregate(list(mean = see$value), 
      by = see[, c("expert", "par", "score", "wt")], FUN = mean)
    see_agg$par <- ifelse(see_agg$par == "o", "match probability", 
      "duplication probability")
    see_agg$score <- factor(see_agg$score, levels = 0:5)
    see_agg <- rbind(see_agg,
      data.frame(expert = "all", par = "match probability",
        score = 1:4, wt = 1, 
        mean=ovrlp_probs[which(ovrlp_probs$ovrlp_score %in% 1:4),"ovrlp_prob"]),
      data.frame(expert = "all", par = "duplication probability",
        score = 1:4, wt = 1, 
        mean = dup_probs[which(dup_probs$dup_score %in% 1:4), "dup_prob"])
    )
    see_agg$expert <- ifelse(see_agg$expert == "all", "all (weighted)", 
      "single expert")
    ggplot(see_agg, aes(x = score, y = mean, colour = score, 
      fill = score, size = wt, shape = expert, linewidth = expert)) +
      geom_point(alpha = 0.7) +
      theme_bw() +
      facet_wrap(. ~ par) +
      theme(legend.position = "top") +
      scale_x_discrete("score") +
      scale_y_continuous("mean expert-elicited probability", labels = percent) +
      scale_color_viridis_d("score") +
      scale_fill_viridis_d("score") +
      scale_shape_manual("expert", values = c(12, 21)) +
      scale_linewidth_manual("expert", values = c(10, 0.5)) +
      scale_size_continuous("expert's weight") +
      guides(colour = "none", fill = "none")
    ggsave(paste0(dir_path, "out/see_distributions.png"), dpi = "print",
      width = 22, height = 12, units = "cm")
    
  # #...................................      
  # ## Compute empirical cumulative distributions for duplication/overlap scores
  # 
  #   # Initialise empirical distribution values for 0.025 intervals from 0 to 1
  #   x_out <- seq(0, 1, 0.025)
  #   see[, paste0("pcum_", x_out)] <- NA
  #   
  #   # Initialise all-expert output
  #   see_all <- expand.grid(par = c("d", "o"), score = 0:5, expert = "all")
  #   see_all[, colnames(see)[! colnames(see) %in% colnames(see_all)]] <- NA
  # 
  # for (i in c("d", "o")) {
  #   
  #   # Select duplication or overlap score answers
  #   see_i <- subset(see, par == i)
  #   
  #   # For each possible score level...
  #   for (j in 0:5) {
  # 
  #     # select answers for this score level
  #     see_ij <- subset(see_i, score == j)
  #           
  #     # linearly interpolate expert quantile distributions
  #     for (k in 1:nrow(see_ij)) {
  #       see_ij[k, paste0("pcum_", x_out)] <- approx(
  #         y = c(0.00, 0.10, 0.50, 0.90, 1.00),
  #         x = c(0, see_ij[k, c("value10", "value50", "value90")], 1),
  #         xout = x_out, ties = "ordered", yright = 2, yleft = 2)$y
  #     }
  #     
  #     # update dataset
  #     see[which(see$par == i & see$score == j), paste0("pcum_",x_out)] <- 
  #       see_ij[, paste0("pcum_",x_out)]
  #     
  #     # compute weighted mean cumulative probability distribution
  #     see_all[which(see_all$par==i & see_all$score==j),paste0("pcum_",x_out)] <- 
  #       apply(see_ij[, paste0("pcum_", x_out)], 2, weighted.mean, 
  #         w = see_ij$wt)
  #   }
  # }  
  #   
  #   # Average expert distributions for each score level
  #   see <- aggregate(see[, paste0("pcum_", x_out)], 
  #     by = see[, c("expert", "par", "score")], mean)
  #   
  #   # Add all-expert means to individual expert distributions
  #   see <- rbind(see, see_all[, colnames(see)])
  #   
  #   # Save
  #   saveRDS(see, paste0(dir_path, "out/see_distributions.rds"))
  # 
  #               
  # #...................................      
  # ## Visualise expert-elicited parameter distributions
  #         
  #   # Compute probability densities
  #   x <- t(apply(see[, grep("pcum_", colnames(see))], 1, diff))
  #   see[, paste0("p_", x_out)] <- cbind(see$pcum_0, x)
  #   
  #   # Reshape long
  #   see_long <- reshape(see, direction = "long", 
  #     varying = grep("p_", colnames(see), value = T),
  #     v.names = "p", timevar = "x", idvar = c("expert", "par", "score"),
  #     times = x_out, drop = grep("pcum_", colnames(see), value = T))
  #   see_long$score <- paste0("score = ", see_long$score)
  #   see_long$par <- factor(see_long$par, levels = c("d", "o"), 
  #     labels = c("potential duplicates (within any list)", 
  #       "potential matches (across lists)"))
  #   
  #   # Plot
  #   ggplot(see_long, aes(x = x, y = p, group = expert, colour = score,
  #     linewidth = expert, linetype = expert)) +
  #     geom_line() +
  #     scale_linewidth_manual(
  #       values = c(1, rep(0.5, unique(length(see_long$expert) - 1))) ) +
  #     scale_linetype_manual(
  #       values = c("solid", rep("21", unique(length(see_long$expert) - 1))) ) +
  #     scale_colour_viridis_d() +
  #     facet_grid(par ~ score) +
  #     theme_bw() +
  #     scale_x_continuous("likelihood that any pair is a duplicate/match") +
  #     scale_y_continuous("expert-elicited probability density") +
  #     theme(legend.position = "none", 
  #       axis.text.x = element_text(angle = 45, hjust = 1))
  #   ggsave(paste0(dir_path, "out/see_distributions.png"),
  #     dpi = "print", units = "cm", width = 18, height = 15)
  # 
  # 
  # #...................................      
  # ## Prepare cumulative expert-elicited distributions for each score level
  #   
  #   # Select only all-expert cumulative distributions
  #   x <- grep("pcum_", colnames(see), value = T)
  #   see_cum <- see[which(see$expert == "all"), c("par", "score", x)] 
  #   
  #   # Reshape long and save
  #   see_cum <- reshape(see_cum, direction = "long", varying = x, 
  #     idvar = c("par", "score"), timevar = "prob", 
  #     times = x_out, v.names ="p_cum")
  #   saveRDS(see_cum, paste0(dir_path, "out/see_cum_dist.rds"))
      
    
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
    
    # Save base datasets
    saveRDS(df, paste0(dir_path, "out/list_data_base.rds"))
    saveRDS(ovrlp, paste0(dir_path, "out/ovrlp_base.rds"))
    
    # For all possible duplication score thresholds...
    for (i in 1:5) {
      print(paste0("now creating dataset for duplication threshold ", i))
      
      # ... and for all possible overlap score thresholds...
      for (j in 1:5) {
        print(paste0("  ...and overlap threshold ", j))
                
        # generate dataset        
        df_out <- f_dup(threshold_dup = i, random_probs = F)
        df_out <- f_ovrlp(threshold_ovrlp = j, random_probs = F)
        
        # save to its own directory
        x <- paste0(dir_path, "out/d", i, "o", j)
        suppressWarnings(dir.create(x))
        saveRDS(df_out, paste0(x, "/list_data.rds"))
      }
    }  


  # #...................................      
  # ## Generate sample of duplication and overlap instances for expert elicitation
  #   
  #   # Eliminate records to be deleted
  #   df <- subset(df, del_score == 5)
  #   
  #   # Number of sets of pairs for each [score_pair1][score_pair2] combination
  #   n_sets <- 5
  #   
  #   # # Set up output
  #   # see_sets <- as.data.frame(rbind(combinations(6,2,0:5), 
  #   #   combinations(6,2,0:5)))
  #   # colnames(see_sets) <- c("score_a", "score_b")
  #   # see_sets$par <- c(rep("dup", nrow(see_sets)/2), 
  #   #   rep("ovrlp", nrow(see_sets)/2))
  #   
  #   # Set up output
  #   see_sets <- expand.grid(par = c("dup", "ovrlp"), set = 1:n_sets,
  #     score_comp = 1:4)
  #   see_sets$score_ref <- 5
  #   see_sets[, c("pair_comp_id1", "pair_comp_id2", "pair_ref_id1", 
  #     "pair_ref_id2")] <- NA
  #       
  #   # Duplication pairs
  #   for (i in 1:4) {
  #     x <- which(df$dup_id != df$parent_id & df$dup_score == i)
  #     x <- sample(x, n_sets, replace = F)
  #     see_sets[which(see_sets$par == "dup" & see_sets$score_comp == i),
  #       c("pair_comp_id1", "pair_comp_id2")] <- 
  #       df[x, c("parent_id", "dup_id")]
  #   }
  #   see_sets <- see_sets[order(see_sets$par,see_sets$set, see_sets$score_comp),]
  #   x <- sort(rep(sample(which(df$dup_id != df$parent_id & df$dup_score == 5), 
  #     n_sets, replace=F), 4))
  #   see_sets[which(see_sets$par=="dup"), c("pair_ref_id1","pair_ref_id2")] <-
  #     df[x, c("parent_id", "dup_id")]
  #    
  #   # Overlap pairs
  #   for (i in 1:4) {
  #     x <- sample(which(ovrlp$ovrlp_score == i), n_sets, replace = F)
  #     see_sets[which(see_sets$par == "ovrlp" & see_sets$score_comp == i),
  #       c("pair_comp_id1", "pair_comp_id2")] <- 
  #       ovrlp[x, c("match1_id", "match2_id")]
  #   }
  #   see_sets <- see_sets[order(see_sets$par,see_sets$set, see_sets$score_comp),]
  #   x <- sort(rep(sample(which(ovrlp$ovrlp_score == 5), n_sets, replace=F), 4))
  #   see_sets[which(see_sets$par=="ovrlp"), c("pair_ref_id1","pair_ref_id2")] <-
  #     ovrlp[x, c("match1_id", "match2_id")]
  #   
  #   # Extract information for each individual in the sets
  #     
  #     # reshape data
  #     x <- reshape(see_sets, direction = "long", varying = c("pair_comp_id1", 
  #       "pair_comp_id2", "pair_ref_id1", "pair_ref_id2"), 
  #       idvar = c("par", "set", "score_comp", "score_ref"), timevar = "type",
  #       times = c("comp_pair_person1", "comp_pair_person2", "ref_pair_person1",
  #         "ref_pair_person2"), v.names = "id")
  #     x$group = ifelse(grepl("comp", x$type), "comparison", "reference")
  #     x$score <- ifelse(x$group == "comparison", x$score_comp, x$score_ref)
  #     x$person <- ifelse(grepl("person1", x$type), "a", "b")  
  #     x$pair_id <- paste(x$par, paste0("set", x$set), x$group, x$score, sep=".")  
  #     x$pair_id <- gsub("reference.5", "reference", x$pair_id)
  #     x <- x[, c("par", "set", "group", "score", "pair_id", "person", "id")]
  #     colnames(x)[ncol(x)] <- "dup_id"
  #     x <- unique(x)
  #     
  #     # read nominal dataset (not uploaded to public Github repository)
  #     private <- data.frame(read_excel(paste0(dir_path, 
  #       "in/sdn_cam_data_final_private.xlsx")))
  #     private$dup_id <- gsub("PB_", "PB", private$dup_id)
  #     private$dup_id <- gsub("PV_", "PV", private$dup_id)
  #     private$dup_id <- gsub("SM_", "SM", private$dup_id)
  #     
  #     # duplicates
  #     x1 <- subset(x, par == "dup")
  #     x1 <- merge(x1, private, by = "dup_id", all.x = T)      
  #   
  #     # matches
  #     x2 <- subset(x, par == "ovrlp")
  #     x2 <- merge(x2, private, by = "dup_id", all.x = T)      
  # 
  #     # all together
  #     x <- rbind(x1, x2)
  #     x <- x[, c("par", "set", "group", "pair_id", "score", "person", 
  #       "dup_id", "gender", "age", "age_cat", "cleaned_name", "name", 
  #       "nickname", "titles", "loc_death", "loc_death_rwd_1", "loc_death_rwd_2", 
  #       "month_death", "year_death", "cod", "cod_rwd", "job", 
  #       "resistance_committees", "hcw", "sudan_leave_year", "sudan_leave_month",
  #       "source", "publication_date")]
  #     write.csv(x, paste0(dir_path, "out/see_sets_data.csv"), row.names = F,
  #       fileEncoding = "UTF-8")
              
#...............................................................................
### ENDS
#...............................................................................
