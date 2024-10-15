#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - SUDAN (2024) ++++++ ###
#...............................................................................

#...............................................................................
## ---- R SCRIPT TO ESTIMATE PROBABILITIES ASSOCIATED WITH EACH SCORE  ------ ##
#...............................................................................


#...............................................................................
### Preparing for estimation
#...............................................................................

  #...................................      
  ## Read and prepare datasets
    
    # Read necessary datasets
    df_base <- readRDS(paste0(dir_path, "out/list_data_base.rds"))    
    ovrlp <- readRDS(paste0(dir_path, "out/ovrlp_base.rds")) 
    df_private <- as.data.frame(
      read_excel(paste0(dir_path, "in/sdn_cam_data_final_private.xlsx")))
    
    # Merge names into dataset and split them into separate columns
    for (i in c("PB", "PV", "SM")) {
      df_private$dup_id <- gsub(paste0(i, "_"), i, df_private$dup_id)
    }
    df_sim <- merge(df_base, df_private[,c("dup_id","cleaned_name","nickname",
      "age", "loc_death_rwd_1", "loc_death_rwd_2", "job", "hcw")],
      by = "dup_id", all.x = T)
    
    # Restrict to Khartoum State
    df_sim <- df_sim[which(df_sim$loc_death == "Khartoum"), ]
    ovrlp_sim <- ovrlp[which(ovrlp$match1_id %in% df_sim$dup_id &
      ovrlp$match2_id %in% df_sim$dup_id), ]
    colnames(ovrlp_sim) <- c("match1_id", "match2_id", "score")
    
    # Omit non-analysable records
    df_sim <- subset(df_sim, del_score == 5)
    
    # Split names into individual variables
    df_sim$cleaned_name <- trimws(gsub("  ", " ", df_sim$cleaned_name))
    x <- strsplit(df_sim$cleaned_name, " ")
    for (i in 1:length(x)) {
      x[[i]] <- c(x[[i]], rep(NA, max(sapply(x, length)) - length(x[[i]])))
    }
    x <- do.call(rbind, x)    
    x <- as.data.frame.matrix(x)
    colnames(x) <- paste0("name_", 1:ncol(x))
    df_sim <- cbind(df_sim, x)
    
    # List of non-missing values for each identifier variable
    df_sim$mmyy <- ((df_sim$year_death - 
      min(df_sim$year_death, na.rm = T)) * 12) + df_sim$month_death
    x <- unlist(df_sim[, paste("name", 1:10, sep = "_")])
    x <- x[! is.na(x)]
    values <- list(
      name_1 = x,
      name_2 = x,
      name_3 = x,
      name_4 = x,
      name_5 = x,
      name_6 = x,
      nickname = na.omit(df_sim$nickname),
      age_cat = na.omit(df_sim$age_cat),
      age = na.omit(df_sim$age),
      mmyy = na.omit(df_sim$mmyy),
      job = na.omit(df_sim$job),
      hcw = df_sim$hcw, # do not omit NA as value = NA matters
      cod = na.omit(df_sim$cod),
      loc_death = na.omit(df_sim$loc_death),
      loc_death_rwd_1 = df_sim$loc_death_rwd_1,
      loc_death_rwd_2 = df_sim$loc_death_rwd_2
    )

    # Dataframe of duplicate pairs
    dup <- df_sim[which(df_sim$dup_score != 0), 
      c("dup_id", "parent_id", "dup_score")]
    colnames(dup) <- c("match1_id", "match2_id", "score")
    dup <- subset(dup, dup$match1_id %in% df_sim$dup_id & 
      dup$match2_id %in% df_sim$dup_id)

    # Initialise output
    n_runs <- 1000
    out_sim <- expand.grid(run = 1:n_runs, par = c("d", "o"), score = 1:5)
    out_sim$prob <- NA
    
    # Progress bar
    pb <- txtProgressBar(min = 1, max = n_runs, style = 3)

    
  #...................................      
  ## Function to compare a pair of records + attribute duplication/overlap score
  f_comp <- function(data_f, df_sim_f = df_sim_i) {
    
    # Recognise pair
    id1 <- df_sim_f[which(df_sim_f$dup_id == data_f[["match1_id"]]), ]
    id2 <- df_sim_f[which(df_sim_f$dup_id == data_f[["match2_id"]]), ]
    
    # Name matches
    n1 <- as.character(id1[paste("name", 1:10, sep = "_")])
    n2 <- as.character(id2[paste("name", 1:10, sep = "_")])
    name_matches <- intersect(n1[! is.na(n1)], n2[! is.na(n2)])
    
    # Age and date matches
    age_matches1 <- id1[["age_cat"]] == id2[["age_cat"]]
    age_matches2 <- (abs(id1[["age"]] - id2[["age"]]) <= 5)
    age_matches <- sum(age_matches1, age_matches2, na.rm = T)
    date_matches <- (abs(id1[["mmyy"]] - id2[["mmyy"]]) <= 6)
    
    # Other variable matches
    other_matches <- c(
      id1[["nickname"]] == id2[["nickname"]],
      id1[["loc_death"]] == id2[["loc_death"]],
      (id1[["loc_death_rwd_2"]] == id2[["loc_death_rwd_2"]]),
      # (id1[["loc_death_rwd_1"]] == id2[["loc_death_rwd_1"]]) |
      #   (id1[["loc_death_rwd_2"]] == id2[["loc_death_rwd_2"]]) |
      #   (id1[["loc_death_rwd_1"]] == id2[["loc_death_rwd_2"]]) |
      #   (id1[["loc_death_rwd_2"]] == id2[["loc_death_rwd_1"]]),
      id1[["cod"]] == id2[["cod"]],
      id1[["job"]] == id2[["job"]],
      id1[["hcw"]] == 1 & id2[["hcw"]] == 1
    )
    
    # Attribute score (result)
    result <- 0
    if (length(name_matches) >= 4) {result <- 5}
    
    if ( (length(name_matches) == 3)) {
      if (age_matches >= 1 & sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 5}
        if (sum(other_matches, na.rm = T) == 2) {result <- 4}
        if (sum(other_matches, na.rm = T) == 1) {result <- 3}
        if (sum(other_matches, na.rm = T) == 0) {result <- 2}
      }
      if (! age_matches >= 1 & sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 4}
        if (sum(other_matches, na.rm = T) == 2) {result <- 3}
        if (sum(other_matches, na.rm = T) == 1) {result <- 2}
        if (sum(other_matches, na.rm = T) == 0) {result <- 1}      
      }
      if (age_matches >= 1 & ! sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 4}
        if (sum(other_matches, na.rm = T) == 2) {result <- 3}
        if (sum(other_matches, na.rm = T) == 1) {result <- 2}
        if (sum(other_matches, na.rm = T) == 0) {result <- 1}      
      }
      if (! age_matches >= 1 & ! sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 3}
        if (sum(other_matches, na.rm = T) == 2) {result <- 2}
        if (sum(other_matches, na.rm = T) == 1) {result <- 1}
        if (sum(other_matches, na.rm = T) == 0) {result <- 0}      
      }
    }

    if ( (length(name_matches) == 2)) {
      if (age_matches >= 1 & sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 4}
        if (sum(other_matches, na.rm = T) == 2) {result <- 3}
        if (sum(other_matches, na.rm = T) == 1) {result <- 2}
        if (sum(other_matches, na.rm = T) == 0) {result <- 1}
      }
      if (! age_matches >= 1 & sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 3}
        if (sum(other_matches, na.rm = T) == 2) {result <- 2}
        if (sum(other_matches, na.rm = T) == 1) {result <- 1}
        if (sum(other_matches, na.rm = T) == 0) {result <- 0}      
      }
      if (age_matches >= 1 & ! sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 3}
        if (sum(other_matches, na.rm = T) == 2) {result <- 2}
        if (sum(other_matches, na.rm = T) == 1) {result <- 1}
        if (sum(other_matches, na.rm = T) == 0) {result <- 0}      
      }
      if (! age_matches >= 1 & ! sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 2}
        if (sum(other_matches, na.rm = T) == 2) {result <- 1}
        if (sum(other_matches, na.rm = T) == 1) {result <- 0}
        if (sum(other_matches, na.rm = T) == 0) {result <- 0}      
      }
    }
      
    if ( (length(name_matches) == 1)) {
      if (age_matches >= 1 & sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 2}
        if (sum(other_matches, na.rm = T) == 2) {result <- 1}
        if (sum(other_matches, na.rm = T) == 1) {result <- 0}
        if (sum(other_matches, na.rm = T) == 0) {result <- 0}
      }
      if (! age_matches >= 1 & sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 1}
        if (sum(other_matches, na.rm = T) == 2) {result <- 0}
        if (sum(other_matches, na.rm = T) == 1) {result <- 0}
        if (sum(other_matches, na.rm = T) == 0) {result <- 0}      
      }
      if (age_matches >= 1 & ! sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 1}
        if (sum(other_matches, na.rm = T) == 2) {result <- 0}
        if (sum(other_matches, na.rm = T) == 1) {result <- 0}
        if (sum(other_matches, na.rm = T) == 0) {result <- 0}      
      }
      if (! age_matches >= 1 & ! sum(date_matches, na.rm = T) == 1) {
        if (sum(other_matches, na.rm = T) >= 3) {result <- 0}
        if (sum(other_matches, na.rm = T) == 2) {result <- 0}
        if (sum(other_matches, na.rm = T) == 1) {result <- 0}
        if (sum(other_matches, na.rm = T) == 0) {result <- 0}      
      }
    }
      
    # Return original score and attributed score
    return(c(data_f[["score"]], result))
  }
        
  
#...............................................................................
### Implementing simulation to estimate probabilities
#...............................................................................

  #...................................      
  ## Run simulation
  for (i in 1:n_runs) {
    
    # Progress
    setTxtProgressBar(pb, i)
    
    # Fresh dataset
    df_sim_i <- df_sim
    
    # Randomly fill in missing variables
    for (j in names(values)) {
      x <- which(is.na(df_sim_i[, j]))
      df_sim_i[x, j] <- sample(values[[j]], length(x), replace = T)
    }
    
    # Compare duplication pairs
    dup_result <- as.data.frame(t(apply(dup, 1, f_comp)))
    colnames(dup_result) <- c("score", "result")
    dup_result$score <- as.numeric(dup_result$score)
    dup_result$result <- as.numeric(dup_result$result)
    dup_result$result <- floor(dup_result$result / 5)
    dup_result <- aggregate(list(result = dup_result$result),
      by = list(score = dup_result$score), FUN = mean)
    dup_result <- dup_result[order(dup_result$score), ]
    
    # Compare overlap pairs
    ovrlp_result <- as.data.frame(t(apply(ovrlp_sim, 1, f_comp)))
    colnames(ovrlp_result) <- c("score", "result")
    ovrlp_result$score <- as.numeric(ovrlp_result$score)
    ovrlp_result$result <- as.numeric(ovrlp_result$result)
    ovrlp_result$result <- floor(ovrlp_result$result / 5)    
    ovrlp_result <- aggregate(list(result = ovrlp_result$result),
      by = list(score = ovrlp_result$score), FUN = mean)
    ovrlp_result <- ovrlp_result[order(ovrlp_result$score), ]
        
    # Record results of run
    out_sim[which(out_sim$run == i & out_sim$par == "d"), c("score", "prob")] <-
      dup_result
    out_sim[which(out_sim$run == i & out_sim$par == "o"), c("score", "prob")] <-
      ovrlp_result
  }  
  close(pb)    

  
  #...................................      
  ## Visualise and store cumulative probability distributions
  
    # Normalise all runs to 1 (not all runs result in 5's being classified
        # as such, as the actual algorithm used relies on more variables
        # or knowledge of individual decedents)
    for (i in 1:n_runs) {
      for (j in unique(out_sim$par))
      x <- which(out_sim$run == i & out_sim$par == j)
      out_sim[x, "prob"] <- out_sim[x, "prob"] / max(out_sim[x, "prob"])
    }
  
    # Visualise distribution of runs, by score and parameter
    df <- subset(out_sim, score != 5)
    df$score <- factor(df$score, levels = 1:4, labels = paste0("score = ", 1:4))
    df$par <- ifelse(df$par == "o", "match probability", 
      "duplication probability")
    ggplot(df, aes(x = prob, colour = score, linetype = par, group = score)) +
      geom_density(linewidth = 1) +
      theme_bw() +
      scale_color_manual("score", values = palette_gen[c(1,5,9,14)]) +
      scale_linetype_manual("par", values = c("11", "solid")) +
      facet_grid(score ~ par) +
      theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(2, "lines")) +
      guides(linetype = "none") +
      scale_x_continuous("estimated probability", limits=c(0,1), expand=c(0,0))
    ggsave(paste0(dir_path, "out/score_probs.png"), units = "cm",
      dpi = "print", width = 20, height = 15)
      
    # Compute cumulative empirical probability distributions
    score_probs <- expand.grid(par = c("d", "o"), score = 0:5,
      p_cum = seq(0, 1, 0.025))
    score_probs$prob <- NA
    for (i in c("d", "o")) {
      out_sim_i <- subset(out_sim, par == i)

      # for each possible score level...
      for (j in 1:4) {

        # select simulation outputs for this score level
        out_sim_ij <- subset(out_sim_i, score == j)

        # take quantiles
        score_probs[which(score_probs$par==i & score_probs$score==j), "prob"] <-
          quantile(out_sim_ij$prob, seq(0, 1, 0.025))
      }
    }    
    score_probs[which(score_probs$score == 0), "prob"] <- 0
    score_probs[which(score_probs$score == 5), "prob"] <- 1
    
    # Save
    saveRDS(score_probs, paste0(dir_path, "out/score_probs.rds"))

            
#...............................................................................
### ENDS
#...............................................................................
