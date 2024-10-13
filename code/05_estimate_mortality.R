#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - SUDAN (2024) ++++++ ###
#...............................................................................

#...............................................................................
## ----- R SCRIPT TO PERFORM CAPTURE-RECAPTURE ESTIMATION OF MORTALITY  ----- ##
#...............................................................................


#...............................................................................
### Estimating deaths by applying random probabilities of duplication / overlap
#...............................................................................

  #...................................      
  ## Read and prepare datasets
    
    # Read necessary datasets
    df_base <- readRDS(paste0(dir_path, "out/list_data_base.rds"))    
    ovrlp <- readRDS(paste0(dir_path, "out/ovrlp_base.rds"))
    score_probs <- readRDS(paste0(dir_path, "out/score_probs.rds"))
    # see_cum <- readRDS(paste0(dir_path, "out/see_cum_dist.rds")) 
    
    # Identify list names
    list_names <- data.frame(list = paste0("list", 1:3), 
      list_name = c("public survey", "private survey", "social media"),
      list_colour = palette_gen[c(3, 8, 13)])
    
    # Start and end dates of study    
    date_start <- as.Date(paste(2023, 4, 15, sep = "-"), "%Y-%m-%d")
    date_end <- as.Date(paste(2024, 6, 4, sep = "-"), "%Y-%m-%d")
    
    # Initialise output
    n_runs <- 10
    out_all_cause <- data.frame(run = 1:n_runs)
    out_all_cause[, c("stratum", "unlisted", "unlisted_lci", "unlisted_uci", 
      "total_deaths_est", "total_deaths_lci", "total_deaths_uci")] <- NA
    out_intl_injury <- out_all_cause
    out_all_cause_sens <- expand.grid(run = 1:n_runs, 
      list = c(list_names$list_name, "all lists"))
    out_all_cause_sens$list <- as.character(out_all_cause_sens$list)
    out_all_cause_sens[, c("n_deaths", "sens_est", "sens_lci", "sens_uci")] <-NA
    out_all_cause_sens <- out_all_cause_sens[order(out_all_cause_sens$run), ]
    out_intl_injury_sens <- out_all_cause_sens
    
    
for (run_i in 1:n_runs) {
  #...................................      
  ## Generate random matched dataset

    # Progress message
    print(paste0(">>> now working on simulation ", run_i, " of ", n_runs))
  
    # Generate random matched dataset based on duplication/overlap probabilities
    df_out <- f_dup(df_f = df_base, random_probs = T)
    df <- f_ovrlp(threshold_ovrlp = NA, random_probs = T)
    # Manage cause of death and add exclusion criterion (not intentional injury)
    df[which(is.na(df$cod)), "cod"] <- "unknown / unclear"
    df$excl_cod <- ifelse(df$cod == "intentional injury", F, T)
    df$cod2 <- ifelse(df$cod == "intentional injury", "intentional injury", 
      "other")
    
    # Manage location of death and add exclusion criterion (not Khartoum State)
    df$loc_death2 <- "other states"
    df[which(df$loc_death == "Khartoum"), "loc_death2"] <- "Khartoum State"
    df[which(df$loc_death == "outside Sudan"), "loc_death2"] <- "outside Sudan"
    df$excl_kht <- ifelse(df$loc_death == "Khartoum", F, T)

    # Manage exclusion date criterion
    df[which(is.na(df$excl_date)), "excl_date"] <- F
    
    # Add/modify missing categories for different variables
    df[which((is.na(df$gender))), "gender"] <- "missing"
    df[which((is.na(df$age_cat))), "age_cat"] <- "missing"
    df[which((is.na(df$resistance_committees))), "resistance_committees"] <- F
    df$resistance_committees <- ifelse(df$resistance_committees, "yes", 
      "no / unknown")
    df[which((is.na(df$year_death))), "year_death"] <- "missing"
    
    # Generate monthly incremental variable as a potential confounder
    df$mmyy <- ((year(df$date_death) - 
      min(year(df$date_death), na.rm = T)) * 12) + month(df$date_death)
    df$mmyy <- df$mmyy - min(df$mmyy, na.rm = T)
    
    # Generate categorical trimester variable as a potential confounder
    df$mmyy_cat <- cut(df$mmyy, c(0,4,8,12,18), include.lowest = T, right = F)

  #...................................      
  ## Estimate mortality
    
    # Select dataset
    df_all <- df[which(!df$excl_del & !df$excl_date & !df$excl_loc_death &
        !df$excl_kht), ]

    # Estimate all-cause deaths and add to output
    out <- f_model(data_f = df_all, confounders = c("mmyy", "n_rep", "cod2"))
    out_all <- f_model_average(confounders = c("mmyy", "n_rep", "cod2"))
    out_all_cause[which(out_all_cause$run == run_i), 
      2:ncol(out_all_cause)] <- 
      out_all$out_est_raw[2:length(out_all$out_est_raw)]
    out_all_cause_sens[which(out_all_cause_sens$run == run_i), 
      2:ncol(out_all_cause_sens)] <- 
      out_all$out_sens_raw[2:ncol(out_all$out_sens_raw)]
    
    # Estimate intentional injury deaths and add to output
    out <- f_model(data_f = df_all[which(!df_all$excl_cod), ],
      confounders = c("mmyy", "n_rep"))
    out_inj <- f_model_average(confounders = c("mmyy", "n_rep"))    
    out_intl_injury[which(out_intl_injury$run == run_i), 
      2:ncol(out_intl_injury)] <- 
      out_inj$out_est_raw[2:length(out_inj$out_est_raw)]
    out_intl_injury_sens[which(out_intl_injury_sens$run == run_i), 
      2:ncol(out_intl_injury_sens)] <- 
      out_inj$out_sens_raw[2:ncol(out_inj$out_sens_raw)]
}  
    
  #...................................      
  ## Average the estimates and save them
        
    # All-cause deaths
    x <- aggregate(out_all_cause[, c("unlisted", "unlisted_lci", "unlisted_uci", 
      "total_deaths_est", "total_deaths_lci", "total_deaths_uci")], 
      by = list(stratum = out_all_cause$stratum), FUN = median, na.rm = T)
    write.csv(x, paste0(dir_path, "out/all_cause_out_est_sim.csv"), row.names=F)
    x <- aggregate(out_all_cause_sens[, c("n_deaths", "sens_est", "sens_lci", 
      "sens_uci")], by = list(list = out_all_cause_sens$list), FUN = median, 
      na.rm = T)
    write.csv(x, paste0(dir_path, "out/all_cause_out_est_sens_sim.csv"), 
      row.names = F)
        
    # Intentional injury deaths
    x <- aggregate(out_intl_injury[, c("unlisted","unlisted_lci","unlisted_uci", 
      "total_deaths_est", "total_deaths_lci", "total_deaths_uci")], 
      by = list(stratum = out_intl_injury$stratum), FUN = median, na.rm = T)
    write.csv(x, paste0(dir_path,"out/intl_injury_out_est_sim.csv"),row.names=F)
    x <- aggregate(out_intl_injury_sens[, c("n_deaths", "sens_est", "sens_lci", 
      "sens_uci")], by = list(list = out_intl_injury_sens$list), FUN = median, 
      na.rm = T)
    write.csv(x, paste0(dir_path, "out/intl_injury_out_est_sens_sim.csv"), 
      row.names = F)
            

#...............................................................................
### Estimating deaths for each combination of duplication and overlap threshold
#...............................................................................

# for each duplication and overlap score combination...
for (dd in 1:5) {
  print(paste0("now estimating mortality for duplication threshold ", dd))

  for (oo in 1:5) {
    print(paste0("   ...and overlap threshold ", oo))

  #...................................      
  ## Read and prepare dataset
    
    # Which directory
    dir_do <- paste0(dir_path, "out/d", dd, "o", oo, "/")
    
    # Read dataset arising from given duplication and overlap thresholds
    df <- readRDS(paste0(dir_do, "list_data.rds"))

    # Identify list names
    list_names <- data.frame(list = paste0("list", 1:3), 
      list_name = c("public survey", "private survey", "social media"),
      list_colour = palette_gen[c(3, 8, 13)])
    
    # Start and end dates of study    
    date_start <- as.Date(paste(2023, 4, 15, sep = "-"), "%Y-%m-%d")
    date_end <- as.Date(paste(2024, 6, 4, sep = "-"), "%Y-%m-%d")

    # Manage cause of death and add exclusion criterion (not intentional injury)
    df[which(is.na(df$cod)), "cod"] <- "unknown / unclear"
    df$excl_cod <- ifelse(df$cod == "intentional injury", F, T)
    df$cod2 <- ifelse(df$cod == "intentional injury", "intentional injury", 
      "other")
    
    # Manage location of death and add exclusion criterion (not Khartoum State)
    df$loc_death2 <- "other states"
    df[which(df$loc_death == "Khartoum"), "loc_death2"] <- "Khartoum State"
    df[which(df$loc_death == "outside Sudan"), "loc_death2"] <- "outside Sudan"
    df$excl_kht <- ifelse(df$loc_death == "Khartoum", F, T)

    # Manage exclusion date criterion
    df[which(is.na(df$excl_date)), "excl_date"] <- F
    
    # Add/modify missing categories for different variables
    df[which((is.na(df$gender))), "gender"] <- "missing"
    df[which((is.na(df$age_cat))), "age_cat"] <- "missing"
    df[which((is.na(df$resistance_committees))), "resistance_committees"] <- F
    df$resistance_committees <- ifelse(df$resistance_committees, "yes", 
      "no / unknown")
    df[which((is.na(df$year_death))), "year_death"] <- "missing"
    

  #...................................      
  ## Estimate mortality
    
    # Select dataset
    df_all <- df[which(!df$excl_del & !df$excl_date & !df$excl_loc_death &
        !df$excl_kht), ]

    # Generate monthly incremental variable as a potential confounder
    df_all$mmyy <- ((year(df_all$date_death) - 
      min(year(df_all$date_death), na.rm = T)) * 12) + month(df_all$date_death)
    df_all$mmyy <- df_all$mmyy - min(df_all$mmyy, na.rm = T)

    # Generate categorical trimester variable as a potential confounder
    df_all$mmyy_cat <- cut(df_all$mmyy, c(0,4,8,12,18), 
      include.lowest = T, right = F)
    
    # Estimate all-cause deaths and save output
    out <- f_model(data_f = df_all,confounders = c("mmyy", "n_rep", "cod2"))
    out_all <- f_model_average(confounders = c("mmyy", "n_rep", "cod2"))
    for (i in names(out_all)) {
      write.csv(out_all[[i]], paste0(dir_do, "all_cause_",i,".csv"),row.names=F)
    }
    
    # Estimate intentional injury deaths and save output
    out <- f_model(data_f = df_all[which(!df_all$excl_cod), ],
      confounders = c("mmyy", "n_rep"))
    out_inj <- f_model_average(confounders = c("mmyy", "n_rep"))    
    for (i in names(out_inj)) {
      write.csv(out_inj[[i]],paste0(dir_do,"intl_injury_",i,".csv"),row.names=F)
    }
  }
}
# closing duplication score and overlap score threshold loops


#...............................................................................
### Combinining estimates for different sensitivity analyses into one table
#...............................................................................

  #...................................      
  ## Construct table
  
    out <- expand.grid(dup_threshold = 1:5, ovrlp_threshold = 1:5)
    out[, paste0("all_cause_", c("conc", "est", "lci", "uci"))] <- NA
    out[, paste0("intl_injury_", c("conc", "est", "lci", "uci"))] <- NA

    # Fill in table for each duplication and overlap score combination...
    for (dd in 1:5) {
      for (oo in 1:5) {
        # concatenated estimates
        x <- read.csv(paste0(dir_path, "out/d", dd, "o", oo, 
          "/all_cause_out_est_pretty.csv"))
        out[which(out$dup_threshold == dd & out$ovrlp_threshold == oo), 
          "all_cause_conc"] <- x$total.deaths..95.CI.
        
        x <- read.csv(paste0(dir_path, "out/d", dd, "o", oo, 
          "/intl_injury_out_est_pretty.csv"))
        out[which(out$dup_threshold == dd & out$ovrlp_threshold == oo), 
          "intl_injury_conc"] <- x$total.deaths..95.CI.

        # point estimates and 95%CI separately 
        x <- read.csv(paste0(dir_path, "out/d", dd, "o", oo, 
          "/all_cause_out_est_raw.csv"))
        out[which(out$dup_threshold == dd & out$ovrlp_threshold == oo), 
          "all_cause_est"] <- x$total_deaths_est
        out[which(out$dup_threshold == dd & out$ovrlp_threshold == oo), 
          "all_cause_lci"] <- x$total_deaths_lci
        out[which(out$dup_threshold == dd & out$ovrlp_threshold == oo), 
          "all_cause_uci"] <- x$total_deaths_uci

        x <- read.csv(paste0(dir_path, "out/d", dd, "o", oo, 
          "/intl_injury_out_est_raw.csv"))
        out[which(out$dup_threshold == dd & out$ovrlp_threshold == oo), 
          "intl_injury_est"] <- x$total_deaths_est
        out[which(out$dup_threshold == dd & out$ovrlp_threshold == oo), 
          "intl_injury_lci"] <- x$total_deaths_lci
        out[which(out$dup_threshold == dd & out$ovrlp_threshold == oo), 
          "intl_injury_uci"] <- x$total_deaths_uci
      }  
    }

    # Save table
    write.csv(out, paste0(dir_path, "out/all_sens_raw.csv"), row.names = F)
    x <-c("dup_threshold","ovrlp_threshold","all_cause_conc","intl_injury_conc")
    write.csv(out[, x],paste0(dir_path,"out/all_sens_pretty.csv"),row.names = F)
    
  #...................................      
  ## Visualise differences among sensitivity estimates
  
    # Prepare for graphing
    out1 <- out[, c("dup_threshold", "ovrlp_threshold", 
      grep("all_cause", colnames(out), value = T))]
    out1$cod = "all cause"
    colnames(out1) <- gsub("all_cause_", "", colnames(out1))
    out2 <- out[, c("dup_threshold", "ovrlp_threshold", 
      grep("intl_injury", colnames(out), value = T))]
    out2$cod = "intentional injury"
    colnames(out2) <- gsub("intl_injury_", "", colnames(out2))
    df <- rbind(out1, out2)
    df$dup_threshold <- paste0("duplication threshold >= ", df$dup_threshold)
    df$ovrlp_threshold <- paste0(">=", df$ovrlp_threshold)
    df$cod <- factor(df$cod, labels = c("all causes", "intentional injury"))
    
    # Graph without 95%CIs
    ggplot(df, aes(x = ovrlp_threshold, y = est, alpha = cod, group = cod,
      colour = dup_threshold, fill = dup_threshold)) +
      geom_point(size = 4) +
      geom_line(linetype = "11", linewidth = 0.75) +
#      geom_errorbar(aes(ymin = lci, ymax = uci)) +
      theme_bw() +
      scale_colour_viridis_d("duplication confidence threshold") +
      scale_fill_viridis_d("match confidence threshold") +
      scale_alpha_manual("cause of death", values = c(0.5, 1)) +
      scale_y_continuous("estimated number of deaths", labels = comma_format(),
        breaks = seq(0, 400000, 20000)) +
      scale_x_discrete("match confidence threshold") +
      facet_grid(cod~dup_threshold, scale = "free_y") +
      theme(legend.position = "none")
    ggsave(paste0(dir_path, "out/all_sens_noci.png"), units = "cm", dpi = "print", 
      height = 15, width = 22)
    
    # Graph without 95%CIs
    ggplot(df, aes(x = ovrlp_threshold, y = est, alpha = cod, group = cod,
      colour = dup_threshold, fill = dup_threshold)) +
      geom_point(size = 4) +
      geom_line(linetype = "11", linewidth = 0.75) +
      geom_errorbar(aes(ymin = lci, ymax = uci)) +
      theme_bw() +
      scale_colour_viridis_d("duplication confidence threshold") +
      scale_fill_viridis_d("match confidence threshold") +
      scale_alpha_manual("cause of death", values = c(0.5, 1)) +
      scale_y_continuous("estimated number of deaths", labels = comma_format(),
        breaks = seq(0, 1000000, 50000)) +
      scale_x_discrete("match confidence threshold") +
      facet_grid(cod~dup_threshold, scale = "free_y") +
      theme(legend.position = "none")
    ggsave(paste0(dir_path, "out/all_sens_ci.png"), units = "cm", dpi = "print", 
      height = 15, width = 22)
    
      
      
#...............................................................................
### ENDS
#...............................................................................
