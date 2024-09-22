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
  dup_probs_f = dup_probs, date_start_f = date_start, date_end_f = date_end,
  random_probs = T) {
    
  #...................................      
  ## Prepare for merging
    
    # Initialise a generic one-row dataframe to hold single merged row
    single <- as.data.frame.matrix(matrix(NA, ncol=nrow(vars_dup_f), nrow=1))
    colnames(single) <- vars_dup_f$var

    # Based on merging option chosen, identify observations that need to merge
    if (!random_probs) {
        
      # apply duplicate threshold and identify the observations to merge
      x <- which((df_f$dup_score >= threshold_dup) | df_f$parent)
      df_merge <- df_f[x, ]
    }

    if (random_probs) {
      
      # attribute probabilities of being duplicate
      df_f <- merge(df_f, dup_probs_f, by = "dup_score", all.x = T)
      
      # randomly decide if each possible duplicate will merge, based on probs
      df_f$dup_yes <- (df_f$dup_prob < runif(nrow(df_f)))
      
      # apply duplicate threshold and identify the observations to merge
      x <- which(df_f$dup_yes | df_f$parent)
      df_merge <- df_f[x, ]
    }
        
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
    if (length(vars_mean) > 1) {
      merged[, vars_mean] <- colMeans(xx[, vars_mean], na.rm = T)
    }
    if (length(vars_mean) == 1) {
      merged[, vars_mean] <- mean(xx[vars_mean], na.rm = T)
    }
    
    # Merge variables by sum
    vars_sum <- vars_dup_f[which(vars_dup_f$how == "sum"), "var"]
    if (length(vars_sum) > 1) {
      merged[, vars_sum] <- colsums(xx[, vars_sum], na.rm = T)
    }
    if (length(vars_sum) == 1) {
      merged[, vars_sum] <- sum(xx[vars_sum], na.rm = T)
    }
    
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
      
      # set dates to NA if they are after end of data collection
      df_out[which(df_out$date_death > date_end_f), "date_death"] <- NA

    # Return
    return(df_out)
}


#...............................................................................
### Function to fit each candidate log-linear model for a 3/4-list system
    # as per Rossi et al. https://rivista-statistica.unibo.it/article/view/9854 
############ WORKS, BUT TO BE REWRITTEN
#...............................................................................


f_logl <- function(data_f = df, n_lists = 3, 
  list_names_f = list_names, exposure = NA, confounders = NA) {
  
  #...................................      
  ## Preparatory steps
  
    # Confirm number of lists
      # stop if fewer than three lists or more than four...
      if (! n_lists %in% c(3, 4)  | is.na(n_lists)) 
        {stop("wrong number of lists: only 3 or 4 lists allowed")}
      
      # list names
      lnames <- list_names_f$list_name
      names(lnames) <- paste0("list", 1:3)
  
    # Identify exposure and confounders    
    exposure <- tolower(trimws(exposure))
    if (! is.na(confounders[1])) {confounders <- 
      tolower(sapply(unlist(strsplit(confounders, ",")), trimws)) }
    
    # Prepare dataset
      # restrict data to eligible observations
      if ("eligible" %in% colnames(data_f)) {data_f <- subset(data_f, eligible == T)}
      
      # create unique id for each observation (ignore any existing id variable)
      data_f$key <- paste("id", 1:nrow(data_f), sep = "")
      
      # create indicator variable for outcome (all = 1)
      data_f$y <- 1
      
      # rename list columns
      for (i in 1:n_lists) {
        colnames(data_f)[colnames(data_f) == paste("list", i, sep = "")] <- 
          paste("x", i, sep = "")
      }
      
      # defactor exposure and confounders
      if (! is.na(exposure) ) {
        if (is.factor(data_f[, exposure]) ) {
          data_f[, exposure] <- as.character(data_f[, exposure])
        }
      }  
      if (! is.na(confounders[1]) ) {
        for (i in confounders) {
          if (is.factor(data_f[, i]) ) {
            data_f[, i] <- as.character(data_f[, i])
          }
        }
      }
    
  #...................................      
  ## Prepare data for modelling (one row = one observation)
  
    # Prepare possible 'profiles' for each individual observation
      # identify all possible list profiles (cells in contingency table)
      profiles <- as.data.frame(permutations(n = 2, r = n_lists, v = c(0, 1), 
        repeats.allowed = T))
      colnames(profiles) <- as.character(1:n_lists)
      
      # columns for which observations appear on which combinations of lists
      for (i in 3:n_lists) {
        x1 <- combinations(n = n_lists, r = i - 1, v = 1:n_lists, 
          repeats.allowed = F)
        for (j in 1:nrow(x1)) { profiles[, paste(x1[j, ], collapse = "")] <- 
          rowSums(profiles[, x1[j, ]]) }
        profiles[, nchar(colnames(profiles)) == (i - 1)] <- 
          ifelse(profiles[, nchar(colnames(profiles)) == (i-1)] == (i-1), 1, 0)
      }
      colnames(profiles) <- paste("x", colnames(profiles), sep = "")
      x2 <- colnames(profiles)
    
    # Set profiles for each observation, based on the data
      # create expanded dataframe with all possible profiles for each obs.
      df <- expand_grid(data_f[, "key"], profiles)
      colnames(df) <- c("key", colnames(profiles))
      df <- as.data.frame(df)
      
      # match each observation to possible profiles 
        # (except for which lists the observation is in)
      x1 <- paste("x", 1:n_lists, sep = "")
      df <- merge(df, data_f[, ! colnames(data_f) %in% c(x1, "y")], 
        by = "key", all.x = T)
      
      # determine which profile the observation has, based on which lists 
         # the observation is in
      # outcome indicator = 1 if observation falls within a given profile, 
         # 0 otherwise and NA for x000(0) profile)
      df <- merge(df, data_f[, c("key", x1, "y")], by = c("key", x1), all.x = T)
      df[, "y"] <- ifelse(is.na(df[, "y"]), 0, df[, "y"])
      df[rowSums(df[, x1]) == 0, "y"] <- NA
    
    # Add columns for interactions between exposure and confounder (if present) 
        # and lists
        # (note: omit interactions among confounders and between 
        # exposure and confounders)
      # profiles:exposure interactions
      if (! is.na(exposure)) { 
        df[, paste(x2, exposure, sep = ":")] <- 
          ifelse(df[, x2] == 1, df[, exposure], 0)
      }
      
      # profiles:confounder(s) interactions
      if (! is.na(confounders[1])) {
        for (i in confounders) {df[, paste(x2, i, sep = ":")] <- 
          ifelse(df[, x2] == 1, df[, i], 0)} 
      }
    
  
  #...................................      
  ## Define candidate models
  
    # Possible terms...
    x1 <- c()
    for (i in 3:n_lists) {
      x1 <- c(x1, apply( combinations(n_lists, i - 1, unlist(lnames) ) , 1 , 
        paste , collapse = " x " ) )
    }
    # ...and their length
    x2 <- lapply(sapply(x1, strsplit, " x "), length)
    
    # Combinations of two-list terms
    x3 <- names(x2[x2 == 2])
    x4 <- c()
    for (i in n_lists:1 ) {
      x4 <- c(x4, apply(combinations(length(x3), i, x3) , 1 , paste , 
        collapse = ", " ) )
    }
      # if 3 lists, stop here and attribute to model
      if (n_lists == 3) {out <- data.frame("model" = c("no interactions", x4) )}
  
    # If 4 lists, continue:
    if (n_lists == 4) {
      # combinations of three-list terms
      x3 <- names(x2[x2 == 3])
      x5 <- c()
      if (length(x3) > 0) {
        for (i in n_lists:1 ) {
          x5 <- c(x5, apply(combinations(length(x3), i, x3) , 1 , paste , 
            collapse = ", " ) )
        }
      }
      
      # matrix to indicate overlap of two- and three-list term combinations
        # (overlap = violation of hierarchy principle)
      x7 <- rep(NA, length(x5))
      for (i in 1:length(x5) ) {
        x8 <- unlist(strsplit(x5[i], ", "))
        x9 <- c()
        for (j in 1:length(x8)) {
          x10 <- unlist(strsplit(x8[j], " x ") )
          x9 <- c(x9, apply(permutations(length(x10), 3, x10), 1, paste, 
            collapse = " x ") )
        }
        x7[i] <- paste(x9, collapse = ", ")
      }
      
      x6 <- as.data.frame(matrix(sapply(gsub(", ", "|", x4), grepl, x7), 
        nrow = length(x5), ncol = length(x4)))
      colnames(x6) <- x4
      rownames(x6) <- x5
      
      # select unique non-overlapping combinations
      x8 <- c()
      for (i in 1:nrow(x6)) {
        for (j in 1:ncol(x6)) {
          if (x6[i, j] == T) {x8 <- c(x8, rownames(x6)[i], colnames(x6)[j])}
          if (x6[i, j] == F) {x8 <- c(x8, rownames(x6)[i], colnames(x6)[j], 
            paste(rownames(x6)[i], colnames(x6)[j], sep = ", "))}
        }
      }
      x8 <- unique(x8)
      
      # attribute to output
      out <- data.frame("model" = c("no interactions", x8) )
  }
 
    # Derive model formulae from the above
    x1 <- out[, "model"]
    for (i in 1:length(lnames) ) {
      x1 <- gsub(lnames[i], gsub("list", "", names(lnames[i]) ), x1)
    }
    x1 <- suppressWarnings(do.call(rbind, strsplit(x1, ", ")))
      # this will result in elements being recycled because of the 
      # uneven n of resulting columns, which the next few lines deal with
    x2 <- apply(x1, c(1, 2), function(x) {
      x3 <- strsplit(x, " x "); x3 <- sort(unlist(x3))
      return(paste(x3, collapse = ""))
    })
    x3 <- apply(x2, 1, function(x) {unique(paste("x", x, sep = ""))})
    x3 <- sapply(x3, paste, collapse = " + ")
    x1 <- paste(paste("x", 1:n_lists, sep = ""), collapse = " + ")
    x2 <- paste("y ~  ", x1, " + ", x3, sep = "")
    out[, "formula"] <- x2
    out[1, "formula"] <- paste("y ~  ", x1, sep = "")
    
    # Lastly, add terms for exposure, confounders and interactions of these 
        # with profiles
      # profiles:exposure interactions
      if (! is.na(exposure)) { 
        for (i in 1:nrow(out)) {  
          x1 <- all.vars(as.formula(out[i, "formula"]))[-1]
          x1 <- paste(x1, exposure, sep = ":")
          x1 <- paste(x1, collapse = " + ")
          out[i, "formula"] <- paste(out[i, "formula"], exposure, x1,sep =" + ")
        }
      }  
      
      # profiles:confounder(s) interactions
      if (! is.na(confounders[1])) {
        for (i in 1:nrow(out)) {  
          for (j in confounders) {
            x1 <- all.vars(as.formula(out[i, "formula"]))[-1]
            if (length(grep(paste(confounders, collapse = "|"), x1)) > 0 ) {
              x1 <- x1[-grep(paste(confounders, collapse = "|"), x1)]
            }
            if (! is.na(exposure)) {x1 <- x1[-grep(exposure, x1)]}
            x1 <- paste(x1, j, sep = ":")
            x1 <- paste(x1, collapse = " + ")
            out[i, "formula"] <- paste(out[i, "formula"], j, x1, sep = " + ")
          }
        }
      }
    
  
  #...................................      
  ## Fit all possible log-linear models and calculate statistics
  
    # Statistics to compute for each model
      # symbol for unlisted observations
      if (n_lists == 3) {unlisted <- "m000"}
      if (n_lists == 4) {unlisted <- "m0000"}
      
      # unlisted observations
      out[, paste(unlisted, "screen", sep = "_")] <- NA
      
      # unlisted observations by exposure category, if exposure is categorical
      if (! is.na(exposure)) {
        if (length(levels(as.factor(df[, exposure]))) < 10 ) {
          out[, paste(unlisted, "screen", 
            levels(as.factor(df[, exposure])), sep = "_")] <- NA
        }
      }
      
      # unlisted observations for model without each of the confounders, 
          # or no adjustment for confounders
      if (! is.na(confounders[1])) {
        out[, paste(unlisted, "without", c(confounders, "adjustment"), 
          sep = "_")] <- NA
      }
      
      # other statistics  
      out[, "lrt_p"] <- NA
      
    # Fit saturated model (needed to perform likelihood-ratio test 
          # for models nested within it)
      # identify saturated model
      x1 <- lapply(out[, "model"], function(x) {unlist(strsplit(x, ","))} )
      x2 <- which.max(unlist(lapply(x1, length)))
      
      # fit saturated model
      fit_sat <- try(glm(formula = as.formula(out[x2, "formula"]), 
        family = "poisson", data = df, maxit = 1000 ), silent = T)
    
    # Fit all other models        
    for (i in 1:nrow(out) ) {
      print(paste0("now fitting candidate model ", i, " of ", nrow(out)))
      
      # fit model (or at least try)
      suppressWarnings(rm(fit) )
      fit <- try(glm(formula = as.formula(out[i, "formula"]), family="poisson", 
        data = df, maxit = 1000 ), silent = T)
      
      # check: if model has not fit or any coefficients is NA, skip to next loop
      if (class("fit")[1] == "try-error" ) {next}
      if (any(is.na(coef(fit))) ) {next}
      
      # estimate expected number unlisted
      # select data with 000[0] profile
      df0 <- subset(df, is.na(y))
      
      # probability of being unlisted for each observation
      df0[, paste(unlisted, "screen", sep = "_")] <- predict(fit, df0)
      
      # number unlisted overall and by exposure level
      out[i, paste(unlisted, "screen", sep = "_")] <- 
        round(sum(exp(na.omit(df0[, paste(unlisted, "screen", sep = "_")]))), 
          digits = 0)
      if (! is.na(exposure)) {
        if( length(levels(as.factor(df[, exposure]))) < 10 ) {
          out[i, paste(unlisted, "screen", levels(as.factor(df[, exposure])), sep = "_")] <- 
            aggregate(df0[, paste(unlisted, "screen", sep = "_")],
              by = list(as.factor(df0[, exposure])), 
              FUN = function(x) {round(sum(exp(na.omit(x))), digits = 0)})[, 2]
        }
      }  
      
      # number unlisted if each confounder is taken out, or all are taken out
      if (! is.na(confounders[1])) {
        
        # formula terms
        x1 <- gsub(" + ", " , ", out[i, "formula"], fixed = T)
        x1 <- gsub("y ~ ", "", x1)
        x1 <- sapply(unlist(strsplit(x1, " , ")), trimws)
        
        # take out each of the confounders
        for (j in confounders) {
          # rewrite formula without confounder
          x2 <- x1[-grep(j, x1)]
          x2 <- as.formula(paste("y ~ ", paste(x2, collapse = " + "), sep = ""))
          
          # update fit with new formula
          x3 <- try(update(fit, formula = x2), silent = T)
          
          # probability of being unlisted for each observation
          x4 <- predict(x3, df0)
          
          # number unlisted overall and by exposure level
          out[i, paste(unlisted, "without", j, sep = "_")] <- 
            round(sum(exp(na.omit(x4))), digits = 0)
        }
        
        # take out all confounders
          # rewrite formula without confounders
          x2 <- x1[-grep(paste(confounders, collapse = "|"), x1)]
          x2 <- as.formula(paste("y ~ ", paste(x2, collapse = " + "), sep = ""))
          
          # update fit with new formula
          x3 <- try(update(fit, formula = x2), silent = T)
        
        # probability of being unlisted for each observation
        x4 <- predict(x3, df0)
        
        # number unlisted overall and by exposure level
        out[i, paste(unlisted, "without_adjustment", sep = "_")] <- 
          round(sum(exp(na.omit(x4))), digits = 0)
      }  
      
    # Likelihood ratio test p-value comparing model to saturated model 
        # (low p = model is better than saturated model)
    x1 <- -2 * (logLik(fit) - logLik(fit_sat))
    out[i, "lrt_p"] <- as.numeric(pchisq(x1, 
      df = fit$df.residual - fit_sat$df.residual, lower.tail = F))
  }
  
  #...................................      
  ## Return output and dataframe used for model fitting
  x1 <- list(out, df)
  names(x1) <- c("out", "df")
  return(x1)    
}      


#...............................................................................
### Function to perform model averaging for all eligible candidate models
############ WORKS, BUT TO BE REWRITTEN
#...............................................................................

f_model_average <- function(f_out = out, n_lists = 3, list_names_f = list_names,
  exposure = NA, confounders = NA, plausibility_1 = 100, 
  plausibility_2 = 0.6) {
  
  #...................................      
  ## Preparatory steps
  
    # Confirm number of lists
      # stop if fewer than three lists or more than four...
      if (! n_lists %in% c(3, 4)  | is.na(n_lists)) 
        {stop("wrong number of lists: only 3 or 4 lists allowed")}
    
      # list names
      lnames <- list_names_f$list_name
      names(lnames) <- paste0("list", 1:3)
      
    # Identify exposure  
    exposure <- tolower(exposure)
    
    # Symbol for unlisted observations
    if (n_lists == 3) {unlisted <- "m000"}
    if (n_lists == 4) {unlisted <- "m0000"}
    
    # Dataframe for model fitting
    df <- f_out[["df"]]
    
      # select data with 000[0] profile
      df0 <- subset(df, is.na(y))
      
    # Set up output of model averaging
    out <- f_out[["out"]]
    
      # model eligibility and statistics
      out[, c("eligible", "aic", "post_prob")] <- NA
      
      # unlisted observations
      out[, unlisted] <- NA
      out[, paste(unlisted, "lci", sep = "_")] <- NA
      out[, paste(unlisted, "uci", sep = "_")] <- NA
      
      # unlisted observations by exposure category, if exposure is categorical
      if (! is.na(exposure)) {
        if (length(levels(as.factor(df[, exposure]))) < 10 ) {
          for (j in levels(as.factor(df[, exposure]))) {
            out[, paste(unlisted, j, sep = "_")] <- NA
            out[, paste(unlisted, j, "lci", sep = "_")] <- NA
            out[, paste(unlisted, j, "uci", sep = "_")] <- NA
          }
        }  
      }
      
    
  #...................................      
  ## Select models for averaging and assign explanations for exclusion    
  out[, "eligible"] <- "yes"
  for (i in 1:nrow(out)) {
    if (is.na(out[i, paste(unlisted, "screen", sep = "_")]) )
    {out[i, "eligible"] <- "no - model error or sparse dataset"}
    if (out[i, paste(unlisted, "screen", sep = "_")] / 
      sum(df[, "y"], na.rm = T) > plausibility_1) 
    {out[i, "eligible"] <- "no - implausible estimate"}
    if (out[i, "lrt_p"] > plausibility_2) 
    {out[i, "eligible"] <- "no - possible over-fitting"}
  }
  
  #...................................      
  ## Re-fit each eligible model, compute AIC and predict n unlisted
  for (i in 1:nrow(out)) {
    if (out[i, "eligible"] == "yes") {
      
      # Refit model      
      fit <- try(glm(formula = as.formula(out[i, "formula"]),family = "poisson", 
        data = df, maxit = 1000 ), silent = T)
      
      # Compute model's AIC
      out[i, "aic"] <- round(AIC(fit), digits = 2)
      
      # Predict unlisted observations
        # contribution to total unlisted for each observation - 
          # point estimate and 95%CI
        x1 <- predict(fit, df0, se.fit = T)
        df0[, unlisted] <- exp(x1[[1]])
        df0[, paste(unlisted, "lci", sep = "_")] <- exp(x1[[1]] - 1.96*x1[[2]])
        df0[, paste(unlisted, "uci", sep = "_")] <- exp(x1[[1]] + 1.96*x1[[2]])
        
        # estimated number unlisted overall... 
        out[i, unlisted] <- round(sum(na.omit(df0[, unlisted])), digits = 0)
        out[i, paste(unlisted, "lci", sep = "_")] <- 
          round(sum(na.omit(df0[, paste(unlisted, "lci", sep = "_")])),digits=0)
        out[i, paste(unlisted, "uci", sep = "_")] <- 
          round(sum(na.omit(df0[, paste(unlisted, "uci", sep = "_")])), digits=0)
      
      # ...and by exposure level
      if (! is.na(exposure)) { 
        if (length(levels(as.factor(df[, exposure]))) < 10 ) {
          out[i, paste0(unlisted, levels(as.factor(df[, exposure])))] <- 
            aggregate(df0[, unlisted], by = list(as.factor(df0[, exposure])), 
              FUN = function(x) {round(sum(na.omit(x)), digits = 0)})[, 2]
          out[i, paste0(unlisted, levels(as.factor(df[, exposure])), "lci")] <- 
            aggregate(df0[, paste(unlisted, "lci", sep = "_")],
              by = list(as.factor(df0[, exposure])), 
              FUN = function(x) {round(sum(na.omit(x)), digits = 0)})[, 2]
          out[i, paste0(unlisted, levels(as.factor(df[, exposure])), "uci")] <- 
            aggregate(df0[, paste(unlisted, "uci", sep = "_")], 
              by = list(as.factor(df0[, exposure])), 
              FUN = function(x) {round(sum(na.omit(x)), digits = 0)})[, 2]
        }
      }
    }
  }
  
  #...................................      
  ## Calculate 'posterior probabilities' (weights) of eligible models 
      # from their AICs (based on Rossi, 2010: 
      # https://rivista-statistica.unibo.it/article/view/3593/2945 )
  
  # Calculate an AIC delta based on lowest one
  out[, "aic_delta"] <- NA
  out[out[, "eligible"] == "yes", "aic_delta"] <- 
    out[out[, "eligible"] == "yes", "aic"] - 
    min(out[out[, "eligible"] == "yes", "aic"])
  
  # Then calculate weights / posterior probabilities 
  tot_prob <- sum(exp(- out[out[, "eligible"] == "yes", "aic_delta"] / 2))
  out[out[, "eligible"] == "yes", "post_prob"] <- 
    exp(- out[out[, "eligible"] == "yes", "aic_delta"] / 2) / tot_prob
  
  #...................................      
  ## Compute and return overall results
  
  # Format model output thus far
  out_raw <- out
  out_raw <- cbind(rep("", nrow(out_raw)), out_raw)
  colnames(out_raw)[1] <- "-"
  
  out_pretty <- out
  out_pretty[, "deaths outside any list (95%CI)"] <- paste(
    out_pretty[, unlisted], " (", 
    out_pretty[, paste(unlisted, "lci", sep = "_")], " to ", 
    out_pretty[, paste(unlisted, "uci", sep = "_")], ")", sep = "")
  if (! is.na(exposure)) { 
    if (length(levels(as.factor(df[, exposure]))) < 10 ) {
      for (i in levels(as.factor(df[, exposure])) ) {
        out_pretty[, paste("deaths outside any list (95%CI)", i, sep = " - ")]<-
          paste(
          out_pretty[, paste(unlisted, i, sep = "_")], " (", 
          out_pretty[, paste(unlisted, i, "lci", sep = "_")], " to ", 
          out_pretty[, paste(unlisted, i, "uci", sep = "_")], ")", sep = "")
      }  
    }
  }
  out_pretty[, "likelihood ratio p-value"] <- out_pretty[, "lrt_p"]
  out_pretty[, "AIC"] <- round(out_pretty[, "aic"], 2)
  out_pretty[, "posterior probability"] <- round(out_pretty[, "post_prob"], 3)
  out_pretty <- out_pretty[, c("model", grep("outside", colnames(out_pretty), 
    value = T), "likelihood ratio p-value", "AIC", "posterior probability")]
  out_pretty <- cbind(rep("", nrow(out_pretty)), out_pretty)
  colnames(out_pretty)[1] <- "-"
  
  # Estimated unlisted and total deaths with 95%CIs
  out_est_raw <-c(
    weighted.mean(out[, unlisted], out$post_prob, na.rm = T),
    weighted.mean(out[, paste(unlisted, "lci", sep = "_")], 
      out$post_prob, na.rm = T),
    weighted.mean(out[, paste(unlisted, "uci", sep = "_")], 
      out$post_prob, na.rm = T)
  )
  out_est_raw <-c(out_est_raw, sum(df[, "y"], na.rm = T) + out_est_raw)
  if (! is.na(exposure)) { 
    if (length(levels(as.factor(df[, exposure]))) < 10 ) {
      for (i in levels(as.factor(df[, exposure])) ) {
        x1 <- c(weighted.mean(out[, paste(unlisted, i, sep = "_")], 
          out$post_prob, na.rm = T),
          weighted.mean(out[, paste(unlisted, i, "lci", sep = "_")], 
            out$post_prob, na.rm = T),
          weighted.mean(out[, paste(unlisted, i, "uci", sep = "_")], 
            out$post_prob, na.rm = T)
        )
        out_est_raw <- 
          rbind(out_est_raw, 
            c(x1, sum(df[which(df[, exposure] == i), "y"], na.rm = T) + x1)
        )
      }  
    }
  }
  
  out_est_raw <- as.data.frame(rbind(out_est_raw))
  out_est_raw <- round(out_est_raw, 0)
  out_est_raw <- cbind(rep(NA, nrow(out_est_raw)), out_est_raw)
  colnames(out_est_raw) <- c("stratum", "unlisted", "unlisted_lci", 
    "unlisted_uci", "total_deaths_est", "total_deaths_lci", "total_deaths_uci")
  out_est_raw[1, "stratum"] <- "overall"
  if (! is.na(exposure)) { 
    if (length(levels(as.factor(df[, exposure]))) < 10 ) {
      for (i in 1:length(levels(as.factor(df[, exposure]))) ) {
        out_est_raw[1+i, "stratum"] <- levels(as.factor(df[, exposure]))[i]
      }  
    }
  }
  out_est_raw <- cbind(rep("", nrow(out_est_raw)), out_est_raw)
  colnames(out_est_raw)[1] <- "-"
  
  out_est_pretty <- as.data.frame(matrix(NA, nrow = nrow(out_est_raw),ncol = 3))
  for (i in 1:nrow(out_est_raw)) {
    out_est_pretty[i, 1] <- out_est_raw[i, 2]
    out_est_pretty[i, 2] <- paste(out_est_raw[i, 3], " (", out_est_raw[i, 4], 
      " to ", out_est_raw[i, 5], ")", sep = "")
    out_est_pretty[i, 3] <- paste(out_est_raw[i, 6], " (", out_est_raw[i, 7], 
      " to ", out_est_raw[i, 8], ")", sep = "") 
  }
  colnames(out_est_pretty) <- c("stratum", "deaths outside any list (95%CI)", 
    "total deaths (95%CI)")
  out_est_pretty <- cbind(rep("", nrow(out_est_pretty)), out_est_pretty)
  colnames(out_est_pretty)[1] <- "-"
  
  # Sensitivity of each list and all lists combined
  out_sens_raw <- as.data.frame(matrix(NA, ncol = 5, nrow = n_lists + 1))
  colnames(out_sens_raw) <- c("list", "n_deaths", "sens_est", "sens_lci", 
    "sens_uci")
  out_sens_raw[, "list"] <- c(unlist(lnames), "all lists")
  
  for (i in 1:n_lists) {
    out_sens_raw[i, "n_deaths"] <- 
      sum(df[df[, grep(as.character(i), colnames(df))] == 1, "y"], na.rm = T)
  }
  
  # Calculate sensitivity for all lists combined
  out_sens_raw[n_lists + 1, "n_deaths"] <- sum(out_sens_raw[1:n_lists, 
    "n_deaths"], na.rm = T)
  
  # Calculate sensitivity estimates, LCI, and UCI
  out_sens_raw[, "sens_est"] <- out_sens_raw[, "n_deaths"] / 
    out_est_raw[1, "total_deaths_est"]
  out_sens_raw[, "sens_lci"] <- out_sens_raw[, "n_deaths"] / 
    out_est_raw[1, "total_deaths_uci"]
  out_sens_raw[, "sens_uci"] <- out_sens_raw[, "n_deaths"] / 
    out_est_raw[1, "total_deaths_lci"]
  
  out_sens_raw <- cbind(rep("", nrow(out_sens_raw)), out_sens_raw)
  colnames(out_sens_raw)[1] <- "-"
  
  out_sens_pretty <- out_sens_raw
  out_sens_pretty[, "number of deaths"] <- out_sens_pretty[, "n_deaths"]
  out_sens_pretty[, "sensitivity (95%CI)"] <- 
    paste(round(out_sens_pretty[, "sens_est"] * 100, 1), "% (",
      round(out_sens_pretty[, "sens_lci"] * 100, 1), "% to ", 
      round(out_sens_pretty[, "sens_uci"] * 100, 1), "%)", sep = "")
  out_sens_pretty <- 
    out_sens_pretty[, c("list", "number of deaths", "sensitivity (95%CI)")]
  out_sens_pretty <- cbind(rep("", nrow(out_sens_pretty)), out_sens_pretty)
  colnames(out_sens_pretty)[1] <- "-"
  
  
  # Return outputs as a list
  x1 <- list("out_raw" = out_raw, "out_pretty" = out_pretty, 
    "out_est_raw" = out_est_raw, "out_est_pretty" = out_est_pretty, 
    "out_sens_raw" = out_sens_raw, "out_sens_pretty" = out_sens_pretty)
  return(x1)
}



#...............................................................................
### Function for first merging overlapping observations and returning data with 
      # average co-variate values and exclusion criteria applied, ready for
      # analysis with generic code
#...............................................................................

f_ovrlp <- function(df_f = df_out, ovrlp_f = ovrlp, vars_ovrlp_f = vars_ovrlp, 
  threshold_ovrlp = 3, ovrlp_probs_f = ovrlp_probs, date_start_f = date_start, 
  date_end_f = date_end, random_probs = T) {   

  #...................................      
  ## Identify a set of unique matches that overlap, based on threshold or
      # random probabilities

    # Based on option chosen, identify observations that need to merge
    if (!random_probs) {

      # apply overlap threshold to overlap list
      ovrlp_yes <- subset(ovrlp_f, ovrlp_score >= threshold_ovrlp)
    }

    if (random_probs) {
      
      # attribute probabilities of being a match
      ovrlp_f <- merge(ovrlp_f, ovrlp_probs_f, by = "ovrlp_score", all.x = T)
      
      # randomly decide who is a match, based on probs
      ovrlp_f$ovrlp_yes <- (ovrlp_f$ovrlp_prob < runif(nrow(ovrlp_f)))
      ovrlp_yes <- ovrlp_f[which(ovrlp_f$ovrlp_yes), ]
    }
  
    # Only retain pairs for which both people actually feature in the dataset
    ovrlp_yes <- ovrlp_yes[which(ovrlp_yes$match1_id %in% df_f$final_id & 
      ovrlp_yes$match2_id %in% df_f$final_id), ]
    
    # Resolve instances in which a person is paired with >1 others 
        # from the same list by choosing the pairing with highest overlap score;
        # if > 1 matches have the same overlap score, choose the one that isn't
        # a potential duplicate, or else choose the first (very, very few cases)
      
      # first find any such instances in col 1 and separate them out for fixing
      ovrlp_yes$match1_list <- paste0(ovrlp_yes$match1_id, "x",
        substr(ovrlp_yes$match2_id, 1, 2))
      x <- unique(names(which(table(ovrlp_yes$match1_list) > 1)))
      x <- which(ovrlp_yes$match1_list %in% x)
      
        # if there are observations to fix...
        if (length(x) > 0) {
          
          # split dataset into part to fix and part that is ok
          ovrlp_tofix <- ovrlp_yes[x, ]
          ovrlp_ok <- ovrlp_yes[-x, ]
  
          # group instances
          x <- duplicated(ovrlp_tofix$match1_list)
          x <- unique(ovrlp_tofix[x, "match1_list"])
          if (length(x) != 0) {
            x <- data.frame(match1_list = x, group1 = 1:length(x))
            ovrlp_tofix <- merge(ovrlp_tofix, x, by = "match1_list", all.x = T)
          } else {ovrlp_tofix$group1 <- NA}

          # fix instances
          x <- by(ovrlp_tofix, ovrlp_tofix$group1, function(xx) {
          
            # is there a single row with the highest overlap score?
            test <- length(which(xx$ovrlp_score == max(xx$ovrlp_score))) == 1
            
            # if yes, adopt that row
            if (test) {return(xx[which.max(xx$ovrlp_score), ])}
            
            # otherwise, adopt row with shortest ids (= no possible duplicates)
            if (! test) {
              return(xx[which.min(nchar(xx$match1_id) + nchar(xx$match2_id)),])}
          })
          ovrlp_fixed <- do.call(rbind, x)
          
          # rebind parts together
          ovrlp_yes <- rbind(ovrlp_ok, ovrlp_fixed[, colnames(ovrlp_ok)])
        }

      # next find any such instances in col 2 and separate them out for fixing
      ovrlp_yes$match2_list <- paste0(ovrlp_yes$match2_id, "x",
        substr(ovrlp_yes$match1_id, 1, 2))
      x <- unique(names(which(table(ovrlp_yes$match2_list) > 1)))
      x <- which(ovrlp_yes$match2_list %in% x)
      
        # if there are observations to fix...
        if (length(x) > 0) {
          
          # split dataset into part to fix and part that is ok
          ovrlp_tofix <- ovrlp_yes[x, ]
          ovrlp_ok <- ovrlp_yes[-x, ]
  
          # group instances
          x <- duplicated(ovrlp_tofix$match2_list)
          x <- unique(ovrlp_tofix[x, "match2_list"])
          if (length(x) != 0) {
            x <- data.frame(match2_list = x, group2 = 1:length(x))
            ovrlp_tofix <- merge(ovrlp_tofix, x, by = "match2_list", all.x = T)
          } else {ovrlp_tofix$group2 <- NA}

          # fix instances
          x <- by(ovrlp_tofix, ovrlp_tofix$group2, function(xx) {
          
            # is there a single row with the highest overlap score?
            test <- length(which(xx$ovrlp_score == max(xx$ovrlp_score))) == 1
            
            # if yes, adopt that row
            if (test) {return(xx[which.max(xx$ovrlp_score), ])}
            
            # otherwise, adopt row with shortest ids (= no possible duplicates)
            if (! test) {
              return(xx[which.min(nchar(xx$match1_id) + nchar(xx$match2_id)),])}
          })
          ovrlp_fixed <- do.call(rbind, x)
          
          # rebind parts together
          ovrlp_yes <- rbind(ovrlp_ok, ovrlp_fixed[, colnames(ovrlp_ok)])
        }
      
    # Restrict column names  
    ovrlp_yes <- ovrlp_yes[, c("match1_id", "match2_id")]      
           
  #...................................      
  ## Identify all overlapping triplets and pairs

    # Identify all overlap triplets
### THIS IS FINE BUT PROBABLY OVERKILL - FEWER LINES OF CODE SHOULD DO IT
    x <- by(ovrlp_yes, ovrlp_yes$match1_id, function(xx) {
      if (nrow(xx) > 1) {return(na.omit(c(unique(xx$match1_id), xx$match2_id)))}
    })
    x1 <- as.data.frame(do.call(rbind, x))
    colnames(x1) <- c("match1_id", "match2_id", "match3_id")
    x <- by(ovrlp_yes, ovrlp_yes$match2_id, function(xx) {
      if (nrow(xx) > 1) {return(na.omit(c(unique(xx$match2_id), xx$match1_id)))}
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
    
    # Skip forward if there are no triplets or pairs
    if (nrow(sets) == 0) {next}
    
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
      if (length(vars_mean) > 1) {
        merged[, vars_mean] <- colMeans(xx[, vars_mean], na.rm = T)
      }
      if (length(vars_mean) == 1) {
        merged[, vars_mean] <- mean(xx[vars_mean], na.rm = T)
      }
      
      # merge variables by sum
      vars_sum <- vars_ovrlp_f[which(vars_ovrlp_f$how == "sum"), "var"]
      if (length(vars_sum) > 1) {
        merged[, vars_sum] <- colsums(xx[, vars_sum], na.rm = T)
      }
      if (length(vars_sum) == 1) {
        merged[, vars_sum] <- sum(xx[vars_sum], na.rm = T)
      }
      
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
    if (exists("out_t") & exists("out_d")) {
      x <- colnames(out_t)
      df_out <- rbind(df_s[, x], out_t[, x], out_d[, x])
    }
    if (! exists("out_t") & exists("out_d")) {
      x <- colnames(out_d)
      df_out <- rbind(df_s[, x], out_d[, x])
    }
    if (exists("out_t") & ! exists("out_d")) {
      x <- colnames(out_t)
      df_out <- rbind(df_s[, x], out_t[, x])
    }
    if (! exists("out_t") & ! exists("out_d")) {
      df_out <- df_s[, as.character(vars_ovrlp_f$var)]
    }
    
    # Reconstruct date of death based on merged values
    x <- grep("day|month|year", colnames(df_out))
    for (i in x) {df_out[, x] <- round(df_out[, x], 0)}
    df_out$date_death <- as.Date(paste(df_out$year_death, df_out$month_death, 
      df_out$day_death, sep="-"), "%Y-%m-%d")
      
      # set dates to NA if they are after end of data collection
      df_out[which(df_out$date_death > date_end_f), "date_death"] <- NA
    
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
  
    # Exclusion criterion for insufficient information to analyse duplication
      # or overlap
    df_out$excl_del <- ifelse(df_out$del_score < 5, T, F)
    
    # Exclusion criterion for dates before analysis period
    df_out$excl_date <- ifelse(df_out$date_death < date_start_f, T, F)

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
