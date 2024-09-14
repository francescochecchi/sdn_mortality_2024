#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - GENERIC CODE ++++++ ###
#...............................................................................

#...............................................................................
## --------- BESPOKE FUNCTIONS TO IMPLEMENT DIFFERENT ANALYSIS TASKS  ------- ##
#...............................................................................


#...............................................................................
### Function for crude capture-recapture estimation based on a two-list system 
    # (m00 = x10 * x01 / x11)
    # Uses Chapman estimator given the expectation of low sample sizes:
    # https://en.wikipedia.org/wiki/Mark_and_recapture
#...............................................................................

f_chapman <- function(data_f, data_name_f, pars_f = pars) {
  
  #...................................      
  ## Preparatory steps   
    # Confirm number of lists
    n_lists <- pars_f[which(pars_f$loc_df == data_name_f), "n_lists"]
    
    # stop if fewer than three lists or more than four...
    if (n_lists != 2  | is.na(n_lists)) 
    {stop("wrong number of lists: only 2 lists allowed for two-list analysis")}
    
    # List names
    list_names <- pars_f[which(pars_f$loc_df == data_name_f), 
      paste0("list", 1:n_lists)]
    
    # Prepare contingency table
    df <- unlist(data_f[! names(data_f) %in% c("variable", "stratum")])
    names(df) <- names(data_f)[! names(data_f) %in% c("variable", "stratum")]
    x1 <- df["x10"] + df["x11"]
    x2 <- df["x01"] + df["x11"]
    x12 <- df["x11"]
  
  #...................................      
  ## Point estimate and confidence interval of total deaths, 
      # unlisted deaths and list sensitivity
  
    # Estimate of total deaths...
    total <- trunc(((x1 + 1) * (x2 + 1) / (x12 + 1)) - 1, 0)
    
    # 95% confidence interval of total deaths
    sigma_0.5 <- sqrt(1 / (x12 + 0.5) + 1 / (x2 - x12 + 0.5) + 1 / 
      (x1 - x12 + 0.5) + (x12 + 0.5) / ( (x1 - x12 + 0.5) * (x2 - x12 + 0.5) ) )
    total_lci <- x2 + x1 - x12 - 0.5 + ( (x2 - x12 + 0.5) * (x1 - x12 + 0.5) / 
      (x12 + 0.5) ) * exp(-1.96 * sigma_0.5)
    total_lci <- round(total_lci, 0)
    total_uci <- x2 + x1 - x12 - 0.5 + ( (x2 - x12 + 0.5) * (x1 - x12 + 0.5) / 
      (x12 + 0.5) ) * exp(1.96 * sigma_0.5)
    total_uci <- round(total_uci, 0)
    
    # Output
    out <- c(
      paste0(total - sum(df), " (95%CI ", total_lci - sum(df), " to ", 
        total_uci - sum(df), ")"),
      paste0(total, " (95%CI ", total_lci, " to ", total_uci, ")"),
      paste0(round(100 * x1 / total, 1), "% (95%CI ", 
        round(100 * x1 / total_uci, 1), "% to ", 
        round(100 * x1 / total_lci, 1), "%)"),
      paste0(round(100 * x2 / total, 1), "% (95%CI ", 
        round(100 * x2 / total_uci, 1), "% to ", 
        round(100 * x2 / total_lci, 1), "%)", sep = ""),
      paste0(round(100 * sum(df) / total, 1), "% (95%CI ", 
        round(100 * sum(df) / total_uci, 1), "% to ", 
        round(100 * sum(df) / total_lci, 1), "%)")
    )
    
    # Return output
    out <- cbind(c(
      "deaths outside any list", "total deaths", 
      paste0("sensitivity of ", list_names[1], "'s list"),
      paste0("sensitivity of ", list_names[2], "'s list"),
      "sensitivity of both lists combined"
    ), out )
    colnames(out) <- c("statistic", "estimate")
    return(out)
}


#...............................................................................
### Function to clean dataset for analysis
#...............................................................................

f_clean <- function(data_f, data_name_f, pars_f = pars, 
  f_clean_month_f = f_clean_month, f_clean_date_f = f_clean_date) {
  
  #...................................      
  ## Clean list variables
  
    # Name of the dataset
    df <- data_name_f
    
    # Figure out number of lists
    n_lists <- pars_f[which(pars_f$loc_df == df), "n_lists"]
    
    # Stop if fewer than two lists or more than four...
    if (n_lists < 2 | n_lists > 4 | is.na(n_lists)) {
      stop(paste("too few or too many lists: only 2, 3 or 4 lists allowed", 
      data_name_f, sep = " - "))
    }
    
    # Make sure each list variable contains only 1, 0 or NA
    for (i in 1:n_lists) {
      x <- paste0("list", i)
      data_f[, x] <- suppressWarnings(as.integer(data_f[, x]))
      data_f[! which(data_f[, x] %in% c(0, 1)), x] <- NA
    }
  
  #...................................      
  ## Clean age variable
    
    # Convert to integer
    data_f$age<- suppressWarnings(as.integer(data_f[, "age"]))
    
    # Warnings for implausible ages
    if (min(data_f[, "age"], na.rm = T) < 0) 
      {warning(paste("at least one age value is negative", df, sep = " - "))}
    if (max(data_f[, "age"], na.rm = T) > 100) 
      {warning(paste("at least one age value is > 100", df, sep = " - "))}
  
  #...................................      
  ## Clean gender variable
    
    # Standardise values
    data_f$gender <- tolower(as.character(data_f$gender))
    data_f[which(data_f$gender %in% c("male", "m")), "gender"] <- "male"
    data_f[which(data_f$gender %in% c("female", "fem", "f")),"gender"]<-"female"

  #...................................      
  ## Clean and generate dates for each death
  
    # Create new clean variable
    data_f$date_clean <- NA
    
    # Identify date / time variables
    date_cols <- colnames(data_f)[which(colnames(data_f) %in% 
      c("date_death", "month_death", "year_death"))]
    
    # If there are no date / time variables...
    if (length(date_cols) == 0) {
      warning(paste0(
        "no month, year or date of death variables found: check dataset ", df))
    }
    
    # If month is one of the variables, clean it to only feature integers (1-12)
    if ("month_death" %in% date_cols) {
      data_f$month_death <- f_clean_month_f(data_f$month_death)
    }
    
    # If year is one of the variables, integer it; # check for unusual values
    if ("year_death" %in% date_cols) {
      data_f$year_death <- suppressWarnings(as.integer(data_f$year_death))
      
      # warnings for implausible years
      if (min(data_f$year_death, na.rm = T) < year(Sys.Date() - 3650) ) {
        warning(paste0("at least one year value is >10y ago: check dataset",df))
      }
      if (max(data_f$year_death, na.rm = T) > year(Sys.Date()) ) {
        warning(paste0("at least one year value in the future: check dataset ", 
          df))
      }
    }
    
    # If year but not month is among the variables, create date from year 
        # (assume 15 June)
    if ((! "month_death" %in% date_cols) & ("year_death" %in% date_cols)) {
      data_f$date_clean <- as.Date(paste(data_f$year_death, 6, 15, sep = "-"))
    }
    
    # If month and year are among the variables, create date from these two
        # (assume 15th day; if month = NA, create date from year only, as above)
    if ("month_death" %in% date_cols & "year_death" %in% date_cols) {
      x <- which(is.na(data_f$month_death))
      data_f[x, "date_clean"] <- as.Date(paste(data_f$year_death, 
        6, 15, sep = "-"))
      data_f[-x, "date_clean"] <- as.Date(paste(data_f$year_death, 
        data_f$month_death, 15, sep = "-"))
    }
    
    # If date is among the variables, clean and override dates created 
        # from month and year
    if ("date_death" %in% data_cols) {
      
      # clean date
      data_f$date_death <- f_clean_date_f(data_f$date_death)

      # override
      x <- which(! is.na(data_f$date_death))
      data_f[x, "date_clean"] <- data_f[x, "date_death"]
    }
    
    # Make sure date is in the right format
    data_f$date_clean <- as.Date(data_f$date_clean)
    
    # Warnings for implausible dates
    if (min(data_f$date_clean, na.rm = T) < (Sys.Date() - 3650) ) {
      warning(paste0("at least one date is > 10y ago: check dataset ", df))
    }
    if (max(data_f$date_clean, na.rm = T) > Sys.Date()) { 
      warning(paste0("at least one date is in the future: check dataset ",df))
    }
  
  #...................................      
  ## Generate eligibility for each death (eligible if all lists have 1/0 value)
  data_f$eligible <- rowSums(is.na(data_f[, paste0("list", c(1:n_lists))]))
  data_f$eligible <- ifelse(data_f$eligible > 0, F, T)
  
  #...................................      
  ## Restrict dataset as desired
  
    # Restrict gender if desired
    x <- pars_f[which(pars_f$loc_df == df), "gender_restrict"]
    if (! is.na(x)) {data_f <- data_f[which(data_f$gender == x), ]}  
    
    # Restrict age if desired
      # minimum age    
      x <- as.integer(pars_f[which(pars_f$loc_df == df), "age_min"])
      if (! is.na(x)) {data_f <- data_f[which(data_f$age >= x), ]}  
      
      # maximum age    
      x <- as.integer(pars_f[which(pars_f$loc_df == df), "age_max"])
      if (! is.na(x)) {data_f <- data_f[which(data_f$age <= x), ]}  
    
    # Restrict analysis period if desired
      # starting date
      x <- f_clean_date_f(pars_f[which(pars_f$loc_df == df), "date_start"])
      if (! is.na(pars_f[which(pars_f$loc_df == df), "date_start"]) & is.na(x)){
        stop(paste0("date_start has been specified but is unreadable: check ", 
          df))        
      }
      if (! is.na(x)) {data_f <- data_f[which(data_f$date_clean >= x), ]}
      
      # ending date
      x <- f_clean_date_f(pars_f[which(pars_f$loc_df == df), "date_end"])
      if (! is.na(pars_f[which(pars_f$loc_df == df), "date_end"]) & is.na(x)){
        stop(paste0("date_end has been specified but is unreadable: check ", 
          df))        
      }
      if (! is.na(x)) {data_f <- data_f[which(data_f$date_clean >= x), ]}
  
  #...................................      
  ## Return prepared dataset
  return(data_f)
}


#...............................................................................
### Function to clean a vector of date values
    # will try to read the dates; if it can't, parse date as a character and
    # adopt the most common and plausible (lowest variance) formats
#...............................................................................

f_clean_date <- function(x) {

  #...................................      
  ## First try to read the date as a date
  out_date <- as.Date(x)

  #...................................      
  ## Try to detect all-digit formats by adopting the one that parses most
    
    # Convert to lower-case character, remove double spaces and non-alphanumeric
    x_conv <- trimws(gsub("[^[:alnum:]]", "", tolower(as.character(x))) )
  
    # Generate all options
    options <- c("%y%m%d", "%Y%m%d", "%y%d%m", "%Y%d%m", "%m%d%y", "%m%d%Y",
      "%d%m%y", "%d%m%Y")
    parses <- as.data.frame(sapply(options, function(xx) as.Date(x_conv, xx)))
    
    # Identify best option
    parses[! is.na(parses)] <- 1
    parses <- colSums(parses, na.rm = T)
    best_format_num <- names(which.max(parses))
    
    # Apply best option
    out_num <- as.Date(x_conv, best_format_num)

  #...................................      
  ## Try to detect formats with char month by adopting the one that parses most
    
    # Generate all options
    options <- c("%y%b%d", "%Y%b%d", "%y%d%b", "%Y%d%b", "%b%d%y", "%b%d%Y",
      "%d%b%y", "%d%b%Y", "%y%B%d", "%Y%B%d", "%y%d%B", "%Y%d%B", "%B%d%y", 
      "%B%d%Y", "%d%B%y", "%d%B%Y")
    parses <- as.data.frame(sapply(options, function(xx) as.Date(x_conv, xx)))
    
    # Set negative values to NA
    parses[parses < 0] <- NA
    
    # Compute variance of each option
    options_var <- scale(apply(parses, 2, var, na.rm = T), center = F)
      
    # Identify best option (most parses and lowest variance)
    parses[! is.na(parses)] <- 1
    parses <- colSums(parses, na.rm = T) - as.vector(options_var)
    best_format_chr <- names(which.max(parses))
      
    # Apply best option
    out_chr <- as.Date(x_conv, best_format_chr)
      
  #...................................      
  ## Return date values, by minimising the number of missing values
  
    # Bind all parses
    out <- data.frame(out_date, out_num, out_char)
    
    # Choose final value: date format if non-NA, if not numeric, otherwise char
    x <- which(is.na(out$out_date))
    out[-x, "out_final"] <- out[-x, "out_date"]
    out[x, "out_final"] <- out[x, "out_num"]
    x <- which(is.na(out$out_final))
    out[x, "out_final"] <- out[x, "out_char"]
    
  # Return date
  return(out$out_final)
}  


#...............................................................................
### Function to clean a vector of month values, returning 1-12 integers or NA
#...............................................................................

f_clean_month <- function(x) {
  
  # Initialise output
  out <- rep(NA, length(x))
    
  # Convert to lower-case character, remove double spaces and non-alphanumeric
  x_conv <- trimws(gsub("[^[:alnum:]]", "", tolower(as.character(x))) )
  
  # Parse values that are full or abbreviated month names 
    # (e.g. "jan", "MAR", "November", "april"..)  
  mmm <- grep(paste(tolower(month.abb), collapse = "|"), x_conv)
  out[mmm] <- match(substr(x_conv[mmm], 1, 3), tolower(month.abb))
  
  # Parse values that can be coerced to integers
  int <- suppressWarnings(which(! is.na(as.integer(x_conv))))
  out[int] <- as.integer(x_conv[int])

  # Return cleaned vector as 1:12 integers or NA values
  out[!out %in% 1:12] <- NA
  return(out)
}


#...............................................................................
### Function to describe patterns in the (clean) data and test for statistical 
    # differences
#...............................................................................

f_describe <- function(data_f, data_name_f, pars_f = pars, 
  palette_f = palette_gen, dir_path_f = dir_path) {
  
  #...................................      
  ## Preparatory steps 
    
    # Name of the dataset
    df <- data_name_f
    
    # Figure out number of lists and their names
    n_lists <- pars_f[which(pars_f$loc_df == df), "n_lists"]
    list_names <- pars_f[which(pars_f$loc_df == df),grep("list",colnames(pars_f))]
    
    # Remove ineligible observations
    if ("eligible" %in% colnames(data_f)) {
      data_f <- data_f[which(data_f$eligible == T), ]}
    
    # Reshape dataset long so as to recreate the 1 row = 1 list entry structure
        # before establishing overlap
    long <- data.frame()
    for (i in 1:n_lists) {
      x <- data_f[which(data_f[, paste0("list", i)] == 1), ]
      x <- x[, -grep("list", colnames(x))]
      x$list <- paste0("list", i)
      long <- rbind(long, x)
    }
    
    # Initialise descriptive output table
    tab <- as.data.frame(matrix(NA, nrow = 0, ncol = (4 + n_lists)))
    colnames(tab) <- c("characteristic", "category", list_names, "p_value" )      
  
  #...................................      
  ## Distribution of deaths by list and year (or month if time-span < two years)   
  
    # Which time unit to use for table
    t_unit <- ifelse(max(long$data_clean, na.rm = T) - 
      min(long$date_clean, na.rm = T) < 730, "month", "year")

    # Tabulate and add additional descriptive columns
    long$t_unit <- ifelse(t_unit == "year", year(long$date_clean),
      paste0(year(long$date_clean), "-", month(long$date_clean)) )
    x <- as.data.frame.matrix(table(long$t_unit, long$list))
    colnames(x) <- c("category", list_names)
    x[, c("characteristic", "p_value")] <- NA
    x$characteristic <- t_unit
    x[1, "p_value"] <- fisher.test(x[, list_names],simulate.p.value = T)$p.value

    # Add to main table
    tab <- rbind(tab, x[, colnames(tab)])
      
    # Warning on p-value
    if (nrow(x) > 5) {
      warning("best not to present p-value: table has too many categories")
    }  

  #...................................      
  ## Distribution of deaths by list and age category
  if ("age" %in% colnames(long)) {
    
    # Create age categories, which must include any chosen age cut-offs
      # identify cutoffs
      age_cuts <- pars_f[pars_f$loc_df == df, "age_cutoffs"]
      age_cuts <- c(0, as.integer(unlist(strsplit(as.character(age_cuts),","))),
        120)
     
      # create labels
      age_labs <- vector("character", length(age_cuts))
      for (i in 1:(length(age_labs) - 1)) {
        age_labs[i] <- paste0("age ", age_cuts[i], " to ", age_cuts[i + 1], "y")
      }

      # create categories
      long$age_cat <- cut(long$age, breaks = age_cuts, labels = age_labs,
        include.lowest = T, right = F)

    # Tabulate and add additional descriptive columns
    x <- as.data.frame.matrix(table(long$age_cat, long$list))
    colnames(x) <- c("category", list_names)
    x[, c("characteristic", "p_value")] <- NA
    x$characteristic <- "age (years)"

    # Add to main table
    tab <- rbind(tab, x[, colnames(tab)])
      
    # Add mean age by list, and p-value
    x <- aggregate(list(mean = long$age), by = list(list = long$list), 
      FUN = mean, na.rm = T)
    x <- c("age (years)", "mean", x$mean, 
      oneway.test(long$age ~ long$list)$p.value)
    tab <- rbind(tab, x)
  }
      
  #...................................      
  ## Distribution of deaths by list and gender
  if ("gender" %in% colnames(long)) {
    
    # Tabulate and add additional descriptive columns
    x <- as.data.frame.matrix(table(long$gender, long$list))
    colnames(x) <- c("category", list_names)
    x[, c("characteristic", "p_value")] <- NA
    x$characteristic <- "sex"
    x[1, "p_value"] <- fisher.test(x[, list_names],simulate.p.value = T)$p.value

    # Add to main table
    tab <- rbind(tab, x[, colnames(tab)])
  }
    
  #...................................      
  ## Write output table   
  write.csv(tab, paste0(dir_path_f, "out/", df, "_table_descr.csv"), 
    row.names = F, na = "")
  
  #...................................      
  ## Create and write plots 
    
    # Choose colours
    if (n_lists == 2) {colours <- palette_f[c(5,11)]}
    if (n_lists == 3) {colours <- palette_f[c(4,8,12)]}
    if (n_lists == 4) {colours <- palette_f[c(2,6,10,14)]}
    chars <- unique(tab$characteristic)
  
    # For each characteristic of interest...
    for (i in c("year", "month", "age", "sex")) {
      if (grepl(i, chars)) {
      
        # reshape table for graphing    
        x <- subset(tab, characteristic == i)
        x <- subset(tab, category != "mean")
        x <- reshape(x, direction = "long", varying = list_names,
          idvar = "category", timevar = "list", times = list_names,
          v.names = "frequency", drop = c("characteristic", "p_value"))
  
        # plot
        pl <- ggplot(x, aes(y = frequency, x = category, 
          group = list, colour = list, fill = list)) +
          geom_bar(position = "dodge", stat = "identity", alpha = 0.4) +
          scale_y_continuous("number") +
          scale_x_discrete(i) +
          theme_bw() +
          scale_color_manual("list", values = colours) +
          scale_fill_manual("list", values = colours) +
          theme(legend.position = "top")
        if (i == time_unit | i = "age") {
          pl <- pl + theme(axis.text.x = element_text(angle = 30, hjust = 1, 
            vjust = 1))
        }  
        ggsave(paste0(dir_path, "out/", df, "_by_list__and_", i, ".png"),  
          width = 15, height = 10, units = "cm", dpi = "print")
        
        # name plot
        assign(paste0("pl_", i), pl)
      }
    }
    
    # Combined graph
    plot_list <- list(mget(grep("pl_", names(.GlobalEnv), value = T)))
    plot_combi <- ggpubr::ggarrange(plotlist = plot_list, ncol = 1)  
    ggsave(paste0(dir_path_f, "out/", df, "_descr_combi_long.png"), 
      width = 15, height = 20, units = "cm", dpi = "print")
    ggsave(paste0(dir_path_f, "out/", df, "_descr_combi_wide.png"), 
      width = 20, height = 13, units = "cm", dpi = "print")
}  


#...............................................................................
### Function to fit each candidate log-linear model for a 3- or 4-list system
    # as per Rossi et al. https://rivista-statistica.unibo.it/article/view/9854 
############ WORKS, BUT TO BE REWRITTEN
#...............................................................................

f_logl <- function(data_f, data_name_f, pars_f = pars) {
  
  #...................................      
  ## Preparatory steps
  
  # Confirm number of lists
  n_lists <- pars_f[pars_f[, "loc_df"] == data_name_f, "n_lists"]
  # stop if fewer than three lists or more than four...
  if (! n_lists %in% c(3, 4)  | is.na(n_lists)) 
  {stop("wrong number of lists: only 3 or 4 lists allowed for a three- or four-list analysis")}
  
  # list names
  list_names <- pars_f[pars_f[, "loc_df"] == data_name_f, paste("list", 1:n_lists, sep = "")]
  
  # Identify exposure and confounders    
  exposure <- tolower(trimws(pars_f[pars_f[, "loc_df"] == data_name_f, "exposure"]))
  confounders <- pars_f[pars_f[, "loc_df"] == data_name_f, "confounders"]
  if (! is.na(confounders[1])) {confounders <- tolower(sapply(unlist(strsplit(confounders, ",")), trimws)) }
  
  # Prepare dataset
  # restrict data to eligible observations
  if ("eligible" %in% colnames(data_f)) {data_f <- subset(data_f, eligible == T)}
  
  # create unique id for each observation (ignore any existing id variable)
  data_f[, "key"] <- paste("id", 1:nrow(data_f), sep = "")
  
  # create indicator variable for outcome (all = 1)
  data_f[, "y"] <- 1
  
  # rename list columns
  for (i in 1:n_lists) {colnames(data_f)[colnames(data_f) == paste("list", i, sep = "")] <- paste("x", i, sep = "")}
  
  # defactor exposure and confounders
  if (! is.na(exposure) ) {
    if (is.factor(data_f[, exposure]) ) {data_f[, exposure] <- as.character(data_f[, exposure]) }
  }  
  if (! is.na(confounders[1]) ) {
    for (i in confounders) {
      if (is.factor(data_f[, i]) ) {data_f[, i] <- as.character(data_f[, i]) }
    }
  }
  
  #...................................      
  ## Prepare data for modelling (note: for individual-level analysis, i.e. one row = one observation)
  
  # Prepare possible 'profiles' for each individual observation
  # identify all possible list profiles (cells in contingency table)
  profiles <- as.data.frame(permutations(n = 2, r = n_lists, v = c(0, 1), repeats.allowed = T))
  colnames(profiles) <- as.character(1:n_lists)
  
  # columns for which observations appear on which combinations of lists
  for (i in 3:n_lists) {
    x1 <- combinations(n = n_lists, r = i - 1, v = 1:n_lists, repeats.allowed = F)
    for (j in 1:nrow(x1)) { profiles[, paste(x1[j, ], collapse = "")] <- rowSums(profiles[, x1[j, ]]) }
    profiles[, nchar(colnames(profiles)) == (i - 1)] <- ifelse(profiles[, nchar(colnames(profiles)) == (i - 1)] == (i - 1), 1, 0)
  }
  colnames(profiles) <- paste("x", colnames(profiles), sep = "")
  x2 <- colnames(profiles)
  
  # Set profiles for each observation, based on the data
  # create expanded dataframe with all possible profiles for each observation
  df <- expand_grid(data_f[, "key"], profiles)
  colnames(df) <- c("key", colnames(profiles))
  df <- as.data.frame(df)
  
  # match each observation to possible profiles (except for which lists the observation is in)
  x1 <- paste("x", 1:n_lists, sep = "")
  df <- merge(df, data_f[, ! colnames(data_f) %in% c(x1, "y")], by = "key", all.x = T)
  
  # determine which profile the observation has, based on which lists the observation is in
  # outcome indicator = 1 if observation falls within a given profile, 0 otherwise and NA for x000(0) profile)
  df <- merge(df, data_f[, c("key", x1, "y")], by = c("key", x1), all.x = T)
  df[, "y"] <- ifelse(is.na(df[, "y"]), 0, df[, "y"])
  df[rowSums(df[, x1]) == 0, "y"] <- NA
  
  # Add columns for interactions between exposure and confounder (if present) and lists
  # (note: omit interactions among confounders and between exposure and confounders)
  # profiles:exposure interactions
  if (! is.na(exposure)) { 
    df[, paste(x2, exposure, sep = ":")] <- ifelse(df[, x2] == 1, df[, exposure], 0)
  }
  
  # profiles:confounder(s) interactions
  if (! is.na(confounders[1])) {
    for (i in confounders) {df[, paste(x2, i, sep = ":")] <- ifelse(df[, x2] == 1, df[, i], 0)} 
  }
  
  
  #...................................      
  ## Define candidate models
  # Possible terms...
  x1 <- c()
  for (i in 3:n_lists) {
    x1 <- c(x1, apply( combinations(n_lists, i - 1, unlist(list_names) ) , 1 , paste , collapse = " x " ) )
  }
  # ...and their length
  x2 <- lapply(sapply(x1, strsplit, " x "), length)
  
  # Combinations of two-list terms
  x3 <- names(x2[x2 == 2])
  x4 <- c()
  for (i in n_lists:1 ) {
    x4 <- c(x4, apply(combinations(length(x3), i, x3) , 1 , paste , collapse = ", " ) )
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
        x5 <- c(x5, apply(combinations(length(x3), i, x3) , 1 , paste , collapse = ", " ) )
      }
    }
    
    # matrix to indicate overlap of two- and three-list term combinations (overlap = violation of hierarchy principle)
    x7 <- rep(NA, length(x5))
    for (i in 1:length(x5) ) {
      x8 <- unlist(strsplit(x5[i], ", "))
      x9 <- c()
      for (j in 1:length(x8)) {
        x10 <- unlist(strsplit(x8[j], " x ") )
        x9 <- c(x9, apply(permutations(length(x10), 3, x10), 1, paste, collapse = " x ") )
      }
      x7[i] <- paste(x9, collapse = ", ")
    }
    
    x6 <- as.data.frame(matrix(sapply(gsub(", ", "|", x4), grepl, x7), nrow = length(x5), ncol = length(x4)))
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
  for (i in 1:length(list_names) ) {
    x1 <- gsub(list_names[i], gsub("list", "", names(list_names[i]) ), x1)
  }
  x1 <- suppressWarnings(do.call(rbind, strsplit(x1, ", ")))
  # this will result in elements being recycled because of the uneven n of resulting columns, which the next few lines deal with
  x2 <- apply(x1, c(1, 2), function(x) {
    x3 <- strsplit(x, " x "); x3 <- sort(unlist(x3)); return(paste(x3, collapse = ""))
  })
  x3 <- apply(x2, 1, function(x) {unique(paste("x", x, sep = ""))})
  x3 <- sapply(x3, paste, collapse = " + ")
  x1 <- paste(paste("x", 1:n_lists, sep = ""), collapse = " + ")
  x2 <- paste("y ~  ", x1, " + ", x3, sep = "")
  out[, "formula"] <- x2
  out[1, "formula"] <- paste("y ~  ", x1, sep = "")
  
  # Lastly, add terms for exposure, confounders and interactions of these with profiles
  # profiles:exposure interactions
  if (! is.na(exposure)) { 
    for (i in 1:nrow(out)) {  
      x1 <- all.vars(as.formula(out[i, "formula"]))[-1]
      x1 <- paste(x1, exposure, sep = ":")
      x1 <- paste(x1, collapse = " + ")
      out[i, "formula"] <- paste(out[i, "formula"], exposure, x1, sep = " + ")
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
      out[, paste(unlisted, "screen", levels(as.factor(df[, exposure])), sep = "_")] <- NA
    }
  }
  
  # unlisted observations for model without each of the confounders, or no adjustment for confounders
  if (! is.na(confounders[1])) { out[, paste(unlisted, "without", c(confounders, "adjustment"), sep = "_")] <- NA }
  
  # other statistics  
  out[, "lrt_p"] <- NA
  
  # Fit saturated model (needed to perform likelihood-ratio test for models nested within it)
  
  # identify saturated model
  x1 <- lapply(out[, "model"], function(x) {unlist(strsplit(x, ","))} )
  x2 <- which.max(unlist(lapply(x1, length)))
  
  # fit saturated model
  fit_sat <- try(glm(formula = as.formula(out[x2, "formula"]), family = "poisson", data = df, maxit = 1000 ), 
                 silent = T)
  
  # Fit all other models        
  for (i in 1:nrow(out) ) {
    print(paste("now fitting candidate model  ", i, "of", nrow(out), sep = " ") )
    
    # fit model (or at least try)
    suppressWarnings(rm(fit) )
    fit <- try(glm(formula = as.formula(out[i, "formula"]), family = "poisson", data = df, maxit = 1000 ), 
               silent = T)
    
    # check: if model has not fit or any of the coefficients is NA, skip to next loop
    if (class("fit")[1] == "try-error" ) {next}
    if (any(is.na(coef(fit))) ) {next}
    
    # estimate expected number unlisted
    # select data with 000[0] profile
    df0 <- subset(df, is.na(y))
    
    # probability of being unlisted for each observation
    df0[, paste(unlisted, "screen", sep = "_")] <- predict(fit, df0)
    
    # number unlisted overall and by exposure level
    out[i, paste(unlisted, "screen", sep = "_")] <- round(sum(exp(na.omit(df0[, paste(unlisted, "screen", sep = "_")]))), digits = 0)
    if (! is.na(exposure)) {
      if( length(levels(as.factor(df[, exposure]))) < 10 ) {
        out[i, paste(unlisted, "screen", levels(as.factor(df[, exposure])), sep = "_")] <- 
          aggregate(df0[, paste(unlisted, "screen", sep = "_")],
                    by = list(as.factor(df0[, exposure])), FUN = function(x) {round(sum(exp(na.omit(x))), digits = 0)})[, 2]
      }
    }  
    
    # number unlisted if each confounder is taken out, or all are taken out of the model
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
        out[i, paste(unlisted, "without", j, sep = "_")] <- round(sum(exp(na.omit(x4))), digits = 0)
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
      out[i, paste(unlisted, "without_adjustment", sep = "_")] <- round(sum(exp(na.omit(x4))), digits = 0)
    }  
    
    # Likelihood ratio test p-value comparing model to saturated model (low p = model is better than saturated model)
    x1 <- -2 * (logLik(fit) - logLik(fit_sat))
    out[i, "lrt_p"] <- as.numeric(pchisq(x1, df = fit$df.residual - fit_sat$df.residual, lower.tail = F))
    
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

f_model_average <- function(f_out, data_name_f, pars_f = pars) {
  
  #...................................      
  ## Preparatory steps
  
  # Confirm number of lists
  n_lists <- pars_f[pars_f[, "loc_df"] == data_name_f, "n_lists"]
  # stop if fewer than three lists or more than four...
  if (! n_lists %in% c(3, 4)  | is.na(n_lists)) 
  {stop("wrong number of lists: only 3 or 4 lists allowed for a three- or four-list analysis")}
  
  # list names
  list_names <- pars_f[pars_f[, "loc_df"] == data_name_f, paste("list", 1:n_lists, sep = "")]
  
  # Identify exposure  
  exposure <- tolower(trimws(pars_f[pars_f[, "loc_df"] == data_name_f, "exposure"]))
  
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
    if (out[i, paste(unlisted, "screen", sep = "_")] / sum(df[, "y"], na.rm = T) > pars_f[pars_f[, "loc_df"] == data_name_f, "plausibility_1"]) 
    {out[i, "eligible"] <- "no - implausible estimate"}
    if (out[i, "lrt_p"] > pars_f[pars_f[, "loc_df"] == data_name_f, "plausibility_2"]) 
    {out[i, "eligible"] <- "no - possible over-fitting"}
  }
  
  #...................................      
  ## Re-fit each eligible model, compute AIC and predict n unlisted
  for (i in 1:nrow(out)) {
    if (out[i, "eligible"] == "yes") {
      
      # refit model      
      fit <- try(glm(formula = as.formula(out[i, "formula"]), family = "poisson", data = df, maxit = 1000 ), 
                 silent = T)
      
      # compute model's AIC
      out[i, "aic"] <- round(AIC(fit), digits = 2)
      
      # predict unlisted observations
      # contribution to total unlisted for each observation - point estimate and 95%CI
      x1 <- predict(fit, df0, se.fit = T)
      df0[, unlisted] <- exp(x1[[1]])
      df0[, paste(unlisted, "lci", sep = "_")] <- exp(x1[[1]] - 1.96 * x1[[2]])
      df0[, paste(unlisted, "uci", sep = "_")] <- exp(x1[[1]] + 1.96 * x1[[2]])
      
      # estimated number unlisted overall... 
      out[i, unlisted] <- round(sum(na.omit(df0[, unlisted])), digits = 0)
      out[i, paste(unlisted, "lci", sep = "_")] <- round(sum(na.omit(df0[, paste(unlisted, "lci", sep = "_")])), digits = 0)
      out[i, paste(unlisted, "uci", sep = "_")] <- round(sum(na.omit(df0[, paste(unlisted, "uci", sep = "_")])), digits = 0)
      
      # ...and by exposure level
      if (! is.na(exposure)) { 
        if (length(levels(as.factor(df[, exposure]))) < 10 ) {
          out[i, paste(unlisted, levels(as.factor(df[, exposure])), sep = "_")] <- aggregate(df0[, unlisted],
                                                                                             by = list(as.factor(df0[, exposure])), FUN = function(x) {round(sum(na.omit(x)), digits = 0)})[, 2]
          out[i, paste(unlisted, levels(as.factor(df[, exposure])), "lci", sep = "_")] <- aggregate(df0[, paste(unlisted, "lci", sep = "_")],
                                                                                                    by = list(as.factor(df0[, exposure])), FUN = function(x) {round(sum(na.omit(x)), digits = 0)})[, 2]
          out[i, paste(unlisted, levels(as.factor(df[, exposure])), "uci", sep = "_")] <- aggregate(df0[, paste(unlisted, "uci", sep = "_")],
                                                                                                    by = list(as.factor(df0[, exposure])), FUN = function(x) {round(sum(na.omit(x)), digits = 0)})[, 2]
        }
      }
    }
  }
  
  #...................................      
  ## Calculate 'posterior probabilities' (weights) of eligible models from their AICs
  # (based on Rossi, 2010: https://rivista-statistica.unibo.it/article/view/3593/2945 )
  
  # Calculate an AIC delta based on lowest one
  out[, "aic_delta"] <- NA
  out[out[, "eligible"] == "yes", "aic_delta"] <- out[out[, "eligible"] == "yes", "aic"] - min(out[out[, "eligible"] == "yes", "aic"])
  
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
        out_pretty[, paste("deaths outside any list (95%CI)", i, sep = " - ")] <- paste(
          out_pretty[, paste(unlisted, i, sep = "_")], " (", 
          out_pretty[, paste(unlisted, i, "lci", sep = "_")], " to ", 
          out_pretty[, paste(unlisted, i, "uci", sep = "_")], ")", sep = "")
      }  
    }
  }
  out_pretty[, "likelihood ratio p-value"] <- out_pretty[, "lrt_p"]
  out_pretty[, "AIC"] <- round(out_pretty[, "aic"], 2)
  out_pretty[, "posterior probability"] <- round(out_pretty[, "post_prob"], 3)
  out_pretty <- out_pretty[, c("model", grep("outside", colnames(out_pretty), value = T),
                               "likelihood ratio p-value", "AIC", "posterior probability")]
  out_pretty <- cbind(rep("", nrow(out_pretty)), out_pretty)
  colnames(out_pretty)[1] <- "-"
  
  # Estimated unlisted and total deaths with 95%CIs
  out_est_raw <-c(
    weighted.mean(out[, unlisted], out$post_prob, na.rm = T),
    weighted.mean(out[, paste(unlisted, "lci", sep = "_")], out$post_prob, na.rm = T),
    weighted.mean(out[, paste(unlisted, "uci", sep = "_")], out$post_prob, na.rm = T)
  )
  out_est_raw <-c(out_est_raw, sum(df[, "y"], na.rm = T) + out_est_raw)
  if (! is.na(exposure)) { 
    if (length(levels(as.factor(df[, exposure]))) < 10 ) {
      for (i in levels(as.factor(df[, exposure])) ) {
        x1 <- c(weighted.mean(out[, paste(unlisted, i, sep = "_")], out$post_prob, na.rm = T),
                weighted.mean(out[, paste(unlisted, i, "lci", sep = "_")], out$post_prob, na.rm = T),
                weighted.mean(out[, paste(unlisted, i, "uci", sep = "_")], out$post_prob, na.rm = T)
        )
        out_est_raw <- rbind(out_est_raw, c(x1, sum(df[which(df[, exposure] == i), "y"], na.rm = T) + x1)
        )
      }  
    }
  }
  
  out_est_raw <- as.data.frame(rbind(out_est_raw))
  out_est_raw <- round(out_est_raw, 0)
  out_est_raw <- cbind(rep(NA, nrow(out_est_raw)), out_est_raw)
  colnames(out_est_raw) <- c("stratum", "unlisted", "unlisted_lci", "unlisted_uci", 
                             "total_deaths_est", "total_deaths_lci", "total_deaths_uci")
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
  
  out_est_pretty <- as.data.frame(matrix(NA, nrow = nrow(out_est_raw), ncol = 3))
  for (i in 1:nrow(out_est_raw)) {
    out_est_pretty[i, 1] <- out_est_raw[i, 2]
    out_est_pretty[i, 2] <- paste(out_est_raw[i, 3], " (", out_est_raw[i, 4], " to ", out_est_raw[i, 5], ")", sep = "")
    out_est_pretty[i, 3] <- paste(out_est_raw[i, 6], " (", out_est_raw[i, 7], " to ", out_est_raw[i, 8], ")", sep = "") 
  }
  colnames(out_est_pretty) <- c("stratum", "deaths outside any list (95%CI)", "total deaths (95%CI)")
  out_est_pretty <- cbind(rep("", nrow(out_est_pretty)), out_est_pretty)
  colnames(out_est_pretty)[1] <- "-"
  
  # Sensitivity of each list and all lists combined
  out_sens_raw <- as.data.frame(matrix(NA, ncol = 5, nrow = n_lists + 1))
  colnames(out_sens_raw) <- c("list", "n_deaths", "sens_est", "sens_lci", "sens_uci")
  out_sens_raw[, "list"] <- c(unlist(list_names), "all lists")
  
  for (i in 1:n_lists) {
    out_sens_raw[i, "n_deaths"] <- sum(df[df[, grep(as.character(i), colnames(df))] == 1, "y"], na.rm = T)
  }
  
  # Calculate sensitivity for all lists combined
  out_sens_raw[n_lists + 1, "n_deaths"] <- sum(out_sens_raw[1:n_lists, "n_deaths"], na.rm = T)
  
  # Calculate sensitivity estimates, LCI, and UCI
  out_sens_raw[, "sens_est"] <- out_sens_raw[, "n_deaths"] / out_est_raw[1, "total_deaths_est"]
  out_sens_raw[, "sens_lci"] <- out_sens_raw[, "n_deaths"] / out_est_raw[1, "total_deaths_uci"]
  out_sens_raw[, "sens_uci"] <- out_sens_raw[, "n_deaths"] / out_est_raw[1, "total_deaths_lci"]
  
  out_sens_raw <- cbind(rep("", nrow(out_sens_raw)), out_sens_raw)
  colnames(out_sens_raw)[1] <- "-"
  
  out_sens_pretty <- out_sens_raw
  out_sens_pretty[, "number of deaths"] <- out_sens_pretty[, "n_deaths"]
  out_sens_pretty[, "sensitivity (95%CI)"] <- paste(round(out_sens_pretty[, "sens_est"] * 100, 1), "% (",
                                                    round(out_sens_pretty[, "sens_lci"] * 100, 1), "% to ", round(out_sens_pretty[, "sens_uci"] * 100, 1), "%)", sep = "")
  out_sens_pretty <- out_sens_pretty[, c("list", "number of deaths", "sensitivity (95%CI)")]
  out_sens_pretty <- cbind(rep("", nrow(out_sens_pretty)), out_sens_pretty)
  colnames(out_sens_pretty)[1] <- "-"
  
  
  # Return outputs as a list
  x1 <- list("out_raw" = out_raw, "out_pretty" = out_pretty, "out_est_raw" = out_est_raw, 
             "out_est_pretty" = out_est_pretty, "out_sens_raw" = out_sens_raw, "out_sens_pretty" = out_sens_pretty)
  return(x1)
  
}



#...............................................................................
### Function to tabulate overlap among lists, create contingency table in vector 
    # form and draw Venn diagram for entire dataset and desired strata
############### WORKS BUT TO BE REWRITTEN
#...............................................................................

f_overlap <- function(data_f, data_name_f, pars_f = pars, 
  palette_f = palette_gen, dir_path_f = dir_path) {
  
  #...................................      
  ## Preparatory steps
  
    # Name of the dataset
    df <- data_name_f
    
    # Figure out number of lists
    n_lists <- pars_f[which(pars_f$loc_df == df), "n_lists"]
    
    # Stop if fewer than two lists or more than four...
    if (n_lists < 2 | n_lists > 4 | is.na(n_lists)) {
      stop(paste0("too few or too many lists: only 2, 3 or 4 lists allowed", 
      ": check dataset", df))
    }
 
    # Remove ineligible observations
    if ("eligible" %in% colnames(data_f)) {
      data_f <- data_f[which(data_f$eligible == T), ]}
  
    # Set up output overlap vectors for overall dataset and each stratum
      # identify strata
      strata <- grep(
        paste(c("gender_cat", "age_cat", "period_cat"), collapse = "|"), 
        colnames(data_f), value = T)
      
      # initialise output
      out <- data.frame(
        variable = c("overall", sapply(strata, function (x) 
          {rep(x, length(levels(data_f[, x])))})),
        stratum = c("overall", sapply(strata, function (x) 
          {levels(data_f[, x])})),
      )

    # Create n-list contingency table
    x <- as.data.frame(permutations(n = 2, r = n_lists, v = c(0, 1), 
      repeats.allowed = T))
    x <- x[! rowSums(x) == 0, ]
    ctab <- apply(cbind("x", x), 1, paste, collapse = "")
    tab <- as.data.frame(matrix(NA, ncol = n_lists + 1, nrow = length(ctab)))
    colnames(tab) <- c("permutation", paste("list", 1:n_lists, sep = "") )
    tab$permutation <- ctab
    tab[, grep("list", colnames(tab))] <- x
    
    # contingency cells in output
    out[, tab$premutation] <- NA
    
    # Choose colours
    if (n_lists == 2) {colours <- palette_f[c(5,11)]}
    if (n_lists == 3) {colours <- palette_f[c(4,8,12)]}
    if (n_lists == 4) {colours <- palette_f[c(2,6,10,14)]}
  
  #...................................      
  ## Compute contingency table and draw Venn diagrams 
      # for all observations and each desired stratum
  for (i in 1:nrow(out)) {
    
    # Select data
    if (out[i, "variable"] == "overall") {data_i <- data_f}
    if (out[i, "variable"] != "overall") {
      data_i <- data_f[which(data_f[,out[i,"variable"]] == out[i, "stratum"]),]}
    
    # Contingency table
      # reset contingency table
      tab <- tab[, colnames(tab) != "Freq"]
    
      # tabulate overlap among lists
      x <- "~ list1"
      for (j in 2:n_lists) {x <- as.formula(paste(x,paste0("list", j),sep="+"))}
      x <- as.data.frame(xtabs(x, data = data_i))
      
      # add to contingency table
      tab <- merge(tab, x, by = paste0("list", 1:n_lists), all.x = T)
      tab[which(is.na(tab$Freq)), "Freq"] <- 0
      
      # transfer to output
      out[i, tab$permutation] <- tab$Freq
    
    # Draw and save Venn diagram
      # create list of sets
      data_i[, "unique_id"] <- paste0("id", 1:nrow(data_i))
      sets <- list()
      for (j in 1:n_lists) {sets[[j]] <- 
        data_i[which(data_i[, paste0("list", j)] == 1), "unique_id"]}
      names(sets) <- list_names
      
      # plot and save
      pl <- ggvenn(sets, fill_color = colours, fill_alpha = 0.3, 
        show_percentage = F, stroke_color = "grey50", 
        stroke_alpha = 0.7, stroke_size = 1, set_name_color = "black", 
        set_name_size = 4, text_color = "black", text_size = 4)
      ggsave(paste0(dir_path_f, "out/", df, "_venn_", 
        gsub("-", "_", paste0(out[i, "stratum"]) ), ".png"),
        width = 25, height = 15, units = "cm", dpi = "print")
      
      # assign plot name
      assign(paste0("venn_", gsub("-", "_", paste(out[i, "stratum"]))), pl)
  }
  
  #...................................      
  ## Format and save output
  
    # Create combined plots of all Venn diagrams
    x <- gsub("-", "_", paste("venn_", out$stratum))
    plot_list <- mget(x)
    plot_labels <- gsub("venn_", "", x)
    plot_labels <- gsub("_", " ", plot_labels)
    
    # Fix date labels
    x <- substr(plot_labels, 1, 1) %in% c(1:2)
    for (i in plot_labels[x]) {
      plot_labels[plot_labels == i] <- paste0(substr(i,9,10), "/",substr(i,6,7), 
        "/",  substr(i,3,4), " to ", substr(i,23,24), "/", substr(i,20,21), "/",
        substr(i,17,18)) 
    }
    
    pl <- ggpubr::ggarrange(plot_list, labels = plot_labels, 
      hjust = 0, vjust = 1, font.label = list(size = 10, face = "plain"))  
    ggsave(paste0(dir_path_f, "out/", ff, "_venn_", "combi_long.png"), 
      width = 18, height = 25, units = "cm", dpi = "print")
    ggsave(paste0(dir_path_f, "out/", ff, "_venn_", "combi_wide.png"), 
      width = 30, height = 20, units = "cm", dpi = "print")
    
    # Return and save output vector
    write.csv(out, paste0(dir_path, "out/",df,"_table_overlap.csv"),row.names=F)
    return(out)
}  


#...............................................................................
### Function to prepare datasets for gender, period and age strata, as desired
#...............................................................................

f_stratify <- function(data_f, data_name_f, pars_f = pars, 
  f_clean_date_f = f_clean_date) {
  
  #...................................      
  ## Preparatory steps  
  
    # Name of the dataset
    df <- data_name_f
    
  #...................................      
  ## Prepare gender stratum variable        
  if (pars_f[which(pars_f$loc_df == df), "gender"] == "Y") {
    
    # Factor gender
    data_f[, "gender_stratum"] <- factor(data_f[, "gender"], 
      levels = c("female", "male"))
    
    # Stop if the dataset contains only one gender
    if (length(unique(data_f[, "gender"])) < 2 ) {
      stop("cannot stratify by gender: data contain only one gender")
    }  
  }
  
  #...................................      
  ## Prepare age stratum variable
  age_cuts <- pars_f[pars_f$loc_df == df, "age_cutoffs"]
  if (! is.na(age_cuts)) {
    
    # Define cutoffs
    age_cuts <- c(0, as.integer(unlist(strsplit(as.character(age_cuts),","))),
      120)
    
    # Stop if the dataset doesn't feature all the age strata
    if (min(data_f$age, na.rm = T) > age_cuts[2] |  
        max(data_f$age, na.rm = T) < age_cuts[length(age_cuts) - 1]) {
      stop("cannot stratify by age: some desired age strata missing in data")
    }
    
    # Create labels
    age_labs <- vector("character", length(age_cuts))
    for (i in 1:(length(age_labs) - 1)) {
      age_labs[i] <- paste0("age ", age_cuts[i], " to ", age_cuts[i + 1], "y")
    }
  
    # Create categories
    data_f$age_cat <- cut(data_f$age, breaks = age_cuts, labels = age_labs,
      include.lowest = T, right = F)
  }

  #...................................      
  ## Prepare period stratum variable     
  period_cuts <- pars_f[pars_f$loc_df == df, "period_cutoffs"]
  if (! is.na(period_cuts)) {

    # Define cutoffs
    period_cuts <- c(min(data_f$date_clean), 
      f_clean_date_f(unlist(strsplit(as.character(period_cuts), ","))),
      max(data_f$date_clean) + 1)
    
    # Stop if the dataset doesn't feature all the period strata
    if (min(data_f$date_clean, na.rm = T) > period_cuts[2] |  
        max(data_f$date_clean, na.rm = T) < period_cuts[length(period_cuts-1)]){
      stop("cannot stratify by period: some desired dates missing in data")
    }
    
    # Create labels
    period_labs <- vector("character", length(period_cuts))
    for (i in 1:(length(period_labs) - 1)) {
      period_labs[i] <- paste0(period_cuts[i], " to ", period_cuts[i + 1])
    }
  
    # Create categories
    data_f$period_cat <- cut(data_f$date_clean, breaks = period_cuts, 
      labels = period_labs, include.lowest = T, right = F)
  }
    
  #...................................      
  ## Return prepared dataset
  return(data_f)
}

#...............................................................................
### ENDS 
#...............................................................................

