#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - SUDAN (2024) ++++++ ###
#...............................................................................

#...............................................................................
## ----------- R SCRIPT TO PRODUCE DESCRIPTIVE ANALYSES OF DATASET  --------- ##
#...............................................................................


#...............................................................................
### Preparing the dataset and tabulating attrition
#...............................................................................

  #...................................      
  ## Read and prepare dataset
    
    # Which directory
    dir_ij <- paste0(dir_path, "out/d3o3/")
    
    # Read dataset arising from given duplication and overlap thresholds
    df <- readRDS(paste0(dir_ij, "list_data.rds"))

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
    df[which(is.na(df$excl_date_imp)), "excl_date_imp"] <- F
    
    # Add/modify missing categories for different variables
    df[which((is.na(df$gender))), "gender"] <- "missing"
    df[which((is.na(df$age_cat))), "age_cat"] <- "missing"
    df[which((is.na(df$age_cat_imp))), "age_cat_imp"] <- "missing"
    df[which((is.na(df$resistance_committees))), "resistance_committees"] <- F
    df$resistance_committees <- ifelse(df$resistance_committees, "yes", 
      "no / unknown")
    df[which((is.na(df$year_death))), "year_death"] <- "missing"
    df[which((is.na(df$year_death_imp))), "year_death_imp"] <- "missing"
    
    
  #...................................      
  ## Tabulate attrition
    
    # Initialise table
    tab <- data.frame(
      criterion = c(
        "total observations",
        "sufficient unique identifiers", 
        "died since 15 April 2023", "died within Sudan (descriptive analysis)",
        "died within Khartoum State (estimation - all causes)", 
        "died within Khartoum State (estimation - intentional injury only)")
    )
    tab[, c("number", "percentage")] <- NA
    
    # Populate number column
    tab$number <- c(
      nrow(df), 
      nrow(df[which(!df$excl_del), ]),
      nrow(df[which(!df$excl_del & !df$excl_date), ]),
      nrow(df[which(!df$excl_del & !df$excl_date & !df$excl_loc_death), ]),
      nrow(df[which(!df$excl_del & !df$excl_date & !df$excl_loc_death &
        !df$excl_kht), ]),
      nrow(df[which(!df$excl_del & !df$excl_date & !df$excl_loc_death &
        !df$excl_kht & !df$excl_cod), ])
    )
    
    # Populate percent column
    tab$percentage <- tab$number / tab[1, "number"]
    tab$percentage <- percent(tab$percentage, accuracy = 0.1)
    
    # Save table
    write.csv(tab, paste0(dir_ij, "attrition.csv"), row.names = F)    

    
#...............................................................................
### Describing and contrasting the different lists
#...............................................................................
    
  #...................................      
  ## Reconstitute 'long' dataset with one row = one record (before matching)
    
    # Apply exclusion criteria for descriptive analysis
    df_desc <- df[which(!df$excl_del & !df$excl_date & !df$excl_loc_death), ]
    
    # Reshape dataset long
    long <- data.frame()
    for (i in 1:3) {
      x <- df_desc[which(df_desc[, paste0("list", i)] == 1), ]
      x <- x[, -grep("list", colnames(x))]
      x$list <- paste0("list", i)
      long <- rbind(long, x)
    }
    
    # Add list name
    long <- merge(long, list_names, by = "list", all.x = T)

      
  #...................................      
  ## Tabulate characteristics of each list, and do significance testing
    
    # Initialise table
    tab <- data.frame(
      characteristic = c(
        c("sex", rep(NA, length(table(long$gender))-1)),
        c("age (years)", rep(NA, length(table(long$age_cat))-1)),
        c("location of death", rep(NA, length(table(long$loc_death2))-1)),
        c("year of death", rep(NA, length(table(long$year_death))-1)),
        c("member of resistance committees", 
          rep(NA, length(table(long$resistance_committees))-1)),
        c("cause of death", rep(NA, length(table(long$cod2))-1))
      ),
      category = c(
        names(table(long$gender)),
        names(table(long$age_cat)),
        names(table(long$loc_death2)),
        names(table(long$year_death)),
        names(table(long$resistance_committees)),
        names(table(long$cod2))
      )      
    )
    
    # Add contingency cells
    x <- rbind(
      as.data.frame.matrix(table(long$gender, long$list_name)),
      as.data.frame.matrix(table(long$age_cat, long$list_name)),
      as.data.frame.matrix(table(long$loc_death2, long$list_name)),
      as.data.frame.matrix(table(long$year_death, long$list_name)),
      as.data.frame.matrix(table(long$resistance_committees, long$list_name)),
      as.data.frame.matrix(table(long$cod2, long$list_name))
    )
    tab <- cbind(tab, x)
    
    # Add column-wise percentages
    x <- rbind(
      as.data.frame.matrix(prop.table(table(long$gender, long$list_name), 2)),
      as.data.frame.matrix(prop.table(table(long$age_cat, long$list_name), 2)),
      as.data.frame.matrix(prop.table(table(long$loc_death2,long$list_name),2)),
      as.data.frame.matrix(prop.table(table(long$year_death,long$list_name),2)),
      as.data.frame.matrix(prop.table(table(long$resistance_committees, 
        long$list_name), 2)),
      as.data.frame.matrix(prop.table(table(long$cod2, long$list_name), 2))
    )
    x <- apply(x, c(1,2), percent, accuracy = 0.1)
    colnames(x) <- list_names$list_name
    for (i in list_names$list_name) {
      tab[, i] <- paste0(tab[, i], " (", x[, i], ")")
    }
    
    # Add Fisher p-value test of significance (non-missing categories only)
    f_fisher <- function(x) {
      x <- x[rownames(x) != "missing", ]
      return(fisher.test(x, simulate.p.value = T)$p.value)
    }
    p_value <- c(
      f_fisher(as.data.frame.matrix(table(long$gender, long$list_name))),
      f_fisher(as.data.frame.matrix(table(long$age_cat, long$list_name))),
      f_fisher(as.data.frame.matrix(table(long$loc_death2, long$list_name))),
      f_fisher(as.data.frame.matrix(table(long$year_death, long$list_name))),
      f_fisher(as.data.frame.matrix(table(long$resistance_committees, 
        long$list_name))),
      f_fisher(as.data.frame.matrix(table(long$cod2, long$list_name)))
    )
    p_value <- scales::pvalue(p_value)
    p_value <- data.frame(characteristic = na.omit(tab$characteristic), p_value)
    tab$row <- 1:nrow(tab)
    tab <- merge(tab, p_value, by = "characteristic", all.x = T)
    tab <- tab[order(tab$row), ]
    tab <- subset(tab, select=-row)
    
    # Add list totals
    tab <- rbind(tab, c("total observations", NA, table(long$list_name), NA))
    
    # Save table
    write.csv(tab, paste0(dir_ij, "descriptive.csv"), row.names = F, na = "")   


  #...................................      
  ## Map number of observations per list, by Sudan state
  
    # Set aside the number of missing
    x <- which(long$loc_death == "unknown / unclear")
    missing <- table(long[x, "list_name"])
    
    # Tabulate observations per state, by list
    tab <- as.data.frame.matrix(table(long[-x, c("loc_death", "list_name")]))
    tab$shapeName <- rownames(tab)
    
    # Merge with admin boundaries
    tab <- merge(sdn_boundaries, tab, by = "shapeName", all.x = T)
    for (i in list_names$list_name) {tab[which(is.na(tab[, i])), i] <- 0}
    
    # Create labels
    for (i in list_names$list_name) {
      x <- unlist(sf::st_drop_geometry(tab[, i]))
      tab[, paste0(i, "_lab")] <- paste0(tab$shapeName, "\n (n = ", x, ")")
    }
    
    # Map each list     
    for (i in list_names$list_name) {
      
      # prepare
      x <- unlist(sf::st_drop_geometry(tab[, paste0(i, "_lab")]))
      tab$n_deaths <- unlist(sf::st_drop_geometry(tab[, i]))
      
      # map
      pl <- ggplot(tab) + 
        geom_sf(aes(fill = n_deaths)) + 
        geom_sf_label(aes(label = x), size = 5, alpha = 0.70) +
        labs(x = "latitude", y = "longitude") +
        theme_bw() +
        theme(legend.position = "none", axis.text = element_text(size = 17),
          axis.title = element_text(size = 17)) +
        scale_fill_gradient(
          low = "white",
          high = list_names[which(list_names$list_name == i), "list_colour"],
          na.value = "grey90"
        ) +
        annotate("text", x = 22.5, y = 22.5, size = 5,
          label = paste0("missing: n = ", missing[i]))

      # assign name to plot
      assign(paste0("map_", i), pl)
    }

    # Combine maps
    map_list <- lapply(grep("map", ls(), value = T), get)
    ggarrange(plotlist = map_list, ncol = 2, nrow = 2, 
      labels = list_names$list_name, font.label = list(size = 20), 
      label.x = 0.35, label.y = 0.95)
    ggsave(paste0(dir_ij, "map_lists.png"), units = "cm", dpi = "print", 
      height = 45, width = 60)

  
  #...................................      
  ## Graph timeline of deaths per month, by list
    
    # Prepare data
    long$date_month <- as.Date(paste(15, month(long$date_death), 
      year(long$date_death), sep = "-"), format = "%d-%m-%Y")

    # Set aside the number of missing
    x <- which(is.na(long$date_month))
    missing <- as.data.frame(matrix(NA, nrow = 1, ncol = 3))
    missing[1, ]  <- as.numeric(table(long[x, "list_name"]))  
    colnames(missing) <-list_names$list_name
    
    # Graph 
    ggplot(long[-x, ], aes(x = date_month, group = list_name, colour =list_name, 
      fill = list_name)) +
      geom_bar(position = "dodge", stat = "count", alpha = 0.75) +
      scale_y_continuous("number of deaths listed", expand = c(0, 0)) +
      scale_x_date("date", breaks = "1 month", date_labels = "%b-%Y") +
      theme_bw() +
      scale_color_manual("list:", values = list_names$list_colour) +
      scale_fill_manual("list:", values = list_names$list_colour) +
      theme(legend.position = "top", panel.grid.major.x = element_blank()) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
      annotate("text", x = date_start + 70, y = 190, label = "missing dates:") +
      annotation_custom(
        tableGrob(missing, rows = NULL,
          theme = ttheme_minimal(
            core = list(fg_params = list(cex = 0.8)),
            colhead = list(fg_params = list(cex = 0.8, fontface = "plain")))), 
        xmin = date_start + 180, xmax = date_start + 270, ymin = 190, ymax =180)
      
    ggsave(paste0(dir_ij, "list_by_date.png"), units = "cm", dpi = "print", 
      height = 12, width = 20)


#...............................................................................
### Describing the matched dataset
#...............................................................................
    
  #...................................      
  ## Graph timeline of deaths per month, by cause
    
    # Prepare data
    df_desc$date_month <- as.Date(paste(15, month(df_desc$date_death), 
      year(df_desc$date_death), sep = "-"), format = "%d-%m-%Y")

    # Set aside the number of missing
    x <- which(is.na(df_desc$date_month))
    missing <- length(x)
    
    # Graph 
    ggplot(df_desc[-x, ], aes(x = date_month, group = cod2, alpha = cod2)) +
      geom_bar(position = "stack", stat = "count",
        fill = palette_gen[10], colour = palette_gen[10]) +
      scale_y_continuous("number of unique deaths", expand = c(0, 0)) +
      scale_x_date("date", breaks = "1 month", date_labels = "%b-%Y") +
      theme_bw() +
      scale_alpha_manual("cause:", values = c(0.75, 0.25)) +
      theme(legend.position = "top", 
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank()) +
      annotate("text", x = date_start + 70, y = 300, 
        label = paste0("missing dates: n = ", missing))
    ggsave(paste0(dir_ij, "cause_by_date.png"), units = "cm", dpi = "print", 
      height = 12, width = 20)
     
    
  #...................................      
  ## Graph cause of death by region of Sudan
    
    # Prepare data
    df_desc$loc_death3 <- df_desc$loc_death
    x <- c("Central Darfur", "East Darfur", "North Darfur", "West Darfur",
      "South Darfur")
    df_desc[which(df_desc$loc_death3 %in% x), "loc_death3"] <- "Darfur states"
    x <- c("West Kordofan", "North Kordofan", "South Kordofan")
    df_desc[which(df_desc$loc_death3 %in% x), "loc_death3"] <- "Kordofan states"
    x <- c("Gedaref", "Kassala", "Northern", "Red Sea", "River Nile",
      "Sennar", "White Nile", "Blue Nile")
    df_desc[which(df_desc$loc_death3 %in% x), "loc_death3"] <- "other states"
    tab <- as.data.frame.matrix(
      prop.table(table(df_desc$loc_death3, df_desc$cod, useNA = "ifany"), 1))
    tab$loc <- rownames(tab)
    tab <- reshape(tab, direction = "long", varying = unique(df_desc$cod),
      idvar = "loc", timevar = "cod", times = unique(df_desc$cod), 
      v.names = "prop")
    tab$lab <- scales::percent(tab$prop, accuracy = 1)
    x <- table(df_desc$loc_death3)
    state_labs <- paste0(names(x), "\n (n = ", as.vector(x), ")")
    
    # Plot
    ggplot(tab, aes(x = loc, y = prop, colour = cod, fill = cod)) +
      geom_bar(position = "fill", stat = "identity", alpha = 0.75) +
      geom_text(aes(label = lab), size = 3.5, colour = "black",
        position = position_fill(vjust = 0.5)) +
      theme_bw() +
      theme(legend.position = "top", panel.grid.major.x = element_blank()) +
      scale_x_discrete("state or grouping of states",
        labels = state_labs) +
      scale_y_continuous("percentage of all unique deaths", labels = percent,
        breaks = seq(0, 1, 0.20), expand = c(0, 0)) +
      scale_colour_manual("cause:", values=c(palette_gen[c(5,15,10)],"grey70"))+
      scale_fill_manual("cause:", values=c(palette_gen[c(5,15,10)],"grey70"))
    ggsave(paste0(dir_ij, "cause_by_region.png"), units = "cm", dpi = "print", 
      height = 12, width = 20)
    



#...............................................................................
### ENDS
#...............................................................................

