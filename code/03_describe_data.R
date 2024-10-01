#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - SUDAN (2024) ++++++ ###
#...............................................................................

#...............................................................................
## ----------- R SCRIPT TO PRODUCE DESCRIPTIVE ANALYSES OF DATASET  --------- ##
#...............................................................................

# for each duplication and overlap score combination...
for (dd in 1:5) {
  print(paste0("now doing descriptive analysis for duplication threshold ", dd))  
  
  for (oo in 1:5) {
    print(paste0("   ...and overlap threshold ", oo))  
  
    
#...............................................................................
### Preparing the dataset and tabulating attrition
#...............................................................................

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
    write.csv(tab, paste0(dir_do, "attrition.csv"), row.names = F)    

    
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
    write.csv(tab, paste0(dir_do, "descriptive.csv"), row.names = F, na = "")   


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
    
    # Map each list     
    for (i in list_names$list_name) {
      
      # prepare
      tab$n_deaths <- unlist(sf::st_drop_geometry(tab[, i]))
      tab$n_deaths_cat <- cut(tab$n_deaths, breaks = c(0,1,20,100,500,100000),
        include.lowest = T, right = F, labels = c("0",
          "1 to 19", "20 to 99", "100 to 499", ">=500"))
      
      # create labels
      tab$labs <- paste0(tab$shapeName, "\n (n = ", tab$n_deaths, ")")
      
      # map
      pl <- ggplot(tab) + 
        geom_sf(aes(alpha = n_deaths_cat), 
          fill = list_names[which(list_names$list_name == i), "list_colour"]) + 
        geom_label_repel(aes(label = labs, geometry = geometry), size = 5, 
          alpha = 0.75, stat = "sf_coordinates") +
        theme_bw() +
        theme(legend.position = "none", axis.text = element_text(size = 17),
          axis.title = element_text(size = 17)) +
        scale_alpha_manual("number of deaths", values=c(0,0.25,0.50,0.75,1)) +
        annotate("text", x = 22.5, y = 22.5, size = 5,
          label = paste0("missing: n = ", missing[i])) +
        labs(y = "latitude", x = "longitude")

      # assign name to plot
      assign(paste0("map_", i), pl)
    }

    # Combine maps
    suppressWarnings(rm(map_list))
    map_list <- lapply(grep("map", ls(), value = T), get)
    ggarrange(plotlist = map_list, ncol = 2, nrow = 2, 
      labels = list_names$list_name, font.label = list(size = 20), 
      label.x = 0.35, label.y = 0.95)
    ggsave(paste0(dir_do, "map_lists.png"), units = "cm", dpi = "print", 
      height = 45, width = 60)

  
  #...................................      
  ## Graph timeline of deaths per month, by list
    
    # Prepare data
    long$date_month <- as.Date(paste("01", month(long$date_death), 
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
      scale_y_continuous("number of deaths listed", 
        expand = expansion(add = c(0, 20))) +
      scale_x_date("date", breaks = "1 month", date_labels = "%b-%Y", 
        expand = c(0,0)) +
      theme_bw() +
      scale_color_manual("list:", values = list_names$list_colour) +
      scale_fill_manual("list:", values = list_names$list_colour) +
      theme(legend.position = "top", panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.7)) +
      annotate("text", x = date_start + 100, y = 200, label = "missing dates:",
        size = 3) +
      annotation_custom(
        tableGrob(missing, rows = NULL,
          theme = ttheme_minimal(
            core = list(fg_params = list(cex = 0.7)),
            colhead = list(fg_params = list(cex = 0.7, fontface = "plain")))), 
        xmin = date_start + 180, xmax = date_start + 270, ymin = 200, ymax =190)
      
    ggsave(paste0(dir_do, "list_by_date.png"), units = "cm", dpi = "print", 
      height = 12, width = 20)


  #...................................      
  ## Graph sex proportion of deaths per month, by list
    
    # Prepare data
    long$date_month <- as.Date(paste("01", month(long$date_death), 
      year(long$date_death), sep = "-"), format = "%d-%m-%Y")
    
    # Set aside the number of missing
    x <- which(is.na(long$date_month) | long$gender == "missing")

    # Graph 
    ggplot(long[-x, ], aes(x = date_month, group = gender, colour = list_name, 
      fill = list_name, alpha = gender)) +
      geom_bar(position = "fill", stat = "count") +
      geom_text(aes(label = after_stat(count)), position = "fill", 
        stat = "count", colour = "black", size = 3, vjust = 1.3, alpha = 1) +
      scale_y_continuous("proportion of deaths listed", 
        expand = expansion(add = c(0, 0)), labels = percent) +
      scale_x_date("date", breaks = "1 month", date_labels = "%b-%Y",
        expand = c(0,0)) +
      theme_bw() +
      scale_color_manual("list:", values = list_names$list_colour) +
      scale_fill_manual("list:", values = list_names$list_colour) +
      scale_alpha_manual("sex", values = c(0.75, 0.25, 1)) +
      facet_grid(list_name ~ .) +
      theme(legend.position = "top", panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.7))

    ggsave(paste0(dir_do, "list_by_date_sex.png"), units = "cm", dpi = "print", 
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
      scale_y_continuous("number of unique deaths", 
        expand = expansion(add = c(0, 20))) +
      scale_x_date("date", breaks = "1 month", date_labels = "%b-%Y") +
      theme_bw() +
      scale_alpha_manual("cause:", values = c(0.75, 0.25)) +
      theme(legend.position = "top", 
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.7),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank()) +
      annotate("text", x = date_start + 70, y = 300, 
        label = paste0("missing dates: n = ", missing), size = 3)
    ggsave(paste0(dir_do, "cause_by_date.png"), units = "cm", dpi = "print", 
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
      theme(legend.position = "top", panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
      scale_x_discrete("state or grouping of states",
        labels = state_labs) +
      scale_y_continuous("percentage of all unique deaths", labels = percent,
        breaks = seq(0, 1, 0.20), expand = c(0, 0)) +
      scale_colour_manual("cause:", values=c(palette_gen[c(5,15,10)],"grey70"))+
      scale_fill_manual("cause:", values=c(palette_gen[c(5,15,10)],"grey70"))
    ggsave(paste0(dir_do, "cause_by_region.png"), units = "cm", dpi = "print", 
      height = 12, width = 20)
    

  #...................................      
  ## Plot overlap among lists
    
    # Select dataset
    df_all <- df[which(!df$excl_del & !df$excl_date & !df$excl_loc_death &
        !df$excl_kht), ]
    
    # Rename list columns and set to logical
    df_all[, list_names$list_name] <- df_all[, paste0("list", 1:3)]
    df_all[, list_names$list_name] <- apply(df_all[, list_names$list_name], 2,
      as.logical)
    
    # Plot all-cause deaths
    pl_all <- ggvenn(df_all[, list_names$list_name], 
      fill_color = list_names$list_colour, stroke_alpha = 0.5, fill_alpha = 0.3,
      stroke_size = 0.5, set_name_size = 4) + 
      theme_void() +
      theme(plot.background = element_rect(fill = "white", colour = NA))
    ggsave(paste0(dir_do, "venn_all_causes.png"), units = "cm", dpi = "print", 
      height = 12, width = 20)

    # Plot intentional injury deaths
    pl_inj <- ggvenn(df_all[which(!df_all$excl_cod), list_names$list_name], 
      fill_color = list_names$list_colour, stroke_alpha = 0.5, fill_alpha = 0.7,
      stroke_size = 0.5, set_name_size = 4) + 
      theme_void() +
      theme(plot.background = element_rect(fill = "white", colour = NA))
    ggsave(paste0(dir_do, "venn_intl_injury.png"), units = "cm", dpi = "print", 
      height = 12, width = 20)
    
    # Combined plot
    ggarrange(plotlist = list(NULL, NULL, pl_all, pl_inj), 
      ncol = 2, nrow = 2, 
      labels = c("all causes", "intentional injury", "", ""), 
      font.label = list(face = "plain"), heights = c(1, 50)) +
      theme_void() +
      theme(plot.background = element_rect(fill = "white", colour = NA))
    ggsave(paste0(dir_do, "venn_combi.png"), units = "cm", dpi = "print", 
      height = 15, width = 22)
    
}
}
# closing duplication score and overlap score threshold loops

#...............................................................................
### ENDS
#...............................................................................

