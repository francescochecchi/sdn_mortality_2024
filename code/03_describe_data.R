#...............................................................................
### +++++ CAPTURE-RECAPTURE ANALYSIS OF MORTALITY DATA - SUDAN (2024) ++++++ ###
#...............................................................................

#...............................................................................
## ----------- R SCRIPT TO PRODUCE DESCRIPTIVE ANALYSES OF DATASET  --------- ##
#...............................................................................


#...............................................................................
### xxx
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
      list_colour = palette_gen[c(4, 8, 12)])
    
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
    ggarrange(plotlist = map_list, ncol = 3, labels = list_names$list_name,
      font.label = list(face = "plain"))
    ggsave(paste0(dir_ij, "map_lists_wide.png"), units = "cm", dpi = "print", 
      height = 30, width = 90)
    
  









#...............................................................................
### ENDS
#...............................................................................

