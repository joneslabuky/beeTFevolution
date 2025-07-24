#!/usr/bin/env Rscript

library(caper)
library(ape)
library(dplyr)
library(tidyr)
library(purrr)
library(phytools)

# Tree
tree <- read.tree("/pscratch/bmjo263_uksr/SAV_TFBS/Results_Nov01/Species_Tree/SpeciesTree_rooted_node_labels.txt")
tree$tip.label <- gsub("_longpeps", "", tree$tip.label)

# Define species and sociality info
info<-data.frame(species = c("AAUR", "ACER", "ADOR", "AFUL", "AMEL", "APLU", "APUR", "AVIR", "BAFF", "BTER", "CCAL", "CGIG", "DNOV", "EMEX", "FVAR", "HANT", "HLAB", "HLIG", "HQUA",
                             "HVOL", "LFIG", "LLAT", "LLEU", "LMAR", "LMOR", "LOEN", "LPAU", "LVIE", "LZEP", "MBIC", "MEUR", "MROT", "NFAB", "NFLA",  "NMEL", "OBIC", "OLIG", "SMON",
                             "SPHA", "STUM", "TDIV", "XVIO"),
                 sociality = as.factor(c("soc", "soc", "sol", "sol", "soc", "sol", "sol_loss", "sol", "soc", "soc", "sol", "sol", "sol", "sol", "soc", "sol", "sol", "soc", "sol_loss",
                                         "sol", "sol_loss", "sol_loss", "sol_loss", "soc", "soc", "sol_loss", "soc", "sol_loss", "soc", "soc", "sol", "sol", "sol", "sol", "sol",
                                         "sol", "sol", "sol","sol", "soc", "sol", "sol")),
                 group = as.factor(c("old","new","old","old","old","new","old","old","new","old","new","new","old","new","new","new","new",
                                    "old","old","new","old","new","new","old","new","old","old","old","old","new","new","old","old","new",
                                    "old","new","old","new","new","new","new","new")),
                 family = as.factor(c("halictid","apid","andrenid","andrenid","apid","apid","halictid","halictid","apid","apid","apid",
                                      "colletid","halictid","apid","apid","colletid","apid","halictid","halictid",
                                      "colletid","halictid","halictid","halictid","halictid","halictid","halictid","halictid","halictid","halictid",
                                      "apid","melittid","megachilid","apid","apid","halictid","megachilid","megachilid","halictid",
                                      "megachilid","halictid","apid","apid")))



ult_tree<-force.ultrametric(tree)

# MAD, MeanAD

handle_outliers <- function(y) {
  x <- y[!is.na(y)]
  med <- median(x)
  mad_val <- median(abs(x - med))
  mean_ad_val <-mean(abs(x- med))

  if(mad_val != 0) {
     threshold <- 2.5 * mad_val
  } else {
     threshold <- 2.5 * mean_ad_val
  }

  x[x > med + threshold] <- med + threshold
  x[x < med - threshold] <- med - threshold
  y[!is.na(y)] <- x
  return(y)
}

# min-max normalizing data
normalize<-function(y) {
    x<-y[!is.na(y)]
    x<-(x - min(x)) / (max(x) - min(x))
    y[!is.na(y)]<-x
    return(y)
}

# Define species for each sociality group and for splititng tree by family 
soc_species<-info[which(info$sociality=="soc"),"species"]
sol_loss_species<-info[which(info$sociality=="sol_loss"),"species"]
sol_species<-info[which(info$sociality=="sol"),"species"]

soc_apid<-info[which(info$sociality=="soc" & info$family=="apid"),"species"]
sol_apid<-info[which(info$sociality=="sol" & info$family=="apid"),"species"]

soc_hal<-info[which(info$sociality=="soc" & info$family=="halictid"),"species"]
sol_hal<-info[which(info$sociality=="sol" & info$family=="halictid"),"species"]
solloss_hal<-info[which(info$sociality=="sol_loss" & info$family=="halictid"),"species"]

# Function for running PGLS with different bounds
run_pgls_with_bound <- function(data, bounds, ult_tree, info) {
  result <- NULL
  for (i in 1:nrow(data)) {
    OG <- rownames(data)[i]
    scores <- data[i,]
    scores_df <- data.frame(species = colnames(data), score = as.numeric(scores))
    merged_data <- merge(scores_df, info, by = "species", all.x = TRUE)
    merged_data$sociality <- droplevels(merged_data$sociality)
    comparative_data <- comparative.data(phy = ult_tree, data = merged_data, species, vcv = TRUE)
    
    # Loop through each bound
    for (bound in bounds) {
      pgls_result <- tryCatch({
        pgls(score ~ sociality, data = comparative_data, lambda = 'ML', bounds = list(lambda = bound))
      }, error = function(e) return(NULL))  # Return NULL if PGLS fails with this bound

      # If PGLS is successful, calculate R^2 and adjusted R^2, store all results and exit the loop
      if (!is.null(pgls_result)) {
        r_squared <- summary(pgls_result)$r.squared
        n <- nrow(merged_data) # Number of observations (species used in the PGLS)
        p <- 1 # Number of predictors 
        r_squared_adj <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))
  
        return(data.frame(Regression_Coefficient = summary(pgls_result)$coefficients[, 1][2],
                           P_value = summary(pgls_result)$coefficients[, 4][2],
                           Lambda = summary(pgls_result)$param[2],
                           Bound_used = paste(bound[1], "to", bound[2]),
                           R_squared = r_squared, 
                           R_squared_adj = r_squared_adj,
                           n = n,
                           OG = OG, Status = "success"))
      }
    }
    
    # If no bound worked, return failure result
    return(data.frame(OG = OG, Regression_Coefficient = NA, P_value = NA, Lambda = NA, Bound_used = "None", R_squared = NA, R_squared_adj = NA, n = NA,  Status = "fail"))
  }
}

# Function to process motif file and apply bounds
process_motif_file <- function(file, tree, info) {
  # Read the data
  data <- read.table(file, header = TRUE, sep = "\t")
  data <- data[-1,]
  rownames(data) <- data$index
  data$index <- NULL
  
  # APPLY MAD MeanAD OUTLIER adjustment
  data [] <- lapply(data, function(col) if (is.numeric(col)) handle_outliers(col) else col)

  # USE NORMALIZED DATA
  data [] <- lapply(data, function(col) if (is.numeric(col)) normalize(col) else col)
  
  data <- as.data.frame(data)

  # Filter data with too many NAs
  soc_count <- rowSums(!is.na(data[soc_species]))
  anc_sol_count <- rowSums(!is.na(data[sol_species]))
  sol_loss_count <- rowSums(!is.na(data[sol_loss_species]))

  
  # Data for PGLS when running on entire 42 species tree
  data_PGLS1 <- data[which(soc_count >= 4 & anc_sol_count >= 7), c(soc_species, sol_species)]
  data_PGLS2 <- data[which(soc_count >= 4 & sol_loss_count >= 2), c(soc_species, sol_loss_species)]
   
  # Remove rows with zero variance
   data_PGLS1 <- data_PGLS1[apply(data_PGLS1, 1, function(x) var(na.omit(x)) >0), ]
   data_PGLS2 <- data_PGLS2[apply(data_PGLS2, 1, function(x) var(na.omit(x)) >0), ]  
                        
  # Define bounds scenarios
  bounds_scenarios <- list(c(1e-05, 1), c(1e-03, 1), c(1e-01, 1), c(1, 1))

  # Initialize results data frames
  results_PGLS1 <- map_dfr(1:nrow(data_PGLS1), function(i) {
    result <- NULL
    # Extract specific OG to process
    single_og_data <- data_PGLS1[i, , drop = FALSE]
    result <- run_pgls_with_bound(single_og_data, bounds_scenarios, ult_tree, info)
    
    # Return result if found, else return a failure
    if (!is.null(result)) {
      result$OG <- rownames(data_PGLS1)[i]
      return(result)
    }
    return(data.frame(OG = rownames(data_PGLS1)[i], Correlation = NA, P_value = NA, Lambda = NA, Bound_used = "None", Status = "fail"))
  })

  results_PGLS2 <- map_dfr(1:nrow(data_PGLS2), function(i) {
    result <- NULL
    # Extract specific OG to process
    single_og_data <- data_PGLS2[i, , drop = FALSE]
    result <- run_pgls_with_bound(single_og_data, bounds_scenarios, ult_tree, info)
    
    # Return result if found, else return a failure
    if (!is.null(result)) {
      result$OG <- rownames(data_PGLS2)[i]
      return(result)
    }
    return(data.frame(OG = rownames(data_PGLS2)[i], Correlation = NA, P_value = NA, Lambda = NA, Bound_used = "None", Status = "fail"))
  })

  # Merge results from PGLS1 and PGLS2
  results_merge <- merge(results_PGLS1, results_PGLS2, by = "OG", all.x = TRUE, all.y = TRUE)
  return(results_merge)
  # return(results_PGLS1)
}

# Run the function on the input file
results <- process_motif_file(snakemake@input[[1]], ult_tree, info)

# Save results to a separate table
write.csv(results, snakemake@output[[1]], row.names = FALSE)

