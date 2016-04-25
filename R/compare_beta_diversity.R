# @ Thomas W. Battaglia

#' Calculate alpha diversity statistics
#'
#' Compute p-values and multiple comparisons adjusted q-values for
#' two-group comparisons across multiple timepoints.
#'
#' @param phylo A phyloseq object
#' @param x The variable that describes Time
#' @param group The variable that describes 2 or more groups
#' @param test The PERMANOVA test to use.
#' @param bdiv Beta diversity metric to calculate significance.
#' @param write Should results be written to a file in the current working directory.
#' @param filename The name of the output file if write is True.
#' @return A dataframe for an PERMANOVA test over each timepoint from each two group comparison.
#' @export
compare_beta_diversity <- function(phylo,
                                   x = "Day",
                                   group = "Treatment",
                                   test = c("adonis", "anosim"),
                                   bdiv = c("unweighted", "weighted"),
                                   write = F,
                                   filename = "results", ...){

  # Sample metadata
  metadata <- as(phyloseq::sample_data(phylo), "data.frame")

  # Get the different levels (time) for comparing
  x_levels <- levels(metadata[ ,x])

  # Get comparing groups
  comparing_groups <- levels(metadata[ ,group])

  # Make dataframes to store results
  if(test == "adonis"){
    results_df <- data.frame(Group1 = as.character(),
                             Group2 = as.character(),
                             x = as.character(),
                             n = as.numeric(),
                             SumsOfSqs = as.numeric(),
                             MeanSqs = as.numeric(),
                             F.Model = as.numeric(),
                             R2 = as.numeric(),
                             pvalue = as.numeric())
  }

  if(test == "anosim"){
    results_df <- data.frame(Group1 = as.character(),
                             Group2 = as.character(),
                             x = as.character(),
                             n = as.numeric(),
                             R_value = as.numeric(),
                             pvalue = as.numeric())
  }


  # Create list of combinations
  comb_list <- combn(comparing_groups, 2, simplify = F)

  # Loop over each comparison
  for(comp in comb_list){

    # Subset metadata
    metadata_comparison_subset <- metadata[metadata[ ,group] %in% comp,]
    metadata_comparison_subset[ ,group] <- droplevels(metadata_comparison_subset[ ,group])

    # Subset taxadata
    taxadata_comparison_subset <- phylo
    sample_data(taxadata_comparison_subset) <- phyloseq::sample_data(metadata_comparison_subset)


    # Loop over each level (x)
    for(i in x_levels){

      # Subset metadata by x
      metadata_comparison_subset_time <- subset(metadata_comparison_subset, metadata_comparison_subset[ ,x] == i)

      # Subset taxonomy by x
      taxadata_comparison_subset_time <- taxadata_comparison_subset
      sample_data(taxadata_comparison_subset_time) <- phyloseq::sample_data(metadata_comparison_subset_time)



      ### Get Unifrac Distances ------------
      # Calculate unweighted values
      if(bdiv == "unweighted"){
        unifrac <- phyloseq::UniFrac(taxadata_comparison_subset_time, weighted = FALSE, normalized = TRUE, parallel = FALSE)
      }

      # Calculate weighted values
      if(bdiv == "weighted"){
        unifrac <- phyloseq::UniFrac(taxadata_comparison_subset_time, weighted = TRUE, normalized = TRUE, parallel = FALSE)
      }



      ### Get Stastical Values ------------
      # perform adonis test
      if(test == "adonis"){

        # Perform Adonis test
        results <- suppressWarnings(try(vegan::adonis(formula = as.dist(unifrac) ~ metadata_comparison_subset_time[ ,group], permutations = 999)))

        # If too little amount of samples are present for either group, result in None.
        if(class(results) == "try-error"){
          pval <- 'NA'
          SumsOfSqs <- "NA"
          MeanSqs <- "NA"
          F.Model <- "NA"
          R2 <- "NA"
        }

        # If no error, assign results to variables
        if(class(results) == "adonis") {
          pval <- results$aov.tab$`Pr(>F)`[1]
          SumsOfSqs <- results$aov.tab$SumsOfSqs[1]
          MeanSqs <- results$aov.tab$MeanSqs[1]
          F.Model <- results$aov.tab$F.Model[1]
          R2 <- results$aov.tab$R2[1]
        }

        # Find which groups are being compared.
        compared <- levels(metadata_comparison_subset_time[ ,group])

        # Place results into dataframe
        results_mat <- data.frame(Group1 = compared[1],
                                  Group2 = compared[2],
                                  x = i,
                                  n = nrow(metadata_comparison_subset_time),
                                  SumsOfSqs = SumsOfSqs,
                                  MeanSqs = MeanSqs,
                                  F.Model = F.Model,
                                  R2 = R2,
                                  pvalue = pval)

      } # EOF Adonis

      # perform anosim test
      if(test == "anosim"){

        # test
        results <- suppressWarnings(try(vegan::anosim(as.dist(unifrac), metadata_comparison_subset_time[,group], permutations = 999)))

        # If too little amount of samples are present for either group, result in None.
        if(class(results) == "try-error"){
          pval <- 'NA'
          R_value <- "NA"
        }

        # If no error, assign results to variables
        if(class(results) == "anosim"){
          pval <- results$signif
          R_val <- results$statistic
        }

        # Find which groups are being compared.
        compared <- levels(metadata_comparison_subset_time[ ,group])

        # Place results into dataframe
        results_mat <- data.frame(Group1 = compared[1],
                                  Group2 = compared[2],
                                  x = i,
                                  n = nrow(metadata_comparison_subset_time),
                                  R_value = R_val,
                                  pvalue = pval)
      } # EOF Anosim

      # Merge results
      results_df <- rbind(results_mat, results_df)

    } # EOF levels

  } # EOF comparisons

  # Write table if chosen
  if (write == T){
    write.csv(x = results_df, file = paste(filename, ".csv", sep = ""))
  }

  #if (write_plots == T){
    #return(list_of_plots)
  #}

  return(results_df)
}

