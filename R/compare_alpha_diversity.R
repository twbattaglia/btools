# @ Thomas W. Battaglia
# @ tb1280@nyu.edu

compare_alpha_diversity <- function(physeq, 
                                    x = "Day", 
                                    group = "Treatment", 
                                    diversity = "Observed", 
                                    test_type = "nonparametric", 
                                    num_perm = 999, 
                                    multiple_corrections = F,
                                    write = F,
                                    filename = "results", ...){

  # Load libraries
  if (!require("phyloseq")) {
    install.packages("phyloseq", repos="http://cran.rstudio.com/") 
  }
  suppressPackageStartupMessages(library(phyloseq))
  
  if(class(physeq)=="phyloseq"){
    # Print message
    message("Recognized input as phyloseq object. Proceeding with analysis.")
    
    # Retrieve data from phyloseq object
    phylo_data <- plot_richness(physeq = physeq, x = x, measures = diversity, color = group)$data
  }
  
  if(class(physeq)=="data.frame"){
    # Print message
    message("Recognized input as data-frame. Proceeding with analysis.")
    
    phylo_data = physeq
    names(phylo_data)[which(names(phylo_data)=="PD_whole_tree_alpha")] <- "value"
  }
  
  # Get xx levels
  x_levels <- levels(phylo_data[ ,x])
  
  # Get comparing groups
  comparing_groups <- levels(phylo_data[ ,group])
  
  # Make dataframe to store results
  results <- data.frame(Group1 = as.character(),
                        Group2 = as.character(),
                        x = as.character(), 
                        n = as.numeric(), 
                        Group1_mean = as.numeric(),
                        Group1_std = as.numeric(),
                        Group2_mean = as.numeric(),
                        Group2_std = as.numeric(),
                        t_stat = as.numeric(),
                        pvalue = as.numeric())
  
  # Create list of combinations
  comb_list <- combn(comparing_groups, 2, simplify = F)
  
  # For each combination of factors
  for(comp in comb_list){
    
    # Subset table for only comparisons on interest.   
    phylo_data_comsub <- phylo_data[phylo_data[ ,group] %in% comp,]
    phylo_data_comsub[ ,group] <- droplevels(phylo_data_comsub[ ,group])
    
    # Iterate over each x level. (Run through each timepoint)
    for(i in x_levels){
      
      # Subset table
      data_table <- subset(phylo_data_comsub, phylo_data_comsub[ ,x] == i)
      
      # Parametric T-test
      if(test_type == "parametric"){
        
        # Perform parametric t.test
        sig_res <- try(t.test(data_table$value ~ data_table[ ,group]), silent = TRUE)
        
        # If too little amount of samples are present for either group, result in None.
        if(class(sig_res) == "try-error"){
          tstat <- "NA"
          pval <- 'NA'
        }
        else{
          # If no error, assign results to variables
          tstat <- sig_res[["statistic"]][[1]]
          pval <- sig_res[["p.value"]][[1]]
        }
        
      }
      
      
      # Non Parametric Monte Carlo Bootstrap
      if(test_type == "nonparametric"){
        
        # Load nonparametric testing library
        if (!require("coin")) {
          install.packages("coin", repos="http://cran.rstudio.com/") 
        }
        suppressPackageStartupMessages(library(coin))
        
        # Perform non-parametric t.test monte carlo simulations
        sig_res <- try(oneway_test(data_table$value ~ data_table[ ,group], distribution = approximate(B = num_perm)), silent = TRUE)
        
        # If too little amount of samples are present for either group, result in None.
        if(class(sig_res) == "try-error"){
          tstat <- "NA"
          pval <- 'NA'
        }
        else{
          # If no error, assign results to variables
          tstat <- sig_res@statistic@teststatistic[[1]]
          pval <- pvalue(sig_res)[[1]]
        }
      }
      
      # Find which groups are being compared.
      compared <- levels(data_table[ ,group])
      
      # Find mean and standard deviations
      means_out <- by(data = data_table$value, data_table[ ,group], mean)
      std_out <- by(data = data_table$value, data_table[ ,group], sd)
      
      # Place results into dataframe
      results_mat <- data.frame(Group1 = compared[1], 
                                Group2 = compared[2], 
                                x = i, 
                                n = nrow(data_table),
                                Group1_mean = means_out[1],
                                Group1_std = std_out[1],
                                Group2_mean = means_out[2],
                                Group2_std = std_out[2],
                                t_stat = tstat,
                                pvalue = pval,
                                row.names = NULL)
      
      # Merge results with existing results
      results <- rbind(results_mat, results)
    
    } # EOF Loop through each value of x
    
  } # EOF Comparing different groups

  
  # Get multiple comparisons qvalue
  if (multiple_corrections == T){
    results$qvalue <- p.adjust(results$pvalue, method = "fdr")
  }
  
  # Write table if chosen
  if (write == T){
    write.csv(x = results, file = paste(filename, ".csv", sep = ""))
  }
  
  # Return Table of results
  return(results)
}



