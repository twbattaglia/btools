# @ Thomas W. Battaglia

#' Remove blanks from OTU table
#'
#' This function will remove the OTU from blanks within the
#' an OTU table, for each plate/run
#'
#' @param phylo A phyloseq object with all samples per run.
#' @param runID The column name used to identify the different runs/plates.
#' @param blankID The column name used to identify the blanks.
#' @param blankName The name of the factor which describes the blanks.
#' @param removeBlank Should blank sample be removed from resulting phyloseq object.
#' @return A phyloseq object with the OTU's removed from each plate per plate-blank
#' @export
remove_blanks <- function(phylo,
                          runID = "Run_number",
                          blankID = "Group",
                          blankName = "blank",
                          removeBlank = F){

  # Make a copy of the main OTU table
  phylo_noblanks = phylo

  # Make a data.frame of metadata
  metadata <- as(sample_data(phylo_noblanks), "data.frame")

  # Split mapping file by runID
  split_PCR <- split(metadata, metadata[[runID]])

  # Run function over each PCR run
  for(map in split_PCR){

    # Create a temporary phyloseq object to store subsetting
    phylo_temp = phylo_noblanks

    # Drop unused levels
    map[[runID]] <- droplevels(map[[runID]])

    # Message
    message(paste0("Working on run: ", unique(map[[runID]])), appendLF = T)

    # Subset map to only include blank samples
    blank_map <- map[which(map[[blankID]] == blankName),]

    # Change metadata to reflect blank sample metadata
    sample_data(phylo_temp) <- blank_map

    # Remove zero sum OTU's
    blanks = prune_taxa(taxa_sums(phylo_temp) > 0, phylo_temp)

    # Print message
    message(paste0("Number of taxa found: ", ntaxa(blanks)), appendLF = T)

    # Get a list of the OTU's founds in the blanks
    blank_otu <- taxa_names(blanks)
    blank_otu_tofilter <- !(row.names(tax_table(phylo_noblanks)) %in% blank_otu)

    # Prune OTU's from original phylo object
    phylo_noblanks <- prune_taxa(blank_otu_tofilter, phylo_noblanks)
  }

  # Remove control samples
  if(removeBlank){
    phylo_noblanks <- phyloseq::subset_samples(phylo_noblanks, metadata[[blankID]] != blankName)
  }

  # Message before and after stats
  message(paste0("Number of OTU's before filtering: ", ntaxa(phylo)))
  message(paste0("Number of OTU's after filtering: ", ntaxa(phylo_noblanks)))
  message(paste0("Difference: ", (ntaxa(phylo) - ntaxa(phylo_noblanks))))

  # Retun filtered dataframe
  return(phylo_noblanks)
}
