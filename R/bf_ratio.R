# @ Thomas W. Battaglia

#' Computer Bacteroidetes/Firmicutes ratio for each sample.
#'
#' This function will computer the B/F ratio based on the relative abundances
#' between the Bacteroidetes and Firmicutes phyla.
#'
#' @param phylo An input phyloseq object.
#' @return A phyloseq object with 'log2_bf_ratio' as a new metadata column
#' @export
bf_ratio <- function(phylo){

  # Check to make sure input is a phyloseq object
  if(class(phylo)[1] != "phyloseq"){
    message("Error, input is not a phyloseq object!")
    return(NULL)
  }

  # Collapse on phyla
  phyla <- phyloseq::tax_glom(phylo, "Phylum")

  # Find relative abundances
  phyla_rel <- phyloseq::transform_sample_counts(phyla, function(x) {x/sum(x)} )

  # Keep B/F taxa
  phyla_rel_bact <- suppressWarnings(phyloseq::otu_table(subset_taxa(phyla_rel, Phylum == "Bacteroidetes")))
  phyla_rel_firm <- suppressWarnings(phyloseq::otu_table(subset_taxa(phyla_rel, Phylum == "Firmicutes")))

  # OTU
  bf_ratio <- log2(phyla_rel_bact /  phyla_rel_firm)

  # Add to sample metadata
  sample_data(phylo)$log2_bf_ratio <- as.numeric(bf_ratio)

  # Return phyllseq object
  return(phylo)
}
