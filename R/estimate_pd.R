# @ Thomas W. Battaglia

#' Estimate phylogenetic diversity from phyloseq object
#'
#' Estimate the Faiths phylogenetic diverstiy from an OTU table and phylogenetic tree
#'
#' @param phylo A phyloseq object with an OTU table and phylogenetic tree slot.
#' @return A data.frame of phylogenetic diversity metrics.
#' @export
# Estimate PD-whole tree
estimate_pd <- function(phylo){

  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }

  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }

  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }

  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")

  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)

  # Print status message
  message("Calculating Faiths PD-index...")

  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }

  # Calculate Faith's PD-index
  pdtable <- picante::pd(otutable, tree, include.root = F)

  # Return data frame of results
  return(pdtable)
}
