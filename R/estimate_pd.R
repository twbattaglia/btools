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

  # Get OTU table data matrix
  otutable <- t(as(otu_table(phylo, taxa_are_rows = F), "matrix"))

  # Get phylogenetic tree from pyloseq object
  tree <- phy_tree(phylo)

  # Calculate Faiths PD
  message("Calculating Faiths PD-index...")
  pdtable <- picante::pd(otutable, tree, include.root = F)

  # Return data frame of results
  return(pdtable)
}
