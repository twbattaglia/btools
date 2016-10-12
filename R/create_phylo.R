# @ Thomas W. Battaglia

#' Generic function for creating a phyloseq object from QIIME intput data.
#'
#' This function loads an OTU table (.biom format) a phylogenetic tree file (.tre) and
#' a sample mapping file (.txt/.tsv) to create a standard phyloseq object. Input BIOM file
#' can be in either JSON format or HDF5.
#'
#' @param biom_fp File location of OTU table in .biom format.
#' @param mappingfile_fp File location of the input sample metadata matching the input OTU table.
#' @param tree_fp File location of the input sample metadata matching the input OTU table.
#' @return A phyloseq object to be used downstream for many different analyses.
#' @export
create_phylo <- function(biom_fp, mappingfile_fp, tree_fp){

  # Import biom
  x = suppressMessages(biomformat::read_biom(biom_fp))
  otutab = otu_table(as(biom_data(x), "matrix"), taxa_are_rows = TRUE)

  # Import taxa
  taxlist = lapply(x$rows, function(i) {parse_taxonomy_default(i$metadata$taxonomy)})
  names(taxlist) = sapply(x$rows, function(i) {i$id})
  taxtab = build_tax_table(taxlist)

  # Import tree
  tree <- phyloseq::read_tree(tree_fp)

  # Import map
  mapping =  phyloseq::import_qiime_sample_data(mapfilename = mappingfile_fp)

  # Create phylodata object
  phylodata <- phyloseq(otutab, taxtab, mapping, tree)

  # Change rank names
  if(phyloseq::rank_names(phylodata) != c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")){
    # Fix rank names
    colnames(tax_table(phylodata)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
  }

  # Return object for phyloseq
  message("Phyloseq object sucessful created")
  return(phylodata)
}
