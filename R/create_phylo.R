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

  # Import BIOM + tree
  data = phyloseq::import_biom(BIOMfilename = biom_fp,
                               treefilename = tree_fp,
                               parseFunction = parse_taxonomy_greengenes)

  # Import mapping file
  mapping =  phyloseq::import_qiime_sample_data(mapfilename = mappingfile_fp)

  # Create phyloseq object
  phylo <- phyloseq::merge_phyloseq(data, mapping)

  # Change rank names
  suppressMessages(colnames(tax_table(phylo)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

  # return object for phyloseq
  message("Phyloseq object sucessful created")
  return(phylo)
}
