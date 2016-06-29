# @ Thomas W. Battaglia

#' Generic function for creating a phyloseq object from QIIME intput data.
#'
#' This function loads an TU table (.biom format) a phylogenetic tree file (.tre) and
#' a sample mapping file (.txt/.tsv) to create a standard phyloseq object. It is important
#' to first convert the OTU table to JSON format, due to an issue currently with the biom package.
#'
#' @param biom_fp File location of OTU table in .biom format.
#' @param mappingfile_fp File location of the input sample metadata matching the input OTU table.
#' @param tree_fp File location of the input sample metadata matching the input OTU table.
#' @return A phyloseq object to be used downstream for many different analyses.
#' @export
create_phylo <- function(biom_fp, mappingfile_fp, tree_fp){

  # If input biom is a text file
  #if(class(biom_fp)=="data.frame"){
  #}

  # Import BIOM + tree
  data = phyloseq::import_biom(BIOMfilename = biom_fp,
                               treefilename = tree_fp,
                               parseFunction = parse_taxonomy_greengenes)

  # Import mapping file
  mapping =  phyloseq::import_qiime_sample_data(mapfilename = mappingfile_fp)

  # Create phyloseq object
  phylo <- phyloseq::merge_phyloseq(data, mapping)

  # Change rank names
  if(phyloseq::rank_names(phylo) != c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")){
    message("Rank name error. Please correct before using!")
  }

  # return object for phyloseq
  message("Phyloseq object sucessful created")
  return(phylo)
}
