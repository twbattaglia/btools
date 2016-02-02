# @ Thomas W. Battaglia

#' Add meta data to contributions file
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param metagenomic_contributions.py output file
#' @param mapping file with sample metadata
#' @param list of column variables to add to new table
#' @return A dataframe with taxa information and sample metadata
#' @export
metagenomic_contributions <- function(input_table = table_input,
                                  mappingfile = mappingfile_input,
                                  metadata_col = metadata_col_input){

  # - - - - - - - - - - - - -
  # Error handling
  # - - - - - - - - - - - - -
  # input is a string
  # input is a data.frame
  # mapping file is a string
  # mapping file is a dataframe
  # column id's are not strings


  # - - - - - - - - - - - - -
  # Import files
  # - - - - - - - - - - - - -
  # Import Neccessary Datatables from given input locations.
  input <- read.delim(input_table, header = TRUE)
  gg_db <- otu_taxonomy_97
  mappingfile <- read.delim(mappingfile, header = TRUE)
  input.df <- dplyr::tbl_df(input)


  # - - - - - - - - - - - - -
  # Process strings
  # - - - - - - - - - - - - -
  message("Converting OTU-IDs's to GreenGenes Taxanomy...")
  otuid_name <- gg_db[match(input.df$OTU, table = gg_db$V1), ]$V2
  otuid_name <- gsub(pattern = ';', replacement = "", x = otuid_name, fixed = T)
  otuid_name_sep <- sapply(otuid_name, function(x) stringr::str_split(x, " "))


  # - - - - - - - - - - - - -
  # Adding Names to Table
  # - - - - - - - - - - - - -
  message("Adding Names to Table..")
  input.df$kingdom <- gsub(pattern = 'k__', replacement = "", x = as.character(lapply(otuid_name_sep, function(x) x[[1]])), fixed = T)
  input.df$phylum <- gsub(pattern = 'p__', replacement = "", x = as.character(lapply(otuid_name_sep, function(x) x[[2]])), fixed = T)
  input.df$class <- gsub(pattern = 'c__', replacement = "", x = as.character(lapply(otuid_name_sep, function(x) x[[3]])), fixed = T)
  input.df$order <- gsub(pattern = 'o__', replacement = "", x = as.character(lapply(otuid_name_sep, function(x) x[[4]])), fixed = T)
  input.df$family <- gsub(pattern = 'f__', replacement = "", x = as.character(lapply(otuid_name_sep, function(x) x[[5]])), fixed = T)
  input.df$genus <- gsub(pattern = 'g__', replacement = "", x = as.character(lapply(otuid_name_sep, function(x) x[[6]])), fixed = T)
  input.df$species <- gsub(pattern = 's__', replacement = "", x = as.character(lapply(otuid_name_sep, function(x) x[[7]])), fixed = T)

  # - - - - - - - - - - - - -
  # Add Metadata infomation
  # - - - - - - - - - - - - -
  message("Collecting Metadata Information and Adding to Table..")
  metadata_names<- unlist(stringr::str_split(metadata_col, ","))
  metadata_tbl <- data.frame(mappingfile[match(input.df$Sample, mappingfile$X.SampleID), metadata_names])


  # - - - - - - - - - - - - -
  # Export table
  # - - - - - - - - - - - - -
  # Error handling::only one metadata column input.
  if(length(metadata_names)==1){
    colnames(metadata_tbl) <- metadata_names
    dplyr::bind_cols(input.df, metadata_tbl) -> input.df
  }
  else{
    dplyr::bind_cols(input.df, metadata_tbl) -> input.df
  }
  return(input.df)
}
