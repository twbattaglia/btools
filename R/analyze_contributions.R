# @ Thomas W. Battaglia

#' Add meta data to contributions file
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param contributions_fp File location of the metagenomic_contributions.py output.
#' @param mappingfile_fp File location of the input sample metadata.
#' @return A very long dataframe with taxa information and sample metadata for each observed OTU.
#' @export
analyze_contributions <- function(contributions_fp, mappingfile_fp){

  # - - - - - - - - - - - - -
  # Error handling
  # - - - - - - - - - - - - -
  # input is a string
  # input is a data.frame
  # mapping file is a string
  # mapping file is a dataframe
  # column id's are not strings
  # First column is not #SampleID


  # - - - - - - - - - - - - -
  # Import files
  # - - - - - - - - - - - - -
  message("Importing files...")
  input <- read.delim(contributions_fp, header = TRUE)
  gg_db <- otu_taxonomy_97
  mappingfile <- read.delim(mappingfile_fp, header = TRUE)
  input.df <- dplyr::tbl_df(input)
  input.df$OTU <- as.factor(input.df$OTU)


  # - - - - - - - - - - - - -
  # Process strings
  # - - - - - - - - - - - - -
  message("Converting identifiers to names...")
  otuid_name <- gg_db[match(input.df$OTU, table = gg_db$V1), ]$V2
  otuid_name <- gsub(pattern = ';', replacement = "", x = otuid_name, fixed = T)
  otuid_name_sep <- sapply(otuid_name, function(x) stringr::str_split(x, " "))


  # - - - - - - - - - - - - -
  # Adding Names to Table
  # - - - - - - - - - - - - -
  message("Adding names to table...")
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
  message("Adding metadata to table...")
  input.df.metadata <- merge(input.df, metadata, by.x = "Sample", by.y = "X.SampleID")


  # - - - - - - - - - - - - -
  # Export table
  # - - - - - - - - - - - - -
  return(input.df.metadata)
}
