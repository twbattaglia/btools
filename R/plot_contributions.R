# @ Thomas W. Battaglia

#' Plot metagenomic contributions data
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param output from metagneomic contributions processing function
#' @param mapping file with sample metadata
#' @param list of column variables to add to new table
#' @return A dataframe with taxa information and sample metadata
#' @export
