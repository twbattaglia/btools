# @ Thomas W. Battaglia

#' Plot taxa abundance barplots with a threshold.
#'
#' Plot a barplot between 2 or more groups, and combining those
#' OTU's with low abundance into a group 'Other', to reduce the amount of
#' taxa to plot.
#'
#' @param phylo A phyloseq object
#' @param group The column variable that explains the groupings/x-axis
#' @param threshold The % abundance threshold to use to classify if a taxa is considered 'Other'
#' @param level The phylogeny level to use for collapsing the OTU table.
#' @param prism Should the data be exported to be used for PRISM.
#' @export
