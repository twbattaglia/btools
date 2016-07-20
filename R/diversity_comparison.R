# @ Thomas W. Battaglia

#' Compute pairwise diversity comparisons of different beta diversity metrics.
#'
#' The function will generate the distance comparisons of various estimates found in the vegan package.
#' The output will be in a dataframe format to generate additionally analyses of the distance
#' comparisons.
#'
#' @param phylo An input phyloseq object.
#' @param distance Distance metrics to use for the pairwise comparisons.
#' @return A dataframe of distance comparions.
#' @export
diversity_comparison <- function(phylo, distance = c("jsd", "unifrac", "wunifrac", "dpcoa")){

  # Get phyloseq object metadata
  metadata <- as(phyloseq::sample_data(phylo), "data.frame")

  # Get Jensen-Shannon Index
  dist = phyloseq::distance(phylo, distance)

  # Convert to dataframe object
  dist.df <- reshape2::melt(as.matrix(dist), varnames = c("sample1", "sample2"))

  # Remove 0 columns
  dist.df <- subset(dist.df, value != 0)

  # Convert to factors
  dist.df$sample1 <- as.factor(dist.df$sample1)
  dist.df$reference <- as.factor(dist.df$sample2)

  # Return dataframe
  return(dist.df)
}
