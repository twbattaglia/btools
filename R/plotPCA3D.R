# @ Thomas W. Battaglia

#' Plot DESeq2's PCA plotting with Plotly 3D scatterplot
#'
#' The function will generate a plot_ly 3D scatter plot image for
#' a 3D exploration of the PCA.
#'
#' @param object a DESeqTransform object, with data in assay(x), produced for example by either rlog or varianceStabilizingTransformation.
#' @param intergroup interesting groups: a character vector of names in colData(x) to use for grouping
#' @param ntop number of top genes to use for principal components, selected by highest row variance
#' @param returnData should the function only return the data.frame of PC1, PC2 and PC3 with intgroup covariates for custom plotting (default is FALSE)
#' @return An object created by plot_ly, which can be assigned and further customized.
#' @export
plotPCA3D <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE){
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1],
                  PC2 = pca$x[, 2],
                  PC3 = pca$x[, 3],
                  group = group,
                  intgroup.df,
                  name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  message("Generating plotly plot")
  return(
    plotly::plot_ly(data = d,
                    x = ~PC1,
                    y = ~PC2,
                    z = ~PC3,
                    mode = "markers",
                    type = "scatter3d")
  )
}
