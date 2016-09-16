# @ Thomas W. Battaglia

#' Import .RCC files
#'
#' Import nanostring RCC raw files and generate a table based on sampleID's from each file.
#'
#' @param rccPath .RCC files from NanoString analysis. Adapted from NanoStringNorm.
#' @param useFolderPath Is the input rccPath a real folder, or from a dataframe.
#' @return A dataframe of count data with columns as samples and rows as genes.
#' @export
#'
import_rcc <- function(rccPath = "", useFolderPath = T){

  if(useFolderPath){
    if (!file.exists(rccPath)) {
      stop(paste("READ.MARKUP.RCC: path was not found.  \n"))
    }
    rcc.files <- list.files(path = rccPath, pattern = "*.RCC|*.rcc", full.names = T)
    if (length(rcc.files) == 0) {
      stop(paste("READ.MARKUP.RCC: no RCC files found.  \n"))
    }
  } else {
    rcc.files <- rccPath
  }
  rcc.header.merged <- NULL
  rcc.data.merged <- NULL
  count = 1
  for (rcc.file in rcc.files) {
    rcc.header <- read.table(rcc.file, nrows = 15, comment.char = "<", sep = ",", as.is = TRUE)
    rcc.data <- read.table(rcc.file, skip = 25, header = TRUE, comment.char = "<", sep = ",", as.is = TRUE, nrows = -1)

    sample.name <- rcc.header$V2[3]
    #sample.name <- gsub(".RCC", "", gsub(" ", "_", rcc.file))
    colnames(rcc.header)[2] <- sample.name
    colnames(rcc.data)[4] <- sample.name

    if (count == 1) {
      rcc.header.merged <- rcc.header
      rcc.data.merged <- rcc.data
    }
    else {
      rcc.header.merged <- data.frame(rcc.header.merged, subset(rcc.header, select = 2))
      rcc.data.merged <- data.frame(rcc.data.merged, subset(rcc.data, select = 4))
    }
    count = count + 1
  }
  rcc.header.merged[3, 1] <- "sample.id"
  rcc.header.merged[9, 1] <- "lane.id"
  rownames(rcc.header.merged) <- rcc.header.merged[, 1]
  rcc.header.merged <- rcc.header.merged[, -1]
  colnames(rcc.data.merged[1]) <- "Code.Count"
  x <- list(x = rcc.data.merged, header = rcc.header.merged)
  class(x) <- "NanoString"
  return(rcc.data.merged)
}
