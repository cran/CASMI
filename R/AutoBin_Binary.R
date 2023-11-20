# return suggested threshold for a quantitative variable
autoBinThresh <- function(data, colIndex){ # column index of X that needs auto binning

  finalThresh <- NA

  if (colIndex >= ncol(data) || colIndex <= 0) {
    return(finalThresh)
  }

  # deal with NAs, equal values, and other problems
  if ( (sum(!is.na(data[,colIndex])) <= 1) || (all(na.omit(data[,colIndex]) == na.omit(data[,colIndex])[1])) || !is.numeric(data[,colIndex]) ) {
    return(finalThresh)
  }


  sortedXY <- na.omit( data[order(data[,colIndex]), c(colIndex,ncol(data))] ) # subset, sort by X
  maxMI <- 0
  for(i in 1:(nrow(sortedXY) - 1)){
    threshold <- (sortedXY[i,1] + sortedXY[i+1,1]) / 2
    tmpMI <- MI.z(table(cut(sortedXY[,1], breaks = c(-Inf, threshold, Inf), labels = c("L", "H")), sortedXY[,ncol(sortedXY)]))
    if(tmpMI > maxMI){
      finalThresh <- threshold
      maxMI <- tmpMI
    }
  }
  return(finalThresh)

}

#' Auto Binning for Quantitative Variables - Binary
#'
#' Automatically suggest the optimal cutting point for categorizing a quantitative variable before using the \pkg{CASMI}-based functions. This function does binary cutting, that is, to convert the quantitative variable into a categorical variable with two levels/categories.
#' @param data data frame. An outcome variable is required. The outcome variable (Y) must be in the last column.
#' @param index index or a vector of indices of the quantitative variable(s) that need(s) to be automatically categorized.
#' @return `autoBin.binary()` returns the entire data frame after automatic binary categorization for the selected quantitative variable(s).
#' @examples
#' ## Use the "iris" dataset embedded in R
#' data("iris")
#' autoBin.binary(iris, c(1,2,3,4))
#'
#' @importFrom EntropyEstimation MI.z
#' @importFrom stats na.omit
#'
#' @export


# return the finalized data frame after auto binning
autoBin.binary <- function(data, index){
  # Check if the inputs are of correct type
  if (!is.data.frame(data)) {
    stop("The 'data' input must be a dataframe.")
  }
  if (!is.numeric(index)) {
    stop("The 'index' input must be an index or a vector of indices.")
  }

  for(i in 1:length(index)){
    colIndex <- index[i]
    finalThresh <- autoBinThresh(data, colIndex)
    if(!is.na(finalThresh)){
      data[,colIndex] <- cut(data[,colIndex], breaks = c(-Inf, finalThresh, Inf)) # cut the original variable by the threshold
    } else {
      sentence1 <- "The automatic categorization is not done for this variable -- Column Index: "
      sentence1 <- paste0(sentence1, colIndex)
      sentence2 <- " -- due to the following possible problems: (1) the column is not of a numeric type; (2) the column index is not valid; (3) there is none or only one (distinct) value in the following variable, so no categorization is needed."
      sentence1 <- paste0(sentence1, sentence2)
      warning(sentence1)
    }
  }

  return(as.data.frame(data))
}

