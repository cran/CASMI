# return suggested threshold for a quantitative variable
autoBinThresh <- function(data, colIndex){ # column index of X that needs auto binning

  # deal with NAs
  if (sum(!is.na(data[,colIndex])) <= 1) {
    stop("There is none or only one available value in a selected variable. No categorization is needed. Please double-check your input variables.")
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
  finalThresh

  return(finalThresh)
}

#' Auto Binning for Quantitative Variables - Binary
#'
#' Automatically suggest the optimal cutting point for categorizing a quantitative variable before using the \pkg{CASMI}-based functions. This function does binary cutting, that is to convert the quantitative variable into a categorical variable with two levels/categories.
#' @param data data frame (features as columns and observations as rows). It requires at least one feature and only one outcome. The outcome variable (Y) must be in the last column.
#' @param index index or a vector of indices of the quantitative variable(s) that need(s) to be automatically categorized.
#' @return `autoBin.binary()` returns the entire data frame after auto binning for the selected quantitative variable(s).
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
    stop("The first input must be a dataframe.")
  }
  if (!is.numeric(index)) {
    stop("The second input must be an index or a vector of indices.")
  }

  for(i in 1:length(index)){
    colIndex <- index[i]
    finalThresh <- autoBinThresh(data, colIndex)
    data[,colIndex] <- cut(data[,colIndex], breaks = c(-Inf, finalThresh, Inf)) # cut the original variable by the threshold
  }

  return(data)
}

