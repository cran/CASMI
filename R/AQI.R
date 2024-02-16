#MI.test tests the independence between a feature and the outcome. It prints out p.value, the smaller the p.value, the stronger evidence of dependence between them.
MI.test=function(feature, outcome, k1, k2){ #feature is the vector of X, outcome is the vector of Y, k1 and k2 are the corresponding number of categories, X and Y must have the same sample size.
  n=length(feature);
  test.stat=2*n*MI.z(table(feature, outcome))+(k1-1)*(k2-1);
  p.value=pchisq(test.stat,(k1-1)*(k2-1),lower.tail = F);
  return(p.value);
}

#' AQI Index
#'
#' A quantitative measure of dataset quality. The AQI Index score indicates the degree that how features are associated with the outcome in a dataset. (synonyms of "feature": "variable" "factor" "attribute") \cr
#' For more information, please refer to the corresponding publication: Shi, J., Zhang, J. and Ge, Y. (2019), "An Association-Based Intrinsic Quality Index for Healthcare Dataset Ranking" <doi:10.1109/ICHI.2019.8904553>
#' @param data data frame (features as columns and observations as rows). The outcome variable (Y) MUST be the last column. It requires at least one features and only one outcome. Both the features (Xs) and the outcome (Y) MUST be discrete (if not naturally discrete, you may try the `autoBin.binary` function in the same package).
#' @param alpha.filter level of significance for the mutual information test of independence in step 2 (<doi:10.1109/ICHI.2019.8904553>). By default, `alpha.filter = 0.2`.
#' @return The AQI Index score.
#' @importFrom EntropyEstimation Entropy.z MI.z
#' @examples
#' ## Generate a toy dataset: "data"
#' n <- 10000
#' set.seed(1)
#' x1 <- rbinom(n, 3, 0.5) + 0.2
#' set.seed(2)
#' x2 <- rbinom(n, 2, 0.8) + 0.5
#' set.seed(3)
#' x3 <- rbinom(n, 5, 0.3)
#' set.seed(4)
#' error <- round(runif(n, min=-1, max=1))
#' y <- x1 + x3 + error
#' data <- data.frame(cbind(x1, x2, x3, y))
#' colnames(data) <- c("feature1", "feature2", "feature3", "Y")
#'
#' ## Calculate the AQI score of "data"
#' AQI(data)
#' @export

AQI <- function(data, alpha.filter=0.2){ #outcome must be in the last column

  data <- as.data.frame(data)

  score=NULL

  # Step 1: return dependent features, step1Index[]

  n=nrow(data)

  u=1/log(log(n))

  step1Index=vector()
  count=0
  k2=length(table(data[length(data)])) #outcome categories
  for(fi in 1:(length(data)-1)){
    k1=length(table(data[fi])) #feature categories
    p_value=MI.test(data[[fi]],data[[length(data)]], k1, k2)
    if(p_value<=alpha.filter){
      count=count+1
      step1Index[count]=fi
    }
  }

  indexCASMI=vector()

  if(count==0){
    warning("No associated variables are found toward the outcome, thus the AQI score is zero. First, please ensure that the outcome variable is the last column. Second, please ensure that a reasonable alpha.filter is set, or use the default value.")
    score=0
  } else{
    # Step 2: select features by joint SMI with Hz estimator
    # select the best feature
    maxKappaStar=0
    maxIndex=0
    for(index in step1Index){
      #valueKappaStar=kappa.star(data[[index]], data[[length(data)]])
      feature=data[[index]]; outcome=data[[length(data)]];

      #valueKappaStar=MI.z(table(feature, outcome))/Entropy.z(table(outcome))*(1-length(which(table(feature)==1))/n)
      valueKappaStar=MI.z(table(feature, outcome))/Entropy.z(table(outcome))*(1-sum(table(feature)==1L)/n)^u

      if(valueKappaStar > maxKappaStar){
        maxKappaStar=valueKappaStar
        maxIndex=index
      }
    }
    indexCASMI=c(indexCASMI, maxIndex)

    if(length(step1Index)==1) {return(100/(1-(log(maxKappaStar))/exp(1)))}


    # select the 2nd, 3rd, ... best features by joint

    maxmaxKappaStar=0
    while(maxKappaStar>maxmaxKappaStar & length(indexCASMI)<count){
      maxmaxKappaStar=maxKappaStar

      step1Index=step1Index[!step1Index==maxIndex]
      maxKappaStar=0
      maxIndex=0
      for(index in step1Index){
        tmpIndex1=c(index, indexCASMI, length(data))
        tmpIndex2=c(index, indexCASMI)
        ftable=ftable(table(data[tmpIndex1]))
        ftableOfFeatures=ftable(table(data[tmpIndex2]))
        #valueKappaStar=kappa.star2(ftable, ftableOfFeature, data[[length(data)]])
        outcome=data[[length(data)]]

        #valueKappaStar=MI.z(ftable)/Entropy.z(table(outcome))*(1-length(which(ftableOfFeatures==1))/n)
        valueKappaStar=MI.z(ftable)/Entropy.z(table(outcome))*(1-sum(ftableOfFeatures==1L)/n)^u

        if(valueKappaStar > maxKappaStar){
          maxKappaStar=valueKappaStar
          maxIndex=index
        }
      }

      if(maxKappaStar>maxmaxKappaStar+10^-14){ # +10^-14 is to solve the problem of precision
        indexCASMI=c(indexCASMI, maxIndex)
      }
    }
    score=maxmaxKappaStar
  }

  return(100/(1-(log(score))/exp(1)))
}
