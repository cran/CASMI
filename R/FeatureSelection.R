# MIz.test tests the independence between a feature and the outcome (with z estimator). It prints out p.value, the smaller the p.value, the stronger evidence of dependence between them.
MIz.test <- function(feature, outcome) {
  #feature is the vector of X, outcome is the vector of Y, k1 and k2 are the corresponding number of categories, X and Y must have the same sample size.

  # considering NA in the feature
  XYtable <- table(feature, outcome)
  n <- sum(XYtable)

  k1 = length(table(feature)) # feature categories
  k2 = length(table(outcome)) # outcome categories

  test.stat = 2 * n * MI.z(XYtable) + (k1 - 1) * (k2 - 1)

  p.value = pchisq(test.stat, (k1 - 1) * (k2 - 1), lower.tail = F)

  return(p.value)

}


g_function <- function(x1, x2, x3) {
  # return a 3*1 matrix
  vec = c(1 / x2,-(x1 - x3) / x2 ^ 2,-1 / x2)
  return(as.matrix(vec, nrow = 3, ncol = 1))
}


# confidence interval (CI) of SMI, whatever the estimator is. return E in CI.
E.CI.smi <- function(XYtable, alpha) {
  # input the frequency table between the feature and the outcome

  K1 = nrow(XYtable)
  K2 = ncol(XYtable)
  K = K1 * K2
  n = sum(XYtable)


  # make sure XYtable[K1,K2] is not zero. If zero, switch.
  if (XYtable[K1, K2] == 0) {
    nonZeroIndex = which(XYtable != 0, arr.ind = T)
    i = nonZeroIndex[1, 1]
    j = nonZeroIndex[1, 2]

    tmp = XYtable[i,]
    XYtable[i,] = XYtable[K1,]
    XYtable[K1,] = tmp

    tmp = row.names(XYtable)[i]
    row.names(XYtable)[i] = row.names(XYtable)[K1]
    row.names(XYtable)[K1] = tmp

    tmp = XYtable[, j]
    XYtable[, j] = XYtable[, K2]
    XYtable[, K2] = tmp

    tmp = colnames(XYtable)[j]
    colnames(XYtable)[j] = colnames(XYtable)[K2]
    colnames(XYtable)[K2] = tmp
  }


  Xmargin = rowSums(XYtable)
  Ymargin = colSums(XYtable)


  # first row of matrix A
  pK1. = Xmargin[K1] / sum(Xmargin)
  derX = vector()
  for (i in 1:(K1 - 1)) {
    pi. =  Xmargin[i] / sum(Xmargin)
    derX[i] = log(pK1.) - log(pi.)
  }

  A1 = vector()
  for (i in 1:(K1 - 1)) {
    for (j in 1:K2) {
      A1 = c(A1, derX[i])
    }
  }
  for (i in ((K1 - 1) * K2 + 1):(K - 1)) {
    A1[i] = 0
  }

  # second row of matrix A
  p.K2 = Ymargin[K2] / sum(Ymargin)
  derY = vector()
  for (i in 1:(K2 - 1)) {
    p.j = Ymargin[i] / sum(Ymargin)
    derY[i] = log(p.K2) - log(p.j)
  }

  A2 = vector()
  for (i in 1:K1) {
    for (j in 1:K2) {
      if (j == K2) {
        A2 = c(A2, 0)
      } else{
        A2 = c(A2, derY[j])
      }
    }
  }
  A2 = A2[-length(A2)]

  # third row of matrix A
  pK = XYtable[K1, K2] / n
  A3 = vector()
  for (i in 1:K1) {
    for (j in 1:K2) {
      pk = XYtable[i, j] / n
      if (pk == 0) {
        A3 = c(A3, 0)
      } else{
        A3 = c(A3, log(pK) - log(pk))
      }
    }
  }
  A3 = A3[-length(A3)]

  A = as.matrix(rbind(A1, A2, A3))


  # the Sigma matrix
  p = c(t(XYtable)) / n
  Sigma = matrix(
    data = NA,
    nrow = K - 1,
    ncol = K - 1,
    byrow = TRUE
  )
  for (i in 1:(K - 1)) {
    for (j in 1:(K - 1)) {
      if (i == j) {
        Sigma[i, i] = p[i] * (1 - p[i])
      } else{
        Sigma[i, j] = -p[i] * p[j]
      }
    }
  }



  # sigma hat
  HXhat = entropy.plugin(Xmargin)
  HYhat = entropy.plugin(Ymargin)
  HXYhat = entropy.plugin(XYtable)
  g = g_function(HXhat, HYhat, HXYhat)
  sigmaHat = sqrt(t(g) %*% A %*% Sigma %*% t(A) %*% g)



  # E in CI
  Zstat = qnorm(alpha / 2, lower.tail = FALSE) # 95% confidence level
  E = Zstat * (sigmaHat / sqrt(n))

  return(E)
}

ftable2xytable <- function(ftable) { # convert an ftable to a regular XYtable where X are joint and Y is the outcome
  return(ftable[rowSums(ftable)!=0,])
}

Kappas <- function(data, idxFeatures){
  # input data, index(es) of feature(s); outcome must be in the last column

  # ftable automatically handles NA values
  tmpIndexFnO = c(idxFeatures, ncol(data)) # Index of features and outcome
  ftable = ftable(table(data[tmpIndexFnO]))

  Xmargin <- rowSums(ftable)
  Ymargin <- colSums(ftable)
  n <- sum(ftable)

  kappa.z = MI.z(ftable) / Entropy.z(Ymargin) # same as smi.z

  kappaStar = kappa.z * (1 - sum(Xmargin == 1L) / n)

  return(list(kappa.z = kappa.z, kappa.star = kappaStar))
}


generateResults <-
  function(data, indexCASMI, kappaStarCASMI, alpha) {
    # return a set of results based on indexCASMI

    nameCASMI = names(data[indexCASMI])

    # 1. SMIziCASMI: standardized mutual information (with z estimator, based on H(Y)) of each selected feature. Not cumulative.
    # 2. PziCASMI: p-value of mutual information test (with z estimator) for each selected feature. (to test independence, prefer large n)
    SMIziCASMI = vector()
    SMIziCASMI.CI = vector()
    PziCASMI = vector()
    for (i in 1:length(indexCASMI)) {
      XYtable = table(data[,indexCASMI[i]], data[,ncol(data)]) # XYtable handled NA automatically
      Ymargin = colSums(XYtable)

      mi_z = MI.z(XYtable) # EntropyEstimation package
      smi_z = mi_z / Entropy.z(Ymargin)
      SMIziCASMI = c(SMIziCASMI, smi_z)

      low = smi_z - E.CI.smi(XYtable, alpha)
      uppr = smi_z + E.CI.smi(XYtable, alpha)
      SMIziCASMI.CI <- rbind(SMIziCASMI.CI, c(low, uppr))

      n = sum(XYtable)
      K1 = nrow(XYtable)
      K2 = ncol(XYtable)

      PziCASMI = c(PziCASMI, pchisq(
        2 * n * mi_z + (K1 - 1) * (K2 - 1),
        df = (K1 - 1) * (K2 - 1),
        lower.tail = FALSE
      )) # same as MIz.test()
    }



    #######

    sample_size <- apply(data[indexCASMI], 2, function(col) sum(!is.na(col)))


    df <-
      data.frame(
        indexCASMI,
        sample_size,
        round(kappaStarCASMI, 4),
        round(SMIziCASMI, 4),
        round(SMIziCASMI.CI, 4),
        sprintf("%.4f", round(PziCASMI, 4)),
        nameCASMI,
        row.names = NULL
      )

    colnames(df) <-
      c(
        "Var.Idx",
        "n",
        "Kappa*",
        "SMIz",
        "SMIz.Low",
        "SMIz.Upr",
        "p.MIz",
        "Var.Name"
      )



    ### overall CASMI score for all selected features
    kappa.star.hat <- Kappas(data, indexCASMI)$kappa.star
    tmpIndexFnO = c(indexCASMI, ncol(data)) # Index of features and outcome
    ftable = ftable(table(data[tmpIndexFnO]))
    tmpXYtable <- ftable2xytable(ftable)
    E <- E.CI.smi(tmpXYtable, alpha)
    low = kappa.star.hat - E
    low = round(low, 6)
    uppr = kappa.star.hat + E
    uppr = round(uppr, 6)
    ci.kappa.star.hat = c(low, uppr)

    kappa.star.hat = round(kappa.star.hat, 6)
    names(kappa.star.hat) = c("Kappa* for all selected features")


    names(ci.kappa.star.hat) = c("Kappa* CI Lower", "Upper")


    if (kappa.star.hat != round(kappaStarCASMI[length(kappaStarCASMI)],6)) {
      warning(
        "Mismatched Kappa* values detected. There may be an issue with the data or the function. Please report the issue to the package developer after ensuring the data has been properly preprocessed."
      )

    }


    ##### return

    return(
      list(
        Outcome.Variable.Name = names(data[ncol(data)]),
        Confidence.Level = 1 - alpha,
        KappaStar = kappa.star.hat,
        KappaStarCI = ci.kappa.star.hat,
        results = df
      )
    )

  }

# Function to check for long values
check_length <- function(column) {
  column = na.omit(column)
  return(any(nchar(as.character(column)) > 100))
}

#' \pkg{CASMI}-Based Feature Selection
#'
#' Selects the most relevant features toward an outcome. It automatically learns the number of features to be selected, and the selected features are ranked. The method automatically handles the feature redundancy issue. (synonyms of "feature": "variable" "factor" "attribute") \cr
#' For more information, please refer to the corresponding publication: Shi, J., Zhang, J. and Ge, Y. (2019), "\pkg{CASMI}—An Entropic Feature Selection Method in Turing’s Perspective" <doi:10.3390/e21121179>
#' @param data data frame (features as columns and observations as rows). The outcome variable (Y) MUST be the last column. It requires at least one features and only one outcome. Both the features (Xs) and the outcome (Y) MUST be discrete (if not naturally discrete, you may try the `autoBin.binary` function in the same package).
#' @param feature.na.handle options for handling NA values in the data. There are three options: `"stepwise", "na.omit", "NA as a category"`. `feature.na.handle = "stepwise"` excludes NA rows only when a particular variable is being calculated. For example, suppose we have data(Feature1: A, NA, B; Feature2: C, D, E; Feature3: F, G, H; Outcome: O, P, Q); the second observation will be excluded only when a particular step includes Feature1, but will not be excluded when a step calculates among Feature2, Feature3, and the Outcome. This option is designed to take advantage of a maximum number of data points. `feature.na.handle = "na.omit"` excludes observations with any NA values at the beginning of the analysis. `feature.na.handle = "NA as a category"` regards the NA value as a new category. This is designed to be used when NA values in the data have a consistent meaning instead of being missing values. For example, in survey data asking for comments, each NA value might consistently mean "no opinion." By default, `feature.na.handle = "stepwise"`.
#' @param alpha.filter level of significance for the mutual information test of independence in step 1 of the features selection (initial screening). The smaller the alpha.filter, the fewer the features sent to step 2 (<doi:10.3390/e21121179>). By default, `alpha.filter = 0.1`.
#' @param alpha level of significance for the confidence intervals in final results. By default, `alpha = 0.05`.
#' @param intermediate.steps output the intermediate process. By default, `intermediate.steps = TRUE`. Set to `FALSE` for showing only summary results.
#' @param kappa.star.cap a threshold of `kappa*` for pausing the feature selection process. The program will automatically pause at the first feature of which the `kappa*` value exceeds the kappa.star.cap threshold. By default, `kappa.star.cap = 1.0`, which is the maximum possible value. A lower value may result in fewer final features but less computing time.
#' @param feature.num.cap the maximum number of features to be selected. A lower value may result in fewer final features but less computing time.
#' @return `CASMI.selectFeatures()` returns selected features and relevant information, including the estimated Kappa* for all selected features (`$KappaStar`) and the corresponding confidence interval (`$KappaStarCI`). The selected features are ranked. The Standardized Mutual Information using the z estimator (`SMIz`) and the corresponding confidence interval (`SMIz.Low` for lower bound, `SMIz.Upr` for upper bound) are given for each selected feature (`Var.Idx` for column index, `Var.Name` for column name). The p-value from the mutual information test of independence using the z estimator (`p.MIz`) is given for each selected feature.
#' @examples
#' ## Generate a toy dataset: "data"
#' ## Features 1 and 3 are associated with Y, while feature 2 is irrelevant.
#' ## The outcome variable Y must be discrete and be the LAST column. Features must be discrete.
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
#' ## Select features and provide relevant results for the toy dataset "data"
#' CASMI.selectFeatures(data)
#'
#' ## For showing only the summary results
#' CASMI.selectFeatures(data, intermediate.steps = FALSE)
#'
#' ## Adjust 'feature.num.cap' for including fewer features.
#' ## A lower 'feature.num.cap' value may result in fewer final features but less computing time.
#' ## For example, if needing only the top one feature based on the toy dataset:
#' CASMI.selectFeatures(data, feature.num.cap = 1)
#'
#'
#' @importFrom EntropyEstimation Entropy.z MI.z
#' @importFrom entropy entropy.plugin mi.plugin
#' @importFrom stats pchisq qnorm
#'
#' @export



CASMI.selectFeatures <- function(data,  # outcome must be in the last column
                                 feature.na.handle = "stepwise",
                                 alpha.filter = 0.1,
                                 alpha = 0.05,
                                 intermediate.steps = TRUE,
                                 kappa.star.cap = 1.0,
                                 feature.num.cap = ncol(data)) {

  # check data type
  if (!is.data.frame(data)) {
    stop("Error: The input is not a data frame type.")
  }


  data <- as.data.frame(data)



  outcome.has.na <- any(is.na(data[,ncol(data)]))
  if (outcome.has.na) {
    stop("Error: The last column has NA values. Please make sure to put the outcome in the last column. NA is not allowed in the outcome.")
  }



  # Check each column
  result <- apply(data, 2, check_length)
  # Print columns with at least one value over 100 characters/digits
  long_columns <- names(result)[result]
  # Asking user
  if(length(long_columns) > 0){
    cat("The following columns have values with over 100 characters/digits:\n")
    print(long_columns)

    # Prompt the user for further action
    user_choice <- readline(prompt="Do you want to continue with this data? (yes/no): ")

    if(tolower(user_choice) == "no"){
      stop("Ok. Please pre-process the data and try again.")
    } else if(tolower(user_choice) != "yes"){
      stop("Invalid choice. Assuming you don't want to continue.")
    }
  }



  # NA handling
  if(feature.na.handle == "stepwise") {
    # do nothing
  }else if(feature.na.handle == "na.omit") {
    data <- na.omit(data)
  }else if(feature.na.handle == "NA as a category") {
    data[] <- apply(data, 2, function(column) {
      ifelse(is.na(column), "=na=", column)
    })
    data <- as.data.frame(data)
  }else {
    stop("Error: Wrong 'feature.na.handle' value. Please check the documentation.")
  }


  if (!is.numeric(alpha.filter)) {
    stop("Error: 'alpha.filter' should be numeric.")
  } else if (alpha.filter < 0 || alpha.filter > 1) {
    stop("Error: 'alpha.filter' should be between 0 and 1.")
  }

  if (!is.numeric(alpha)) {
    stop("Error: 'alpha' should be numeric.")
  } else if (alpha < 0 || alpha > 1) {
    stop("Error: 'alpha' should be between 0 and 1.")
  }

  if (!is.logical(intermediate.steps) || length(intermediate.steps) != 1) {
    stop("Error: 'intermediate.steps' must be either TRUE or FALSE.")
  }

  if (!is.numeric(kappa.star.cap)) {
    stop("Error: 'kappa.star.cap' should be numeric.")
  } else if (kappa.star.cap < 0 || kappa.star.cap > 1) {
    stop("Error: 'kappa.star.cap' should be between 0 and 1.")
  }

  if (!is.numeric(feature.num.cap)) {
    stop("Error: 'feature.num.cap' should be numeric.")
  } else if (feature.num.cap < 1) {
    stop("Error: 'feature.num.cap' should be at least 1.")
  }




  start_time <- proc.time()


  # Step 1: return dependent features, step1Index[]

  step1Index = vector()
  count = 0
  for (fi in 1:(ncol(data) - 1)) {
    p_value = MIz.test(data[,fi], data[,ncol(data)])
    if (p_value <= alpha.filter) {
      count = count + 1
      step1Index[count] = fi
    }
  }

  indexCASMI = vector()
  kappaStarCASMI = vector()

  if (count == 0) {
    warning(
      "No associated variables are found toward the outcome. Please double-check: 1. the outcome variable is the last column; 2. a reasonable 'alpha.filter' is set (by default it's 0.1)."
    )

    indexCASMI = NULL

  } else{

    if(intermediate.steps) {
      cat(count, "feature(s) passed Step 1 (a test of independence).\n")
    }


    # Step 2: select features by joint SMI with Hz estimator

    # select the best feature (the first feature)
    maxKappaStar = 0
    maxIndex = 0
    for (index in step1Index) {
      valueKappaStar = Kappas(data, index)$kappa.star

      if (valueKappaStar > maxKappaStar) {
        maxKappaStar = valueKappaStar
        maxIndex = index
      }
    }
    indexCASMI = c(indexCASMI, maxIndex)
    kappaStarCASMI = c(kappaStarCASMI, maxKappaStar)

    if(intermediate.steps) {
      intermediateCount = 1
      cat("Selecting Feature:", intermediateCount, "- selected column (idx.):", maxIndex, ", Kappa*:", maxKappaStar, ", elapsed (sec.):", (proc.time() - start_time)["elapsed"], "\n")
    }

    if (length(step1Index) == 1) {
      return(generateResults(data, indexCASMI, kappaStarCASMI, alpha))
    }



    # select the 2nd, 3rd, ... best features by joint
    maxmaxKappaStar = 0
    while (maxKappaStar > maxmaxKappaStar &
           length(indexCASMI) < count &
           maxKappaStar < kappa.star.cap &
           length(indexCASMI) < feature.num.cap) {

      maxmaxKappaStar = maxKappaStar

      step1Index = step1Index[!step1Index == maxIndex]
      maxKappaStar = 0
      maxIndex = 0
      for (index in step1Index) {
        tmpIndexFeatures = c(indexCASMI, index) # Index of selected features + the new feature

        valueKappaStar = Kappas(data, tmpIndexFeatures)$kappa.star

        if (valueKappaStar > maxKappaStar) {
          maxKappaStar = valueKappaStar
          maxIndex = index
        }
      }

      if (maxKappaStar > maxmaxKappaStar + 10 ^ -14) { # +10^-14 is to solve the problem of precision
        indexCASMI = c(indexCASMI, maxIndex)
        kappaStarCASMI = c(kappaStarCASMI, maxKappaStar)

        if(intermediate.steps) {
          intermediateCount = intermediateCount + 1
          cat("Selecting Feature:", intermediateCount, "- selected column (idx.):", maxIndex, ", Kappa*:", maxKappaStar, ", elapsed (sec.):", (proc.time() - start_time)["elapsed"], "\n")
        }
      }

      if (maxKappaStar >= kappa.star.cap) {
        warning(
          "The feature selection process paused automatically because kappa* of currently selected features reached the 'kappa.star.cap' value."
        )
      }

      if (length(indexCASMI) >= feature.num.cap) {
        warning(
          "The feature selection process paused automatically because the number of selected features reached the 'feature.num.cap' value."
        )
      }
    }

    if(intermediate.steps) {
      cat("---End of intermediate process.---\n")
      cat("In progress of generating summary.\n\n")
    }

    list = generateResults(data, indexCASMI, kappaStarCASMI, alpha)

    return(list)
  }
}

