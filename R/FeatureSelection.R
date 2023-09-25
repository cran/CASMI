# MIz.test tests the independence between a feature and the outcome (with z estimator). It prints out p.value, the smaller the p.value, the stronger evidence of dependence between them.
MIz.test=function(feature, outcome, k1, k2){ #feature is the vector of X, outcome is the vector of Y, k1 and k2 are the corresponding number of categories, X and Y must have the same sample size.
  n=length(feature);
  test.stat=2*n*MI.z(table(feature, outcome))+(k1-1)*(k2-1);
  p.value=pchisq(test.stat,(k1-1)*(k2-1),lower.tail = F);
  return(p.value);
}


g_function = function(x1, x2, x3){ # return a 3*1 matrix
  vec = c(1/x2, - (x1-x3)/x2^2, - 1/x2)
  return(as.matrix(vec, nrow = 3, ncol = 1))
}


# confidence interval (CI) of SMI, whatever the estimator is. return E in CI.
E.CI.smi <- function(XYtable, alpha=0.05){ # input the frequency table between the feature and the outcome

  K1 = nrow(XYtable)
  K2 = ncol(XYtable)
  K = K1*K2
  n = sum(XYtable)


  # make sure XYtable[K1,K2] is not zero. If zero, switch.
  if(XYtable[K1,K2] == 0){
    nonZeroIndex = which(XYtable != 0, arr.ind = T)
    i = nonZeroIndex[1,1]
    j = nonZeroIndex[1,2]

    tmp = XYtable[i,]
    XYtable[i,] = XYtable[K1,]
    XYtable[K1,] = tmp

    tmp = row.names(XYtable)[i]
    row.names(XYtable)[i] = row.names(XYtable)[K1]
    row.names(XYtable)[K1] = tmp

    tmp = XYtable[,j]
    XYtable[,j] = XYtable[,K2]
    XYtable[,K2] = tmp

    tmp = colnames(XYtable)[j]
    colnames(XYtable)[j] = colnames(XYtable)[K2]
    colnames(XYtable)[K2] = tmp
  }


  Xmargin = rowSums(XYtable)
  Ymargin = colSums(XYtable)


  # first row of matrix A
  pK1. = Xmargin[K1]/sum(Xmargin)
  derX = vector()
  for(i in 1:(K1-1)){
    pi. =  Xmargin[i]/sum(Xmargin)
    derX[i] = log(pK1.) - log(pi.)
  }

  A1 = vector()
  for(i in 1:(K1-1)){
    for(j in 1:K2){
      A1 = c(A1, derX[i])
    }
  }
  for(i in ((K1-1)*K2+1):(K-1)){
    A1[i] = 0
  }

  # second row of matrix A
  p.K2 = Ymargin[K2]/sum(Ymargin)
  derY = vector()
  for(i in 1:(K2-1)){
    p.j = Ymargin[i]/sum(Ymargin)
    derY[i] = log(p.K2) - log(p.j)
  }

  A2 = vector()
  for(i in 1:K1){
    for(j in 1:K2){
      if(j==K2){
        A2 = c(A2, 0)
      }else{
        A2 = c(A2, derY[j])
      }
    }
  }
  A2 = A2[-length(A2)]

  # third row of matrix A
  pK = XYtable[K1,K2]/n
  A3 = vector()
  for(i in 1:K1){
    for(j in 1:K2){
      pk = XYtable[i,j]/n
      if(pk == 0){
        A3 = c(A3, 0)
      }else{
        A3 = c(A3, log(pK) - log(pk))
      }
    }
  }
  A3 = A3[-length(A3)]

  A = as.matrix(rbind(A1, A2, A3))


  # the Sigma matrix
  p = c(t(XYtable))/n
  Sigma = matrix(data = NA, nrow = K-1, ncol = K-1, byrow = TRUE)
  for(i in 1:(K-1)){
    for(j in 1:(K-1)){
      if(i == j){
        Sigma[i,i] = p[i] * (1-p[i])
      }else{
        Sigma[i,j] = - p[i] * p[j]
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
  Zstat = qnorm(alpha/2, lower.tail = FALSE) # 95% confidence level
  E = Zstat * (sigmaHat / sqrt(n))

  return(E)
}


generateResults <- function(data, indexCASMI, kappaStarCASMI, alpha = 0.05){ # return a set of results based on indexCASMI
  nameCASMI = names(data[indexCASMI])

  # 1. SMIziCASMI: standardized mutual information (with z estimator, based on H(Y)) of each selected feature. Not cumulative.
  # 2. PziCASMI: p-value of mutual information test (with z estimator) for each selected feature. (to test independence, prefer large n)
  SMIziCASMI = vector()
  SMIziCASMI.CI = vector()
  PziCASMI = vector()
  for(i in 1:length(indexCASMI)){
    XYtable = table(data[[indexCASMI[i]]], data[[length(data)]])
    Ymargin = colSums(XYtable)

    mi_z = MI.z(XYtable) # EntropyEstimation package
    smi_z = mi_z/Entropy.z(Ymargin)
    SMIziCASMI = c(SMIziCASMI, smi_z)

    low = smi_z - E.CI.smi(XYtable, alpha)
    uppr = smi_z + E.CI.smi(XYtable, alpha)
    SMIziCASMI.CI <- rbind(SMIziCASMI.CI, c(low, uppr))

    n = nrow(data)
    K1 = nrow(XYtable)
    K2 = ncol(XYtable)


    PziCASMI = c(PziCASMI, pchisq(2*n*mi_z + (K1-1)*(K2-1), df = (K1-1)*(K2-1), lower.tail = FALSE)) # same as MIz.test()
  }



  #######


  df <- data.frame(indexCASMI, nameCASMI, round(kappaStarCASMI,6), round(SMIziCASMI,6), round(SMIziCASMI.CI, 6), round(PziCASMI,6))
  colnames(df) <- c("Var.Index", "Var.Name", "Kappa*", "SMIz", "CI.SMIz.Lower", "CI.SMIz.Upper",  "P-value.MIz")



  ### overall CASMI score for all selected features
  combinedFeatures = do.call(paste0, data[indexCASMI])
  jointTable = cbind(combinedFeatures, data[length(data)])
  kappa.z.hat = MI.z(table(jointTable)) / Entropy.z(table(data[length(data)]))
  kappa.star.hat = kappa.z.hat * (1 - sum(table(combinedFeatures)==1L)/n)
  low = kappa.star.hat - E.CI.smi(table(jointTable), alpha)
  low = round(low, 6)
  uppr = kappa.star.hat + E.CI.smi(table(jointTable), alpha)
  uppr = round(uppr, 6)
  ci.kappa.star.hat = c(low, uppr)

  kappa.star.hat = round(kappa.star.hat,6)
  names(kappa.star.hat) = c("Kappa* for all selected features")


  names(ci.kappa.star.hat) = c("Kappa* CI Lower", "Upper")



  ##### return

  return(list(Confidence.Level = 1 - alpha, KappaStar = kappa.star.hat, KappaStarCI = ci.kappa.star.hat, results = df))

}

#' \pkg{CASMI}-Based Feature Selection
#'
#' Selects the most relevant features toward an outcome. It automatically learns the number of features to be selected, and the selected features are ranked. The method automatically handles the feature redundancy issue. (synonyms of "feature": "variable" "factor" "attribute") \cr
#' For more information, please refer to the corresponding publication: Shi, J., Zhang, J. and Ge, Y. (2019), "\pkg{CASMI}—An Entropic Feature Selection Method in Turing’s Perspective" <doi:10.3390/e21121179>
#' @param data data frame (features as columns and observations as rows). It requires at least one feature and only one outcome. The features must be discrete. The outcome variable (Y) must be in the last column.
#' @param alpha.filter level of significance for the mutual information test of independence in step 1 of the features selection. The smaller the alpha.filter, the fewer the features sent to step 2 (<doi:10.3390/e21121179>). By default, `alpha.filter = 0.1`.
#' @param alpha level of significance for the confidence intervals in final results. By default, `alpha = 0.05`.
#' @return `CASMI.selectFeatures()` returns selected features and relevant information, including the estimated Kappa* for all selected features (`$KappaStar`) and the corresponding confidence interval (`$KappaStarCI`). The selected features are ranked. The Standardized Mutual Information using the z estimator (`SMIz`) and the corresponding confidence interval (`CI.SMIz.Lower` and `CI.SMIz.Upper`) are given for each selected feature. The p-value from the mutual information test of independence using the z estimator (`P-value.MIz`) is given for each selected feature.
#' @examples
#' ## Generate a toy dataset: "data"
#' ## Features 1 and 3 are associated with Y, while feature 2 is irrelevant.
#' ## The outcome variable Y must be discrete and in the LAST column. Features must be discrete.
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
#' @importFrom EntropyEstimation Entropy.z MI.z
#' @importFrom entropy entropy.plugin mi.plugin
#' @importFrom stats pchisq qnorm
#'
#' @export



CASMI.selectFeatures <-function(data, alpha.filter=0.1, alpha=0.05){
  # outcome must be in the last column


  # Step 1: return dependent features, step1Index[]

  n=nrow(data)

  step1Index=vector()
  count=0
  k2=length(table(data[length(data)])) #outcome categories
  for(fi in 1:(length(data)-1)){
    k1=length(table(data[fi])) #feature categories
    p_value=MIz.test(data[[fi]],data[[length(data)]], k1, k2)
    if(p_value<=alpha.filter){
      count=count+1
      step1Index[count]=fi
    }
  }

  indexCASMI=vector()
  kappaStarCASMI = vector()

  if(count==0){

    warning("No associated variables are found toward the outcome. First, please ensure that the inputted data is in the data frame format (data.frame()) as the format may affect the result. Second, please ensure that the outcome variable is in the last column. Third, please ensure that a reasonable alpha.filter is set, or use the default value.")

    indexCASMI=NULL
  } else{
    # Step 2: select features by joint SMI with Hz estimator
    # select the best feature
    maxKappaStar=0
    maxIndex=0
    for(index in step1Index){
      #valueKappaStar=kappa.star(data[[index]], data[[length(data)]])
      feature=data[[index]]; outcome=data[[length(data)]];

      smi = MI.z(table(feature, outcome))/Entropy.z(table(outcome))
      #valueKappaStar=smi*(1-length(which(table(feature)==1))/n)
      valueKappaStar=smi*(1-sum(table(feature)==1L)/n)

      if(valueKappaStar > maxKappaStar){
        maxKappaStar=valueKappaStar
        maxIndex=index
      }
    }
    indexCASMI=c(indexCASMI, maxIndex)
    kappaStarCASMI = c(kappaStarCASMI, maxKappaStar)

    if(length(step1Index)==1){

      return(generateResults(data, indexCASMI, kappaStarCASMI, alpha))

    }

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

        smi = MI.z(ftable)/Entropy.z(table(outcome))
        #valueKappaStar=smi*(1-length(which(ftableOfFeatures==1))/n)
        valueKappaStar=smi*(1-sum(ftableOfFeatures==1L)/n)

        if(valueKappaStar > maxKappaStar){
          maxKappaStar=valueKappaStar
          maxIndex=index
        }
      }

      if(maxKappaStar>maxmaxKappaStar+10^-14){ # +10^-14 is to solve the problem of precision
        indexCASMI=c(indexCASMI, maxIndex)
        kappaStarCASMI = c(kappaStarCASMI, maxKappaStar)
      }
    }

    list = generateResults(data, indexCASMI, kappaStarCASMI, alpha)

    return(list)
  }
}


