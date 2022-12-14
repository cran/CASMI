% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeatureSelection.R
\name{CASMI.selectFeatures}
\alias{CASMI.selectFeatures}
\title{\pkg{CASMI}-Based Feature Selection}
\usage{
CASMI.selectFeatures(data, alpha.filter = 0.1, alpha = 0.05)
}
\arguments{
\item{data}{data frame (features as columns and observations as rows). It requires at least one feature and only one outcome. The features must be discrete. The outcome variable (Y) must be in the last column.}

\item{alpha.filter}{level of significance for the mutual information test of independence in step 1 of the features selection. The smaller the alpha.filter, the fewer the features sent to step 2 (<doi:10.3390/e21121179>). By default, `alpha.filter = 0.1`.}

\item{alpha}{level of significance for the confidence intervals in final results. By default, `alpha = 0.05`.}
}
\value{
`CASMI.selectFeatures()` returns selected features and relevant information, including the estimated Kappa* for all selected features (`$KappaStar`) and the corresponding confidence interval (`$KappaStarCI`). The selected features are ranked. The Standardized Mutual Information using the z estimator (`SMIz`) and the corresponding confidence interval (`CI.SMIz.Lower` and `CI.SMIz.Upper`) are given for each selected feature. The p-value from the mutual information test of independence using the z estimator (`P-value.MIz`) is given for each selected feature.
}
\description{
Selects the most relevant features toward an outcome. It automatically learns the number of features to be selected, and the selected features are ranked. The method automatically handles the feature redundancy issue. (synonyms of "feature": "variable" "factor" "attribute") \cr
For more information, please refer to the corresponding publication: Shi, J., Zhang, J. and Ge, Y. (2019), "\pkg{CASMI}—An Entropic Feature Selection Method in Turing’s Perspective" <doi:10.3390/e21121179>
}
\examples{
## Generate a toy dataset: "data"
## Features 1 and 3 are associated with Y, while feature 2 is irrelevant.
## The outcome variable Y must be in the last column. Features must be discrete.
n=10000
x1=rbinom(n,3,0.5)+0.2
x2=rbinom(n,2,0.8)+0.5
x3=rbinom(n,5,0.3)
error=round(runif(n,min=-1,max=1))
y=x1+x3+error
data=data.frame(cbind(x1,x2,x3,y))
colnames(data) = c("feature1", "feature2", "feature3", "Y")

## Select features and relevant results from the toy dataset "data"
CASMI.selectFeatures(data)

}
