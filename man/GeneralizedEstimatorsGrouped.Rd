\encoding{UTF-8}
\name{GeneralizedEstimatorsGrouped}
\alias{GeneralizedEstimatorsGrouped}
\title{Weighted multiple hypothesis testing under discrete and heterogeneous null distributions.}
\description{
Implement weighted multiple testing procedure of Chen, X., Doerge, R. and Sanat, S. K. (2019) for independent p-values whose null distributions are super-uniform but not necessarily identical or continuous,
where groups are formed by the infinity norm for functions, p-values weighted by data-adaptive weights, and multiple testing conducted. The weights are not
constructed using plug-in null proportion estimators. For multiple testing based on p-values of Binomial tests or Fisher's exact tests, grouping using quantiles
of observed counts is recommended both for fast implementation and excellent power performance of the weighted FDR procedure.
}
\usage{
GeneralizedEstimatorsGrouped(data_in = NULL,
 grpby = c("quantileOfRowTotal","kmeans","InfNorm"),ngrp_in = NULL,
 GroupMergeSize = 50,test_in = NULL,FET_via_in = NULL,OneSide_in = NULL,
 FDRlevel_in = NULL, eNetSize = NULL, unif_tol= 10^-3,
 TuningPar = c(0.5,100),lambda = 0.5) 
}
\arguments{
  \item{data_in}{Data to be analyzed in the form of a matrix for which observations for a single entity are in a row. 
                 Format of data will be checked by this function automatically and the functions stops execution if the format is wrong.}
   \item{grpby}{The method to be used to form the groups. It should be exactly one entry from the string c("quantileOfRowTotal","kmeans","InfNorm").  Grouping by "quantileOfRowTotal" is a good choice as demonstrated by simulation stuides and it is very fast.}
    \item{ngrp_in}{The number of groups to be formed from the orginal data. It refers to the number of groups that the rows of the data matrix
                 will be formed, and also to the number of groups that the discrete null distributions and their associated p-values will be formed.}
\item{GroupMergeSize}{When the grouping method is "InfNorm", the default minimal group size is "GroupMergeSize", which is 50 by default.}
  \item{test_in}{The type of test to be conducted. It should be exactly one entry from the string 
                 c("Binomial Test", "Fisher's Exact Test"). Currently no other type of test is
                 supported by the package.}
  \item{FET_via_in}{When the type of test is the Fisher's exact test, how the marginal counts are formed should be specified to be
                    exactly one entry from the string 'c("PulledMarginals", "IndividualMarginals")'. When "PulledMarginals" is used, the data matrix 
                 should have only two clumns, each row of which contains the observed counts for the two binomial distributions, 
                 whereas when "IndividualMarginals" is used the data matrix should have four columns, each row of which has the first and third entries
                 as the observed count and total number of trials of one binomial distribution, and the second and fourth entries as the observed 
                 count and total number of trials of the other binomial distribution. For other types of test, this argument need not to be specified.}
  \item{OneSide_in}{Specify if one-sided p-value is to be computed from the test. If "OneSide_in= NULL", then two-sided p-value
  will be computed; if `OneSide_in="Left"', then the p-value is computed using the left tail of the CDF of the test statistics;
  if `OneSide_in="Right" ', then the p-value is computed using the right tail of the CDF of the test statistics.}
    \item{FDRlevel_in}{The nominal false discovery rate (FDR) no larger than which the method to be applied is to have.}
  \item{eNetSize}{The argument is needed only when the argument ``InfNorm'' is used. It specifies the size of the metric balls
                  to be used to partition the set of discrete cdf's to form the groups.}
  \item{unif_tol}{The argument is needed only when the argument ``InfNorm'' is used. It specifies the tolerance on the infinity norm
                  under which a discrete cdf of a p-value will be considered approximately uniform on [0,1]. 
                  By default, it is set to be 0.001.}
  \item{TuningPar}{A vector of 2 scalars (a,b), used to implement the generalized proportion estimator. Let rho be the maximum of the minimum of each support whose minimum is smaller than 1. If rho is smaller 0.5, then the smallest guiding value is set as a times (0.5-rho) and the biggest guiding value as 0.5, and b determines the number of equally spaced guiding values. If rho is at least 0.5, then all guiding values are set to be rho and b=1.}
  \item{lambda}{A scalor in (0,1) that is used as a tuning parameter to construct data-adaptive weights. By default, it is set to be 0.5.}
}

\value{
It returns estimated proportion of true nulls:
\item{pi0estAll}{Estimated proportion of true nulls.}
The above quantity is a vector and contains the following:
   \item{pi0E_GE}{Estimated proportion of true nulls, obtained by the generalized estimator.}
   \item{pi0E_gGE}{Estimated proportion of true nulls, obtained by grouping and weighting and the generalized estimator.}
   \item{pi0Est_gp*}{Estimated proportion of true nulls for each group by the generalized estimator, where * is a group number.}

It returns the results on multiple testing that are returned by \code{\link{GeneralizedFDREstimators}}, plus the following list:
\item{wFDR}{Results from the weighted false discovery rate procedure; these results are stored using the same list structure as multiple testing
results returned by \cr \code{\link{GeneralizedFDREstimators}}.}
}
\references{
 Chen, X., Doerge, R. and Sanat, S. K. (2020). A weighted FDR procedure under discrete and heterogeneous null distributions.
 \url{https://arxiv.org/abs/1502.00973v5}; Biometrical Journal (in press).

 Lister, R., O'Malley, R., Tonti-Filippini, J., Gregory, B. D., Berry, Charles C. Millar,
 A. H. and Ecker, J. R. (2008). Highly integrated single-base resolution maps of the
 epigenome in arabidopsis. Cell 133(3): 523-536.
}
\seealso{
\code{\link{GeneralizedFDREstimators}}
}
\examples{
library(fdrDiscreteNull)
library(qvalue)
data(listerdata)
ResTmp = GeneralizedEstimatorsGrouped(listerdata[1:500,], 
  grpby= "quantileOfRowTotal", ngrp_in= 3, GroupMergeSize = 50,
  test_in= "Fisher's Exact Test", FET_via_in = "PulledMarginals",OneSide_in = NULL, 
  FDRlevel= 0.05,TuningPar = c(0.5,20),lambda = 0.5)
}
\keyword{GeneralizedEstimatorsGrouped}
