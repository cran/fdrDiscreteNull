\encoding{UTF-8}
\name{GeneralizedEstimatorsGrouped}
\alias{GeneralizedEstimatorsGrouped}
\title{Weighted multiple hypothesis testing under discrete and heterogeneous null distributions.}
\description{
Implement weighted multiple testing using the generalized estimator of the proportion of true null hypotheses,
where groups are formed by a new divergence on bounded cadlag functions, weights obtained from the groups
using the generalized estimator of the proportion of true null hypotheses, p-values weighted, and multiple testing conducted.
}
\usage{
GeneralizedEstimatorsGrouped(data_in = NULL,
 grpby = c("quantileOfRowTotal","kmeans","divergence"),ngrp_in = NULL,
 GroupMergeSize = 150, minGroupSize = 50,
 test_in = NULL,FET_via_in = NULL,OneSide_in = NULL, FDRlevel_in = NULL,
 eNetSize = NULL, unif_tol= 10^-3, Tunings = c(0.5,100)) 
}
\arguments{
  \item{data_in}{Data to be analyzed in the form of a matrix for which observations for a single entity are in a row. 
                 Format of data will be checked by this function automatically and the functions stops execution if the format is wrong.}
   \item{grpby}{The method to be used to form the groups. It should be exactly one entry from the string c("quantileOfRowTotal","kmeans","divergence").  Grouping by "quantileOfRowTotal" is a good choice as demonstrated by simulation stuides and it is very fast.}
    \item{ngrp_in}{The number of groups to be formed from the orginal data. It refers to the number of groups that the rows of the data matrix
                 will be formed, and also to the number of groups that the discrete null distributions and their associated p-values will be formed.}
\item{minGroupSize}{When the grouping method is "divergence", the default minimal group size "minGroupSize" is 50, which means at least 50 hypotheses are needed to form a group.}
\item{GroupMergeSize}{When the grouping method is "divergence", the last group may merge what is left from forming "ngrp_in" groups or from forming the first "group_in - 1" groups, if the number of hypotheses left is less than or equal to "GroupMergeSize". This can speed up grouping by divergence.}
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
  \item{eNetSize}{The argument is needed only when both the arguments ``divergence'' and `RefDivergence="No"' are used. It specifies the size of the metric balls
                  to be used to partition the set of discrete cdf's to form the groups.}
  \item{unif_tol}{The argument is needed only when the argument ``divergence'' is used. It specifies the tolerance on the infinity norm
                  under which a discrete cdf of a p-value will be considered approximately uniform on [0,1]. 
                  By default, it is set to be 0.001.}
  \item{Tunings}{A vector of 2 scalars (a,b). Let rho be the maximum of the minimum of each support whose minimum is smaller than 1. If rho is smaller 0.5, then the smallest guiding value is set as a times (0.5-rho) and the biggest guiding value as 0.5, and b determines the number of equally spaced guiding values. If rho is at least 0.5, then all guiding values are set to be rho and b=1.}
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
 Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a
 practical and powerful approach to multiple testing, J. R. Statist. Soc. Ser. B
 57(1): 289-300.
 
 Chen, X. and Doerge, R. (2017). A weighted FDR procedure under discrete and heterogeneous null distributions,
 \url{https://arxiv.org/abs/1502.00973v4}.

Chen, X., Doerge, R. and Heyse, J. F. (2017). Multiple testing with discrete data: proportion of true null hypotheses and two adaptive {FDR} procedures, \url{https://arxiv.org/abs/1410.4274v2}.
 
 Heyse, J. F. (2011). A false discovery rate procedure for categorical data, in M. Bhattacharjee,
S. K. Dhar and S. Subramanian (eds), Recent Advances in Biostatistics: False Discovery
Rates, Survival Analysis, and Related Topics, chapter 3.

 Lister, R., O'Malley, R., Tonti-Filippini, J., Gregory, B. D., Berry, Charles C. Millar,
 A. H. and Ecker, J. R. (2008). Highly integrated single-base resolution maps of the
 epigenome in arabidopsis, Cell 133(3): 523-536.
}
\seealso{
\code{\link{GeneralizedFDREstimators}}
}
\examples{
library(fdrDiscreteNull)
library(qvalue)
data(listerdata)
ResTmp = GeneralizedEstimatorsGrouped(listerdata[1:500,], 
  grpby= "quantileOfRowTotal", ngrp_in= 3,GroupMergeSize = 150, minGroupSize = 50,
  test_in= "Fisher's Exact Test", FET_via_in = "PulledMarginals",OneSide_in = NULL, 
  FDRlevel= 0.05,Tunings = c(0.5,20))
}
\keyword{GeneralizedEstimatorsGrouped}
