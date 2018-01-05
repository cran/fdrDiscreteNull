\encoding{UTF-8}
\name{GeneralizedFDREstimators}
\alias{GeneralizedFDREstimators}
\title{Adaptive false discovery rate procedure using generalized estimator of the null proportion.}
\description{
Implement false discovery rate procedures, inlcuding the Adaptive Benjamini-Hochberg procedure, the Adaptive Benjamini-Hochberg-Heyse procedure, using the generalized estimator of the proportion of true nulls, for discrete p-values distributions.
}
\usage{
GeneralizedFDREstimators(data=NULL,
 Test=c("Binomial Test", "Fisher's Exact Test"),
 FET_via = c("PulledMarginals","IndividualMarginals"),
 OneSide = NULL,FDRlevel=NULL,TuningRange = c(0.5,100)) 
}
\arguments{
  \item{data}{Data to be analyzed in the form of a matrix for which observations for a single entity are in a row. 
                 Format of data will be checked by this function automatically and the functions stops execution if the format is wrong.}
  \item{Test}{The type of test to be conducted. It should be exactly one entry from the string 
                 c("Binomial Test", "Fisher's Exact Test"). Currently no other type of test is supported by the package.}
  \item{FET_via}{When the type of test is the Fisher's exact test, how the marginal counts are formed should be specified to be
                 exactly one entry from the string c("PulledMarginals", "IndividualMarginals"). When "PulledMarginals" is used, the data matrix 
                 should have only two clumns, each row of which contains the observed counts for the two binomial distributions, 
                 whereas when "IndividualMarginals" is used the data matrix should have four columns, each row of which has the first and third 
                 entries as the observed count and total number of trials of one binomial distribution, and the second and fourth entries as the
                 observed count and total number of trials of the other binomial distribution. For other types of test, this argument need not to
                 be specified.}
 \item{OneSide}{Specify if one-sided p-value is to be computed from the test. If "OneSide= NULL", then two-sided p-value
  will be computed; if `OneSide="Left"', then the p-value is computed using the left tail of the CDF of the test statistics;
  if `OneSide="Right" ', then the p-value is computed using the right tail of the CDF of the test statistics.}
 \item{FDRlevel}{The nominal false discovery rate (FDR) no larger than which the method to be applied is to have.}
 \item{TuningRange}{A vector of 2 scalars (a,b). Let rho be the maximum of the minimum of each support whose minimum is smaller than 1. If rho is smaller 0.5, then the smallest guiding value is set as a times (0.5-rho) and the biggest guiding value as 0.5, and b determines the number of equally spaced guiding values. If rho is at least 0.5, then all guiding values are set to be rho and b=1.}
}

\value{
It returns the following lists:
 \item{BH}{Restuls obtained by the Benjamini-Hochberg (BH) procedure.}
 \item{aBH}{Results obtained by the adaptive BH procedure using the generalized estimator of the proportion.}
 \item{BHH}{Results obtained by the Benjamini-Hochberg-Heyse (BHH) procedure.}
 \item{aBHH}{Results obtained by the adaptive BHH (aBHH) procedure using the generalized estimator of the proportion.}
 
Each of the above contains:
 \item{pi0Est}{The estimated proprtion of true nulls, where for the BH
              procedure, it is set to be 1.}
 \item{Threshold}{The threshold below which p-values and their associated hypotheses are rejected.}
 \item{NumberOfDiscoveries}{The number of rejections.}
 \item{IndicesOfDiscoveries}{The row indices of the data matrix for the rejections.}

It also returns the following:
  \item{pvalues}{Vector of p-values of the individual tests without grouping.}
  \item{pvalSupp}{It is a list. For binomial test, each entry of the list is a vector, whose first element is the mean of the p-value under the null,
                  second element the p-value itself, and the rest the values at the support of the discrete cdf of the p-value without grouping; for
                  Fisher's exact test, the structure of the list is the same except that in the vector the element denoting the p-value itself 
                  is removed.}
                  
Finally, it also returns randomized p-values (as "RndPval") and results (as "SARP") of the procedure in  Habiger (2015) that is exactly the procedure of Storey et al. (2004) applied to the randomized p-values, and mid p-values (as "MidPval") and "aBHmidP" as the adaptive BH procedure of Benjamini and Hochberg (1995) applied to these mid p-values together with the estimated proportion of true null hypotheses obtained by Storey's estimator in Storey et al. (2004) applied to these mid p-values.
}
\references{
 Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a
 practical and powerful approach to multiple testing, J. R. Statist. Soc. Ser. B
 57(1): 289-300.
 
  Chen, X. and Doerge, R. (2017). A weighted FDR procedure under discrete and heterogeneous null distributions,
 \url{https://arxiv.org/abs/1502.00973v4}.
 
  Chen, X., Doerge, R. and Heyse, J. F. (2017). Multiple testing with discrete data: proportion of true null hypotheses and two adaptive {FDR} procedures, \url{https://arxiv.org/abs/1410.4274v2}.

Habiger, J. D. (2015). Multiple test functions and adjusted p-values for test statistics with
discrete distributions, J. Stat. Plan. Inference 167: 1-13.
 
Heyse, J. F. (2011). A false discovery rate procedure for categorical data, in M. Bhattacharjee,
S. K. Dhar and S. Subramanian (eds), Recent Advances in Biostatistics: False Discovery
Rates, Survival Analysis, and Related Topics, chapter 3.

Lister, R., O'Malley, R., Tonti-Filippini, J., Gregory, B. D., Berry, Charles C. Millar,
 A. H. and Ecker, J. R. (2008). Highly integrated single-base resolution maps of the
 epigenome in arabidopsis, Cell 133(3): 523-536.
 
Storey, J. D., Taylor, J. E. and Siegmund, D. (2004). Strong control, conservative point estimation
in simultaneous conservative consistency of false discover rates: a unied approach, J.
R. Statist. Soc. Ser. B 66(1): 187-205.
}
\seealso{
\code{\link{GeneralizedEstimatorsGrouped}}
}
\examples{
library(qvalue)
library(fdrDiscreteNull)
data(listerdata)
ResG = GeneralizedFDREstimators(listerdata[1:100,],
  Test= "Fisher's Exact Test", FET_via = "PulledMarginals",
  OneSide = NULL,FDRlevel=0.05,TuningRange = c(0.5,20))
}
\keyword{GeneralizedEstimators}
