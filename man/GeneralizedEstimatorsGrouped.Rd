\encoding{UTF-8}
\name{GeneralizedEstimatorsGrouped}
\alias{GeneralizedEstimatorsGrouped}
\title{Weighted multiple hypothesis testing by generalized estimators.}
\description{
Implement weighted multiple testing by generalized estimators, where groups are formed
by a new divergence on cadlag functions, weights obtained from the groups
using the generalized estimator of the proportion of true null
hypotheses, p-values weighted, and multiple testing conducted.
}
\usage{
GeneralizedEstimatorsGrouped(data_in = NULL, test_in = NULL, FET_via_in = NULL,
         ss=NULL, gpids_in = NULL, cdisp_in = NULL, nopsinestdisp_in = 5,
         grpby = c("quantileOfRowTotal", "kmeans", "divergence"),
         ngrp_in = NULL, RefDivergence = NULL, eNetSize = NULL, unif_tol= 10^-3,
         FDRlevel_in = 0.05, lambda_in = 0.5, epsilon_in = 1)  
}
\arguments{
  \item{data_in}{Data to be analyzed in the form of a matrix for which observations for a single entity are in a row. 
                 Format of data will be checked by this function automatically and the functions stops execution if the format is wrong.}
  \item{test_in}{The type of test to be conducted. It should be exactly one entry from the string 
                 c("Binomial Test", "Fisher's Exact Test","Exact Negative Binomial Test"). Currently no other type of test is
                 supported by the package. When the test is the Exact Negative Binomial Test (ENT), normalization of counts will
                 be done automatically by R package edgeR. Futher, when the test is ENT, the user has the option to specify
                 if common dispersion or tagwise dispersions will be used.}
  \item{FET_via_in}{When the type of test is the Fisher's exact test, how the marginal counts are formed should be specified to be
                    exactly one entry from the string 
                    c("PulledMarginals",
                    "IndividualMarginals"). When "PulledMarginals" is used, the data matrix 
                 should have only two clumns each row of which contains the observed counts for the two binomial distributions, 
                 whereas when "IndividualMarginals" is used the data matrix should have four columns each row of which has the first and third entries
                 as the observed count and total number of trials of one binomial distribution, and the second and fourth entries as the observed 
                 count and total number of trials of the other binomial distribution. For other types of test, this argument need not to be specified.}
  \item{ss}{Number of replicates for each condition. It is only needed when the test ENT is used, and it should be at least 3. 
            When it is needed, in a row of the data matrix, the first ss replicates should be for one condition, and the rest should be for the other.}
  \item{gpids_in}{When the exact negative binomial test (ENT) is used, the group id's should be specified. Since the ENT needs at least
                  three replicates for each of the two treatments or conditions, the ground id's should be a row vector of
                  at least length six; e.g., gpids_in = rep(1:2,each = 3). For other types of test, this argument need not to be specified.}
  \item{cdisp_in}{This argument should be exactly one entry from the string c("Yes","No"). It is needed only when the ENT is used, and specifies if
                  common dispersion will be used by the ENTs.}
  \item{nopsinestdisp_in}{The number of iteration used by R package edgeR to estimate the common dispersion or tagwise dispersions. This arguments
                          is needed only when the ENT is used, and by default it is set to be 5.}
  \item{grpby}{The method to be used to form the groups. It should be exactly one entry from the string c("quantileOfRowTotal","kmeans","divergence").
               When the grouping method is "divergence", the default minimal group size and the merging size are both 30, which means at least 30 hypotheses
               are needed to implement the grouping strategy. For the three types of test the package supports, grouping by "quantileOfRowTotal" is
               a good choice as demonstrated by simulation stuides and theory.}
  \item{ngrp_in}{The number of groups to be formed from the orginal data. It refers to the number of groups the rows of the data matrix
                 will be formed, and also to the number of groups of the discrete null distributions and their associated p-values will be formed.}
  \item{RefDivergence}{The argument should be exactly one entry from the string c("Yes","No"). When the grouping method is "divergence", 
                       the reference distribution, i.e., the uniform distribution on [0,1], can be used to group the discrete and heterogenenous null
                       p-values distributions with respect to their ``distance'' from this uniform distribution. The argument is needed only when
                       ``divergence'' is used to group the distributions.}
  \item{eNetSize}{The argument is needed only when both the arguments ``divergence'' and `RefDivergence="No"' are used. It specifies the size of the metric balls
                  to be used to partition the set of discrete cdf's to form the groups.}
  \item{unif_tol}{The argument is needed only when both the arguments ``divergence'' and `RefDivergence="Yes"' are used. It specifies the tolerance
                  under which a discrete cdf of a p-value should be considered approximately uniform on [0,1]. 
                  By default, it is set to be 0.001.}
  \item{FDRlevel_in}{The nominal false discovery rate (FDR) no larger than which the method to be applied is to have.}
  \item{lambda_in}{The first tuning parameter of the generalized estimator of the proportion of true nulls. It is used to assess if p-values greater
                   than this value should be considered as from the null hypothesis. By default, it is set to be 0.5.}
  \item{epsilon_in}{The second tuning parameter of the generalized estimator of the proportion of true nulls. It is used to adjust the jumps of the
                    discrete p-value cdf's when estimating the proportion of true nulls, so as to produce less upwardly biased estimate. 
                    By default, it is set to be 1.}
}

\value{
It returns estimated proportion of true nulls:
\item{pi0ests}{Estimated proportion of true nulls by various estimators of this proportion.}
The above quantity is a vector and contains the following:
   \item{pi0E_Gen}{Estimated proportion of true nulls, obtained by the generalized estimator.}
   \item{pi0E_Storey}{Estimated proportion of true nulls, obtained by Storey's estimator.}
   \item{pi0E_gwGen}{Estimated proportion of true nulls, obtained by grouping and weighting, and the generalized estimator.}
   \item{pi0E_gwStorey}{Estimated proportion of true nulls, obtained by grouping and weighting, and Storey's estimator.}
   \item{pi0Est_gp*}{Estimated proportion of true nulls for each group by the generalized estimator, where * is a group number.}
   \item{piEst_st_gp*}{Estimated proportion of true nulls for each group by Storey's estimator, where * is a group number.}

It returns the following results on multiple testing: 
  \item{Gen}{The row indices of the data matrix for the discoveries, i.e., rejected null hypotheses, made by the generalized FDR procedure.}
  \item{Storey}{The row indices of the data matrix for the discoveries made by Storey's FDR procedure.}
  \item{BH}{The row indices of the data matrix for the discoveries made by Benjamini-Hochberg procedure.}
  \item{GWGen}{The row indices of the data matrix for the discoveries made by weighted generalized FDR procedure using grouping.}
  \item{GWStorey}{The row indices of the data matrix for the discoveries made by weighted Strorey's FDR procedure using grouping.
                  This method is just a special case of GWGen.}
  \item{DispEst}{The estimated common dispersion or the median of the tagwise dispersions when the ENT is used.}

It returns the information on the discrete cdf of the p-values:
  \item{pval_ungrped}{Vector of two-sided p-values of the individual tests without grouping.}
  \item{pSupp_ungrped}{It is a list. For binomial test and exact negative binomial test, each entry of the list is a vector related to a p-value, 
  whose first element is the mean of the p-value under the null,
                  second element the p-value itsefl, and the rest the values at the support of the discrete cdf of the p-value without grouping; for
                  Fisher's exact test, the structure of the list is the same except that in the vector the element denoting the p-value itself 
                  is removed.}
  \item{deltas_ungrped}{Vector of the final adjustment involved in the generalized estimator of the proportion of true nulls for each discrete cdf of
                        the pvalues without grouping.}
}
\references{
 Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a
 practical and powerful approach to multiple testing, J. R. Statist. Soc. Ser. B
 57(1): 289-300.
 
 Chen, X. and Doerge, R. (2014). Generalized estimators formultiple testing: pro559
 portion of true nulls and false discovery rate, \url{http://arxiv.org/abs/1410.4274}.
 
 Chen, X. and Doerge, R. (2015). False discovery rate control under discrete and hetero-
 geneous null distributions, \url{http://arxiv.org/abs/1502.00973}.
 
 Di, Y., Schafer, D. W., Cumbie, J. S. and Chang, J. H. (2011). The NBP negative
 binomial model for assessing differential gene expression from RNA-Seq, Statistical
 Applications in Genetics and Molecular Biology 10(1): 24 pages.

 Lister, R., O'Malley, R., Tonti-Filippini, J., Gregory, B. D., Berry, Charles C. Millar,
 A. H. and Ecker, J. R. (2008). Highly integrated single-base resolution maps of the
 epigenome in arabidopsis, Cell 133(3): 523-536.

 Robinson, M. D., McCarthy, D. J. and Smyth, G. K. (2009). edgeR: a Bioconductor
 package for differential expression analysis of digital gene expression data, Bioin-
 formatics 26(1): 139-140.
 
 Storey JD. (2002). A direct approach to false discovery rates. J. R. Statist. Soc. Ser. B 64(3): 479-498.
 
 Storey, J. D., Taylor, J. E. and Siegmund, D. (2004). Strong control, conservative
 point estimation in simultaneous conservative consistency of false discover rates: a
 unified approach, J. R. Statist. Soc. Ser. B 66(1): 187-205.
}
\seealso{
\code{\link{GeneralizedFDREstimators}}
}
\examples{
# due to the "few seconds" policy of CRAN, the follow example
# has been tested successfully before submission to CRAN
\dontrun{
library(fdrDiscreteNull)
data(arabidopsisE)
tmp = arabidopsisE
# extract data with specified range of row total counts
tmpA = tmp[rowSums(tmp[,1:3])> 0 & 
rowSums(tmp[,4:6])> 0 & rowSums(tmp[,1:6])<= 500,]
dim(tmpA)  

# implement the method

ResArab_tagdisp = GeneralizedEstimatorsGrouped(tmpA,
  test_in = "Exact Negative Binomial Test",ss=3,
  gpids_in = rep(1:2,each = 3),cdisp_in = "No",
  grpby = "quantileOfRowTotal",ngrp_in = 10,
  lambda_in = 0.5,epsilon_in = 0.8) 
# view results
names(ResArab_tagdisp) 
sapply(ResArab_tagdisp,length)                                               
}
}
\keyword{GeneralizedEstimatorsGrouped}