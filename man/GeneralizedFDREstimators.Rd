\encoding{UTF-8}
\name{GeneralizedFDREstimators}
\alias{GeneralizedFDREstimators}
\title{Multiple hypothesis testing by generalized estimators.}
\description{
Implement multiple testing by generalized estimators, i.e., the generalized estimator of the proportion of true nulls
and the generalized FDR estimator, for discrete p-values distributions.
}
\usage{
GeneralizedFDREstimators(data=NULL, 
       Test=c("Binomial Test", "Fisher's Exact Test",
       "Exact Negative Binomial Test"),
       FET_via = c("PulledMarginals", "IndividualMarginals"),
       ss=NULL,grp.ids = NULL,CommonDisp = NULL,
       nOptimize = 5,FDRlevel=0.05,lambda=0.5,epsilon=1) 
}
\arguments{
  \item{data}{Data to be analyzed in the form of a matrix for which observations for a single entity are in a row. 
                 Format of data will be checked by this function automatically and the functions stops execution if format is wrong.}
  \item{Test}{The type of test to be conducted, which should be exactly one entry from the string 
                 `c("Binomial Test", "Fisher's Exact Test","Exact Negative Binomial Test")'. Currently no other type of test is
                 supported by the package. When the test is "Exact Negative Binomial Test" (ENT), the counts will be 
                 automatically normalized by R package edgeR. Futher, when the test is ENT, the user has the option to specify
                 if common dispersion or tagwise dispersions will be used.}
  \item{FET_via}{When the type of test is the Fisher's exact test, how the marginal counts are formed should be specified to be
                 exactly one entry from the string 
                 c("PulledMarginals",
                 "IndividualMarginals"). When "PulledMarginals" is used, the data matrix 
                 should have only two clumns each row of which contains the observed counts for the two binomial distributions, 
                 whereas when "IndividualMarginals" is used the data matrix should have four columns each row of which has the first and third entries
                 as the observed count and total number of trials of one binomial distribution, and the second and fourth entries as the observed 
                 count and total number of trials of the other binomial distribution. For other types of test, this argument need not to be specified.}
  \item{ss}{Number of replicated for each condition. It is only needed when the test ENT is used, and it should be at least 3. 
            When it is needed, in a row of the data matrix, the first ss replicates should be for one condition, and the rest should be for the other.}
  \item{grp.ids}{When the exact negative binomial test (ENT) is used, the group id's should be specified. Since the ENT needs at least
                  three replicates for each of the two treatments or conditions, the group id's should be a row vector of
                  at least length six; e.g., grp.ids = rep(1:2,each = 3). For other types of test, this argument need not to be specified.}
  \item{CommonDisp}{This argument should be exactly one entry from the string c("Yes","No"). It is needed only when the ENT is used, and specifies if
                  common dispersion will be used by the ENTs.}
  \item{nOptimize}{The number of iteration used by R package edgeR to estimate the common dispersion or tagwise dispersions. This argument
                          is needed only when the ENT is used, and by default it is set to be 5.}
  \item{FDRlevel}{The nominal false discovery rate (FDR) no larger than which the method to be applied is to have.}
  \item{lambda}{The first tuning parameter of the generalized estimator of the proportion of true nulls. It is used to assess if p-values greater
                   than this value should be considered as from the null hypothesis. By default, it is set to be 0.5.}
  \item{epsilon}{The second tuning parameter of the generalized estimator of the proportion of true nulls. It is used to adjust the jumps of the
                    discrete p-value cdf's when estimating the proportion of true nulls, so as to produce less upwardly biased estimate. 
                    By default, it is set to be 1.}
}

\value{
It returns results on multiple testing, each of which is a list:
 \item{GeneralizedEstimator}{Restuls obtained by the generalized FDR procedure.}
 \item{StoreyEstimator}{Results obtained by Storey's FDR procedure.}
 \item{BHProcedure}{Results obtained by Benjamini-Hochberg (BH) procedure.}

Each of the above contains:
 \item{pi0Est}{The estimated proprtion of true nulls, where for the BH
              procedure, this proportion of estimated to be 1.}
 \item{Threshold}{The threshold below which p-values and their associated hypotheses are rejected.}
 \item{NumberOfDiscoveries}{The number of rejections.}
 \item{IndicesOfDiscoveries}{The row indices of the data matrix for the rejections.}
 \item{FDPEst}{The realized false discovery proportion, where for the BH procedure it is set to be the speficied
              nominal FDR level.}

The returned list "GeneralizedEstimator" contains information on the discrete cdf of the p-values:
  \item{pvalues}{Vector of two-sided p-values of the individual tests without grouping.}
  \item{pvalSupp}{It is a list. For binomial test and exact negative binomial test, each entry of the list is a vector related to a p-value, 
  whose first element is the mean of the p-value under the null,
                  second element the p-value itsefl, and the rest the values at the support of the discrete cdf of the p-value without grouping; for
                  Fisher's exact test, the structure of the list is the same except that in the vector the element denoting the p-value itself 
                  is removed.}
  \item{Deltas}{Vector of the final adjustment involved in the generalized estimator of the proportion of true nulls for each discrete cdf of
                        the pvalues without grouping.}
                      
For the exact negative binomial test, it returns:
 \item{DispEst}{The estimated common dispersion or the median of the tagwise dispersion when the ENT is used.}
}
\references{
 Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a
 practical and powerful approach to multiple testing, J. R. Statist. Soc. Ser. B
 57(1): 289-300.
 
 Chen, X. and Doerge, R. (2014a). Generalized estimators formultiple testing: pro559
 portion of true nulls and false discovery rate, \url{http://arxiv.org/abs/1410.4274}.
 
 Chen, X. and Doerge, R. (2014b). False discovery rate control under discrete and hetero-
 geneous null distributions, Submitted.
 
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
\code{\link{GeneralizedEstimatorsGrouped}}
}
\examples{
# due to the "few seconds" policy of CRAN, the follow example
# has been tested successfully before submission to CRAN
\dontrun{
# load library (please install them if needed)
library(fdrDiscreteNull)
library(MCMCpack)
data(arabidopsisM)  # load data
tmp = arabidopsisM
# filter: take counts between 1 and 25
 tmpA = rep(0,2)
for (j in 1:nrow(tmp)) {
   if (sum(tmp[j,])> 0 & max(tmp[j,])<= 25 )
   {
    tmpA = rbind(tmpA,tmp[j,])
   }
  }

tmpA = tmpA[-1,]

m = dim(tmpA)[1]   # check dimensions and get # of hypothesis to test

# apply method
Arab_Lister_BT = GeneralizedFDREstimators(data=tmpA,
  Test= "Fisher's Exact Test",FET_via = "PulledMarginals",
  FDRlevel=0.05, lambda=0.5, epsilon=1)

names(Arab_Lister_BT)
names(Arab_Lister_BT$Generalized_Estimator)
ResultGenEst = Arab_Lister_BT$Generalized_Estimator
ResultGenEst$NumberOfDiscoveries                                               
}
}
\keyword{GeneralizedEstimators}