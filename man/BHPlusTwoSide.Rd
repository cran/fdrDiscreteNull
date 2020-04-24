\encoding{UTF-8}
\name{BHPlusTwoSide}
\alias{BHPlusTwoSide}
\title{FDR control for multiple testing based on p-values with cadlag distribution functions.}
\description{
Implement the BH+ procedures of Chen, X. (2019) for FDR control for multiple testing based on p-values whose distributions are cadlag, i.e., right-continuous with left-limits. This includes conventional
p-values and mid p-values. Currently, the methods are implemented for two-sided p-values.
}
\usage{
BHPlusTwoSide(data=NULL, Test=c("Binomial Test", "Fisher's Exact Test"),
                                  FET_via = c("PulledMarginals","IndividualMarginals"),
                                  FDRlevel=NULL, espilon=NULL) 
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
 \item{FDRlevel}{The nominal false discovery rate (FDR) no larger than which the method to be applied is to have.}
 \item{espilon}{A scalar used to determine the guiding value for the estimator of the proportion of true null hypotheses. It is usually set to be 0.01.}
}

\value{
It returns the following lists:
 \item{BH}{Restuls obtained by the Benjamini-Hochberg (BH) procedure when applied to conventional p-values.}
 \item{BHplus}{Results obtained by the BH+ procedure when applied to conventional p-values.}
 \item{aBHplus}{Results obtained by the adaptive BH+ procedure when applied to conventional p-values.}
 \item{MidpBHplus}{Results obtained by the BH+ procedure when applied to mid p-values.}
 \item{aMidpBHplus}{Results obtained by the adaptive BH+ procedure when applied to mid p-values.}
 
Each of the above contains:
 \item{pi0Est}{The estimated proprtion of true nulls, where for non-adaptive
              procedure, it is set to be 1.}
 \item{Threshold}{The threshold below which p-values and their associated hypotheses are rejected.}
 \item{IndicesOfDiscoveries}{The row indices of the data matrix whose corresponding hypotheses are rejected.}

For an adaptive procedure, each of the above lists contains:
 \item{Tuning}{The guiding value for the estimator of the proportion.}
 \item{TuningVon}{If TuningVon=0, the guiding value is chosen by theory, meaning that the adaptive procedure is conservative; if TuningVon=1, the guiding value
 is chosen approximately, meaning that the adaptive procedure may not be conservative. A user may get the "Warning message:
In max(DevCount) : no non-missing arguments to max; returning -Inf", which can be safely ignored.}

It also returns the following:
  \item{pval}{Vector of conventional p-values, whose indices matches those of the hypotheses.}
  \item{pvalSupp}{It is a list, each of whose element is the support of a conventional p-value. The indice of pvalSupp match those of pval.}
   \item{Midpvals}{Vector of mid p-values, whose indices matches those of the hypotheses.}
  \item{pvalMidSupp}{It is a list, each of whose element is the support of a mid p-value. The indice of pvalSupp match those of Midpvals.}
                                  
}
\references{
 Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a
 practical and powerful approach to multiple testing. J. R. Statist. Soc. Ser. B
 57(1): 289-300.

 Hwang, J. T. G. and Yang, M.-C. (2001). An optimality theory for mid pcvalues in 2 x 2
contingency tables. Statistica Sinica 11(3): 807-826.

 Chen, X. (2019). False discovery rate control for multiple testing based on discrete p-values. \url{https://arxiv.org/abs/1803.06040}; Biometrical Journal (in press).

 Gilbert, P. B. (2005). A modified false discovery rate multiple-comparisons procedure for discrete
data, applied to human immunodeficiency virus genetics, J. R. Statist. Soc. Ser. C 54(1): 143-158.
}

\examples{
library(fdrDiscreteNull)
data(hivdata)
m = dim(hivdata)[1]
hivdata = cbind(hivdata,rep(73,m),rep(73,m))
hivdataA = hivdata[rowSums(hivdata[,1:2])>=5,]
ResTmp = BHPlusTwoSide(data = hivdataA, Test = "Fisher's Exact Test",
 FET_via = "IndividualMarginals",FDRlevel = 0.05)
}
\keyword{BHPlusTwoSide}
