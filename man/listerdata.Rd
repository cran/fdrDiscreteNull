\encoding{UTF-8}
\name{listerdata}
\alias{listerdata}
\title{Methylation data for Arabidopsis thaliana}
\description{
This data set has been analyzed and provided by the listed reference. There are around
22000 cytosines, each of which is under two conditions. For each cytosine under
each condition, there is only one replicate. The discrete count for each replicate
can be modelled by binomial distribution, and Fisher'e exact test can be applied 
to assess if a cytosine is differentially methylated. The filtered data "listerdata.RData" contains
cytosines whose total counts for both lines are greater than 5 and whose count for each line does not exceed 25.
}
\references{ 
 Lister, R., O'Malley, R., Tonti-Filippini, J., Gregory, B. D., Berry, Charles C. Millar,
 A. H. and Ecker, J. R. (2008). Highly integrated single-base resolution maps of the
 epigenome in arabidopsis, Cell 133(3): 523-536.
}