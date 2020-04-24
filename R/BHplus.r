#########################################################################
# BH+ procedures by xiongzhi.chen@wsu.edu
####################################################################

BHPlusTwoSide = function(data=NULL, Test=c("Binomial Test", "Fisher's Exact Test"),
                                  FET_via = c("PulledMarginals","IndividualMarginals"),
                                          FDRlevel=NULL,espilon=NULL) 
 { # start of main function
  ###### check data
  if (is.null(data))      stop("^^Please provide data.")
   chktest = identical(Test,"Binomial Test") + identical(Test,"Fisher's Exact Test") 
  if (chktest ==0 )      stop("^^Unsupported type of test.")    
  
  # Fisher's exact test
  if (Test == "Fisher's Exact Test")   {
   if (is.null(FET_via)) stop("^^Please specify type of marginals for Fisher's Exact Test","\n")
   else {
      if (FET_via == "PulledMarginals" & ncol(data) != 2)          stop("Please reformat data into an m-by-2 matrix.")
      if (FET_via == "IndividualMarginals" & ncol(data) != 4)      stop("Please reformat data into an m-by-4 matrix.")
    }
  }
    
  # Binomial test       
  if (Test == "Binomial Test" & ncol(data) != 2)      stop("^^Please reformat data into an m-by-2 matrix.")
  if (is.null(FDRlevel))       stop("^^Please specify FDR level.")
   
  # get m, number of rows of data
  data = as.matrix(data)
  m = nrow(data)
  
  # Section 2: Get p-values and their supports from Binomial Test 
  if (Test == "Binomial Test") {      
       cat("^^Using Binomial Test.","\n")
       # check for zero rows 
       if (any(rowSums(data) == 0)) stop("^^Please remove zero rows from data matrix.","\n")
       
       sizePoi = data[,1] + data[,2]
       marginalsAndsizes = cbind(data[,1], sizePoi)
       
       pValAndSupp = apply(marginalsAndsizes,1,pvalBT)           
     } # end of case 1, if (Test == "Binomial Test")  
  
  # Section 3: Get p-values and their supports from Fisher's Exact Test 
  if (Test == "Fisher's Exact Test") {
        if (any(rowSums(data[,1:2]) == 0)) stop("^^Please remove rows with zero observed total count in data matrix.","\n")
        
        if (FET_via == "PulledMarginals" & ncol(data) == 2) {  
           cat("^^Fisher's Exact Test is based on pulled marginal","\n")

          # get cell counts and marginals for Fisher's exact test
          countveccontrol = data[,1];    countvectreat = data[,2]
          cellcountsmarginals <- getcellcountsandmarginals(countveccontrol,countvectreat)
          }
         
        if (FET_via == "IndividualMarginals" & ncol(data) == 4) {
            cat("^^Fisher's Exact Test is based on individual marginals","\n")
             cellcountsmarginals = getcellcountsandmarginals_DE(data)
           }
        
        simallcellcounts <- cellcountsmarginals[[1]]  # each element of cellcountsmarginals[[1]] is a 2-by-2 matrix
#       simallmarginals <- cellcountsmarginals[[2]]   # each element of cellcountsmarginals[[2]] is a 3 vector
        pValAndSupp = lapply(simallcellcounts,pvalFET)
     }       # end of case 2, if (Test == "Fisher's Exact Test")


  # Section 4: extract p-values and their supports 
   pvals = double(m); pvalsMid = double(m)      
   ## extract supports
   pCDFSupp = vector("list",m); pCDFSuppMidp = vector("list",m)
   for (ix in 1:m)  {
      pvasuppTmp = pValAndSupp[[ix]]
      pvals[[ix]] = pvasuppTmp[[1]][1]
      pvalsMid[[ix]] = pvasuppTmp[[1]][2]
      pCDFSupp[[ix]] = pvasuppTmp[[3]] # conventional pvalue supp
      pCDFSuppMidp[[ix]] = pvasuppTmp[[2]] # values of mid p-val
      
      if (length(pCDFSupp[[ix]]) != length(pCDFSuppMidp[[ix]]))
        cat("^^ Mis-alignment of mid p-value and its support","\n")
   }
  
  
  ## Section 5: Estimate pi0 from p-values and their supports
  pCDFEval = ExtremeCDF(pCDFSupp, pCDFSupp, espilon)
  pCdfMax = pCDFEval[[1]]
  pCdfMin = pCDFEval[[2]] 
  tauA = pCDFEval[[3]]      # conventional p
  tauAVon = pCDFEval[[4]]   # conventional p
  cat("^^^ Estimating pi0 based on conventional p-values ...","\n")
  pi0GConv = GenEstProp(pvals,pCDFSupp,c(tauA,1))  # used default tuning in sim B
   
  #BH  based on conventional pval
  BHOrig <- BH(pvals,FDRlevel)
  # non-adaptive BH plus for conventional pval
  BHPlus = BHplus(pvals, pCdfMax,1,FDRlevel)
  # adaptive BH plus for conventional pval
   aBHPlus = BHplus(pvals, pCdfMax, pi0GConv,FDRlevel) 
  
  ###################### BHplus based on mide p-values
  MidpCDFEval = ExtremeCDF(pCDFSuppMidp, pCDFSupp,espilon)  # set epsi here
  MidpCdfMax = MidpCDFEval[[1]]
  MidpCdfMin = MidpCDFEval[[2]]
  tau = MidpCDFEval[[3]]
  tauVon = MidpCDFEval[[4]]
  cat("^^^ Estimating pi0 based on mid p-values ...","\n")
  pi0GMidp = GenEstProp(pvalsMid,pCDFSupp,c(tau,1))
    cat("^^ pi0estConvP: ",pi0GConv, "; pi0estMIdp: ", pi0GMidp,"\n")
          
  # non-adaptive BHPlus for mid p
  BHPlusMidp = BHplus(pvalsMid, MidpCdfMax,1,FDRlevel)
  # adaptive BHPlus for mid p
  aBHPlusMidp = BHplus(pvalsMid, MidpCdfMax,pi0GMidp,FDRlevel)
  
  # Section 6: put estimates in a dataframe and return it 
  # names for entries: thresh, FDP, number of discoveries, indices of discoveries
  BH = list("pi0Est" = 1, "Threshold" = BHOrig[[2]],"IndicesOfDiscoveries" = BHOrig[[1]][,2])
  BHp = list("pi0Est" = 1, "Threshold" = BHPlus[[2]],"IndicesOfDiscoveries" = BHPlus[[1]][,2])
  aBHp = list("pi0Est" = pi0GConv, "Threshold" =aBHPlus[[2]],"IndicesOfDiscoveries" = aBHPlus[[1]][,2],"Tuning"=tauA,"TuningVon"=tauAVon)
 
 MidpBHp = list("pi0Est" = 1, "Threshold" =BHPlusMidp[[2]],"IndicesOfDiscoveries" = BHPlusMidp[[1]][,2])
 aMidBHp = list("pi0Est" = pi0GMidp, "Threshold" =aBHPlusMidp[[2]],"IndicesOfDiscoveries" = aBHPlusMidp[[1]][,2],"Tuning"=tau, "TuningVon"=tauVon)
 
  ## gathering results  ###########
 AnalysisResults = list("BH"=BH,"BHplus"=BHp,"aBHplus"=aBHp,"MidpBHplus"=MidpBHp,"aMidpBHplus"=aMidBHp,"pval" = pvals,"pvalSupp" = pCDFSupp,"Midpvals"=pvalsMid,"pvalMidSupp"=pCDFSuppMidp)
 
 # rerturn results 
 return(AnalysisResults)  
      
} # end of main function 