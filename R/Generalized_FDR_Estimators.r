#########################################################################
# Generalized Estimator of Proportion of True Nulls
# This file contains the main functions for Generalized FDR Estimators
# Contact: xiongzhi.chen@wsu.edu
####################################################################

### Start of main function  
GeneralizedFDREstimators = function(data=NULL, Test=c("Binomial Test", "Fisher's Exact Test"),
                                  FET_via = c("PulledMarginals","IndividualMarginals"), OneSide = NULL,
                                          FDRlevel=NULL,TuningRange = c(0.5,100)) 
  { # start of main function

  #########################################################################
  # Section 1: check if arguments are correctly provided
  ####################################################################
  
  if (is.null(data))      stop("^^Please provide data.")
  # check test type
   chktest = identical(Test,"Binomial Test") + identical(Test,"Fisher's Exact Test")
   
  if (chktest ==0 )      stop("^^Unsupported type of test.")    
  
  # For Fisher's exact test, two column data is ok
  if (Test == "Fisher's Exact Test")   {
   if (is.null(FET_via)) stop("^^Please specify type of marginals for Fisher's Exact Test","\n")
   else {
        if (FET_via == "PulledMarginals" & ncol(data) != 2)          stop("Please reformat data into an m-by-2 matrix.")
        if (FET_via == "IndividualMarginals" & ncol(data) != 4)      stop("Please reformat data into an m-by-4 matrix.")
        } 
   }
  
  # For Binomial Test, two column data is ok      
  if (Test == "Binomial Test" & ncol(data) != 2)      stop("^^Please reformat data into an m-by-2 matrix.")
  
  if (is.null(FDRlevel))       stop("^^Please specify FDR level.")
   
  # get m, number of rows of data
  data = as.matrix(data)
  m = nrow(data);   n = ncol(data)
  
  ######################################################################### ##########
  # Section 2: Get p-values and their supports from Binomial Test 
  #################################################################### ########## ##########
  if (Test == "Binomial Test") {      
       cat("^^Using Binomial Test.","\n")
       # check for zero rows 
       if (any(rowSums(data) == 0)) stop("^^Please remove zero rows from data matrix.","\n")
        
       sizePoi = data[,1] + data[,2]
       marginalsAndsizes = cbind(data[,1], sizePoi)
       
       ################## two-sided p-values #####################
       if(is.null(OneSide))  {
        cat("^^^Computing two-sided p-value of binomial test and its support","\n")
        } else {
        cat("^^^Computing one-sided", OneSide, "Tail p-value of binomial test and its support","\n")
        }
       simpvaluesupports  <- apply(marginalsAndsizes,1,pvalueByBinoSupport)
       simpvalues <- unlist(sapply(simpvaluesupports,'[',1))
       MidPval  = unlist(sapply(simpvaluesupports,'[',2))
       randPval = unlist(sapply(simpvaluesupports,'[',3))
       
       ## extract supports
        pCDFlist = vector("list",m)
       for (ix in 1:m)  {
          pCDFlist[[ix]] = simpvaluesupports[[ix]][-(1:3)]
       }
      

      #################  one sided p-values  #####################
     if (! is.null(OneSide)) { 
       # one sided p-values
       simpvalues = double(m); MidPval = double(m); randPval = double(m)
       pCDFlist = vector("list",m)
       for (ia in 1:m){
            obsved = marginalsAndsizes[ia,]
            pvsupp = pvalOneSideBTSupport(obsved, Side = OneSide)
            simpvalues[ia] = pvsupp[1]   #pO, pMid, pRnd
            MidPval[ia] =  pvsupp[2]
            randPval[ia] =  pvsupp[3]
            pCDFlist[[ia]] = pvsupp[-(1:3)] # support
       }
      } # end of if for one-sided p-value              
     } # end of case 1, if (Test == "Binomial Test")  
  
  ######################################################################### ##########
  # Section 3: Get p-values and their supports from Fisher's Exact Test 
  #################################################################### ########## ##########
  # different test will be used as per method chosen
  if (Test == "Fisher's Exact Test") {
  
      #  cat("^^Using Fisher's Exact Test.","\n")
        # check for zero rows
        if (any(rowSums(data) == 0)) stop("^^Please remove zero rows in data matrix.","\n")
        
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
        
        ############## two-sided p-values and supports  ####################   
        # assign cell counts and marginals
        simallcellcounts <- cellcountsmarginals[[1]]  # each element of cellcountsmarginals[[1]] is a 2-by-2 matrix
 #       simallmarginals <- cellcountsmarginals[[2]]   # each element of cellcountsmarginals[[2]] is a 3 vector
        
        # compute two-sides pvalues and their supports
        if (is.null(OneSide))   {
         cat("^^^Computing two-sided p-value of Fisher's exact test and its support","\n") 
        } else {
         cat("^^^Computing one-sided", OneSide, "Tail p-value of Fisher's exact test and its support","\n")
        }   

      simpvalues = double(m); MidPval = double(m); randPval = double(m)
        pCDFlist = vector("list",m);         simpvaluesupports = vector("list",m)
       for (ia in 1:m)  {
          obsved = simallcellcounts[[ia]]
          pvsupp = pvalFETSupport(obsved)
           simpvalues[ia] = pvsupp[1]   #pO, pMid, pRnd
            MidPval[ia] =  pvsupp[2]
            randPval[ia] =  pvsupp[3]
          pCDFlist[[ia]] = pvsupp[-(1:3)]
       }
         
        
      ############## one-sided p-values and supports  ####################
      if (! is.null(OneSide)) {
       simpvalues = double(m); MidPval = double(m); randPval = double(m)
       pCDFlist = vector("list",m)
       for (ia in 1:m){
            obsved = cellcountsmarginals[[1]][[ia]]
            pvsupp = pvalOneSideFETSupport(obsved, Side = OneSide)
            simpvalues[ia] = pvsupp[1]   #pO, pMid, pRnd
            MidPval[ia] =  pvsupp[2]
            randPval[ia] =  pvsupp[3]
            pCDFlist[[ia]] = pvsupp[-(1:3)] # support
        }
       } # end of if one-sided p-value 
     }       # end of case 2, if (Test == "Fisher's Exact Test")


   #########################################################################
  # Section 5: Estimate pi0 and FDP from p-values and their supports 
  #################################################################### 
  
   cat("^^Estimating pi0 by generalized estimator ...","\n")
   pi0G = GenEstProp(simpvalues,pCDFlist,TuningRange)
  cat("^^ pi0 estimated by generalized estimator as: ","Gen",pi0G,"\n")
  
 ### FDP estimators
  cat("^^Implementing multiple testing procedures ...","\n")
  #BH 
  BHOrig <- BHFDRApp(simpvalues,FDRlevel)   # BH original
  #BH adaptive
  BHAdap <- BHFDRApp(simpvalues,FDRlevel/pi0G)   # BH adaptive with new pi0 est
  # heyse adjusted p-values
   pDBH <- HeyseAdjFunc(simpvalues, pCDFlist)
  # adaptive heyse   
  DBHAdap <- FDRAdjustedPval(pDBH,FDRlevel/pi0G) 
  # heyse   
  DBH <- FDRAdjustedPval(pDBH,FDRlevel) 
  ##################################################################
  # Section 7: put estimates in a dataframe and return it 
  ####################################################################
  cat("^^Gathering analysis results ... ","\n")

  # names for entries: thresh, FDP, number of discoveries, indices of discoveries
  BHH = list("pi0Est" = 1, "Threshold" = DBH[[2]], "NumberOfDiscoveries" = length(DBH[[1]][,2]),"IndicesOfDiscoveries" = DBH[[1]][,2])
  aBHH = list("pi0Est" = pi0G, "Threshold" = DBHAdap[[2]], "NumberOfDiscoveries" = length(DBHAdap[[1]][,2]),"IndicesOfDiscoveries" = DBHAdap[[1]][,2])
  
  BH = list("pi0Est" = 1, "Threshold" = BHOrig[[2]],"NumberOfDiscoveries" = length(BHOrig[[1]][,2]),"IndicesOfDiscoveries" = BHOrig[[1]][,2])
  aBH = list("pi0Est" = pi0G, "Threshold" = BHAdap[[2]],"NumberOfDiscoveries" = length(BHAdap[[1]][,2]),"IndicesOfDiscoveries" = BHAdap[[1]][,2])

   ################ estimators and procedures based on randomized and mid p-values
 # storey FDR proc to randomized p-value
    #  pi0SrandP = pi0est(randPval,pi0.method="smoother")$pi0    # this estimate sereverly underestimates after averaing in randp
      pi0SrandP = storeyPi0Est(0.5,randPval) # this less understimates pi0 than smoother and boostrap in qvaule pi0est
  StoryRndP <- StoreyFDREst(randPval,FDRlevel,pi0SrandP)
    SARP = list("pi0Est" = pi0SrandP,"Threshold" = StoryRndP[[2]][1],"NumberOfDiscoveries" = length(StoryRndP[[1]]),"IndicesOfDiscoveries" = StoryRndP[[1]])
   
  # adaptive BH based on mid p-values
   pi0SMidp = storeyPi0Est(0.5,MidPval)  #pi0est based mid pval
  BHAdapMidP <- BHFDRApp(MidPval,FDRlevel/pi0SMidp)
   aBHMidp =   list("pi0Est" = pi0SMidp, "Threshold" = BHAdapMidP[[2]],"NumberOfDiscoveries" = length(BHAdapMidP[[1]][,2]),"IndicesOfDiscoveries" = BHAdapMidP[[1]][,2])
  
  ## gathering results  ###########
AnalysisResults = list("aBHH"=aBHH,"BHH"=BHH,"BH"=BH,"aBH"=aBH,"pvalues" = simpvalues,"pvalSupp" = pCDFlist,"adjpval"=pDBH,"pi0Est"=pi0G, "RndPval" = randPval, "MidPval" = MidPval, "SARP" = SARP, "aBHmidP" = aBHMidp)
 
 # rerturn results 
 return(AnalysisResults)
    
      
} # end of main function