###################################################
####### Subsection 19: BH FDR procedure        ########
######################################################

BH <- function(Pvals,FDRlevel)
  {   
    if (any(is.na(Pvals)) | any(is.nan(Pvals))) { cat("^^NA or NAN in pvalues... ","\n")
    print(Pvals[is.na(Pvals)])
    }
    
    if (any(is.na(FDRlevel)) | any(is.nan(FDRlevel))) cat("^^NA or NAN FDRlevel... ","\n")
       
   if (is.vector(Pvals))  lgh3 <- length(Pvals)
    if (is.matrix(Pvals) | is.data.frame(Pvals))  lgh3 <- dim(Pvals)[1] 
    PvalAndIdx <- cbind(Pvals,seq(1:lgh3))
    PvalAndIdxOrd <- PvalAndIdx[order(PvalAndIdx[,1]),]
    
 #   cat("^^Dims of PvalAndIdxOrd is:",as.vector(dim(PvalAndIdxOrd)),"\n")

    BHstepups <- seq(1:lgh3)*(FDRlevel/lgh3)
    cmp3 <- PvalAndIdxOrd[,1] <= BHstepups
    scmp3 <- sum(cmp3)

   # collect rejections if any
    if (scmp3 == 0) {   
    print ("No rejections made by BH procedure")        # when there are no rejections
    rejAndTresh <- list(matrix(numeric(0), ncol=2,nrow=0),0)       
    }  else   {  r <- max(which(cmp3))        
        #  cat("^^^ Minimax index in BH is:",r,"\n")
         # cat("^^^ Minimax threshold in BH is:",BHstepups[r],"\n")
               if (r==1) {
                # when r =1, a row is chosen, so should change it into a matrix of two columns
                 BHrej = as.matrix(t(PvalAndIdxOrd[1:r,]))  
               #  print(BHrej)
               } else {
                  BHrej <- PvalAndIdxOrd[1:r,]
               }
           rejAndTresh <- list(BHrej,BHstepups[r])    
          }
        return(rejAndTresh)
  }

####################################################
#### BH plus procedure                 ############
#############################
BHplus = function(pvals, fmax, pi0hat, alpha) {
    m = length(pvals)
    if (min(fmax)==1) stop("^^ All p-value CDF's are Dirac mass")
    stepConst = double(m)
  
    stepUB = alpha*(1:m)/(pi0hat*m)  # stepup bounds at balpha[ia] 
    # compute critical values
    for (i in 1:m) {
        chk = (fmax <= stepUB[i])
        if (sum(chk)>0) {
          stepConst[i] = fmax[max(which(chk))+1]-10^(-7)
        } else {
         stepConst[i] = 0
        }
      }  # out of compute critical values 

    ## compare
    PvalAndIdx <- cbind(pvals,1:m)
    PvalAndIdxOrd <- PvalAndIdx[order(PvalAndIdx[,1]),]
    cmp <- PvalAndIdxOrd[,1] <= stepConst
    scmp <- sum(cmp)
    
    # collect rejections if any
    if (scmp == 0) {   
        print ("No rejections made by BH+ procedure")        # when there are no rejections
        rejAndTresh <- list(matrix(numeric(0), ncol=2,nrow=0),0)       
      }  else   {  r <- max(which(cmp))        
           if (r==1) {
            # when r =1, a row is chosen, so should change it into a matrix of two columns
             BHPlusrej = as.matrix(t(PvalAndIdxOrd[1:r,]))  
           #  print(BHrej)
           } else {
              BHPlusrej <- PvalAndIdxOrd[1:r,]
           }
           rejAndTresh <- list(BHPlusrej,stepConst[r])   
      } 
    return(rejAndTresh)                          
 }

#################################################################
#### Function to get two-sided p-val and its support for FET
pvalFET <- function(cellmatrix)
  {

      ns = cellmatrix[1,1]; nw = sum(cellmatrix[,1]); nb = sum(cellmatrix[,2]); nd = sum(cellmatrix[1,])
#          x <- x[1, 1]; m <- sum(x[, 1]);    n <- sum(x[, 2]); k <- sum(x[1, ])
#        lo <- max(0, k - n);   hi <- min(k, m);       support <- lo:hi
#        logdc <- dhyper(support, m, n, k, log = TRUE)
      lo <- max(0, nd - nb);   hi <- min(nd, nw)
      support <- lo:hi
      # all probability masses for this distribution
      mass = dhyper(support, nw, nb, nd)

      ### compute p value support
      lgh2 = length(mass)
      temp <- double(lgh2)
      tempA <- double(lgh2)
      for (i in 1:lgh2)
      {
        # temp[i] is the two-sided mid pvalue for i-th table given marginals
        etmp = sum(mass[which(mass == mass[i])])
        temp[i] <- sum(mass[which(mass < mass[i])]) + 0.5*etmp
        tempA[i] = temp[i]+0.5*etmp
      }
       # match mid pval and cdf evaluated at mid pval
      ui = unique(which(unique(temp) %in% temp))
      temp = temp[ui]
      tempA = tempA[ui]
      ## order and matching them
      orMidp = order(temp)
      psupport <- temp[orMidp] # mid p-values
      psuppA = tempA[orMidp] # F(midp) = convp
                
      # realized mass
      realizedmass = dhyper(ns, nw, nb, nd)
       # two-sided pvalue
       pvalue <- sum(mass[which(mass <= realizedmass)])
        # add: sum_y P(y) over y such that P(y) < P(X)
       lessProb = sum(mass[which(mass < realizedmass)])
       # eqProb: sum_y P(y) over y such that P(y) = P(x)
       eqProb = pvalue - lessProb
       
      # randomized p-value
      randPval = pvalue - mean(runif(1))*eqProb  # definition from Dickhaus et al Lemma 2
      randPval = min(1,randPval)

      # mid p-value
      MidPval = lessProb + 0.5*eqProb # definition from HWang et al 2001
      MidPval = min(1,MidPval)
      # return all three p-values
      pvals = c(pvalue,MidPval,randPval)

      # save probless to the fist entry
      support<- list(pvals, psupport, psuppA)
      return(support)
      }

##################################################################################
####  Function to get two-sided pval and its support for binomial test    #########
pvalBT <- function(marginal)
  {
      # get all possible masses; ## correction (suggested by Florian Junge) to R's intrinsic numerical value inprecision
      # mass <- dbinom(seq(from=0,to=marginal[2]),marginal[2],0.5)
        spTmp = seq(from=0,to=marginal[2])
        d <- dbinom(spTmp,marginal[2],0.5,TRUE) ### using log of probs provides more accurate values
        ord <- order(d)
        d.ordered <- d[ord]
        
        # when marginal[2] is odd, the null CDF is symmetric and each prob appears twice
        relDev <- abs(diff(d.ordered)) ### check differences
        idx <- which(relDev > 0 & relDev <= .Machine$double.eps*2) ### "2" has proven to be a good choice
        # slight numerical correction
        ler = length(idx)  
        if (ler >0) {  
          for (jx in 1:ler) {
            d.ordered[c(idx[jx], idx[jx] + 1)] <- sum(d.ordered[c(idx[jx], idx[jx] + 1)])/2        
          } } 
        
        d.ordered <- exp(d.ordered)
        s <- sum(d.ordered) # check total prob
        if(s != 1) d.ordered <- d.ordered/s     
        d <- d.ordered[order(ord)]              ### restore original order (if necessary)
        #### end of correction
         mass = d      # re-assign
                  
      ### compute p value support
      lgh2 = length(mass)
      temp <- double(lgh2)
      tempA <- double(lgh2)
      for (i in 1:lgh2)
      {
        # temp[i] is the two-sided mid pvalue for i-th table given marginals
        etmp = sum(mass[which(mass == mass[i])])
        temp[i] <- sum(mass[which(mass < mass[i])])+ 0.5*etmp
        tempA[i] = temp[i]+0.5*etmp
      }
      # match mid pval and cdf evaluated at mid pval
      ui = unique(which(unique(temp) %in% temp))
      temp = temp[ui]
      tempA = tempA[ui]
      ## order and matching them
      orMidp = order(temp)
      psupport <- temp[orMidp] # mid p-values
      psuppA = tempA[orMidp] # F(midp) = convp
      

      ## compute pvalue
      # realizedmass <- dbinom(marginal[1],marginal[2],0.5)
      realizedmass <- mass[which(spTmp==marginal[1])]
      
       # two-sided pvalue
       pvalue <- sum(mass[which(mass <= realizedmass)])
       # add: sum of probabilities of y such that P(y) < P(X)
       lessProb = sum(mass[which(mass < realizedmass)])
       # eqProb: sum_y P(y) st P(y) = P(x)
       eqProb = pvalue - lessProb
       # check p-value
       if (is.na(pvalue)) print("pvalue is NaN")
       
      # ranomized p-value
      randPval = pvalue - mean(runif(1))*eqProb
      randPval = min(1,randPval)

      # mid p-value
      MidPval = lessProb + 0.5*eqProb
      MidPval = min(1,MidPval)
      # return all three p-values
      pvals = c(pvalue,MidPval,randPval)

      # save mean to the fist entry
      support<- list(pvals, psupport,psuppA)
      return(support)
  }


### maximum and minum of of p-value cdfs
ExtremeCDF = function(pMidSupps, pSupps,epsi){
  # step 1: take the union of all mid pval supports and 1000 equally spaceed as grid 
      m = length(pMidSupps) # number of supports
   midpAll = sort(unique(unlist(pMidSupps)))
   manualPoints = seq(from=0,to=1,by=0.001) # add 1000 equally spaced points
   midpAll = sort(unique(c(midpAll,manualPoints))) # combine
   midpAll = midpAll[midpAll<1] # remove 1 
   lgh3 = length(midpAll)
  # step 2: evaluate each mid pval cdf at this union
   evalCDFMat = matrix(0,m,lgh3) # store evaluation in a row
   DeviationMat = matrix(0,m,lgh3) # each row is F_j(t) -t
    for (j in 1:m) {
      # pick range of a mid p and its cdf
      currentConvPSupp = pSupps[[j]] # cdf of mid p maps mid p to conv p
      currentMidpSupp = pMidSupps[[j]]
      lgh = length(currentMidpSupp) # check where evaluation point falls into successive mid p-vals
      # evaluate current cdf for each element in midpAll
      y <- double(lgh3)  # to store evaluation
      for (i in 1:lgh3)  {
        currnetEvalPoint = midpAll[i]
        cpv <- currnetEvalPoint >= currentMidpSupp  # compare current eval point wrt range of mid p-val
        # the key is to compare where t falls in the support
        if (sum(cpv)==0) {y[i]=0} # because t is strictly less than minimal val of mid p
        else if (sum(cpv)==lgh) {y[i]=1} # because t is no less than maximal val of mid p
        else # t falls in one interval inside (min,max) of range of mid p
        {  y[i] <- currentConvPSupp[max(which(cpv))]   }  #since P(midp <= a given midp) = the given conv p
      } # end of evaluating a cdf of mid pval
      # store evaluation for j-th mid pval cdf 
      evalCDFMat[j,] = y
      DeviationMat[j,] = y - midpAll
    } # end of evaluating all cdfs 
  
 # step 3: take vector max of all rows of evalCDF matrix to get extremal cdf's
 cdfMax = apply(evalCDFMat, 2, max)
 cdfMin = apply(evalCDFMat, 2, min)
 
 # step 4: determine tau
  LookForTau = pmax(cdfMax, 1-cdfMin) <= midpAll
  if (sum(LookForTau)>0) {
        tau = min(midpAll[which(LookForTau)])
        tauVon = 0    # theoretical
        cat("^^Theoretical tau ",tau,"\n")
      } else {
         # find tau such that most of F_i at tau are close enough to tau
        DevCount = apply(DeviationMat, 2, function(x) {sum(abs(x)) <= epsi})
        if (max(DevCount)>0) {
          tau = max(midpAll[which(DevCount==max(DevCount))])
          tauVon = 1  # approximate
          cat("^^Approximate tau ",tau,"\n")
          } else {
            tau = 0.5
            tauVon = 2  # prespecified
            cat("^^Prespecified tau ",tau,"\n")
          } 
      }
 return(list(cdfMax,cdfMin,tau,tauVon))
}   

############################################################
  storeyEst <- function(lambda,pvector)
      {
          # m is the coherent length of the argument
           m <- length(pvector)

          over <- m*(1-lambda)
          stunmin <- sum(pvector > lambda)/over 
          st <- min(1,stunmin)
          if (is.nan(st))  print("Storey's estimator of pi0 is NaN")
          if (is.na(st))  print("Storey's estimator of pi0 is NA")

          return(st)
      }  