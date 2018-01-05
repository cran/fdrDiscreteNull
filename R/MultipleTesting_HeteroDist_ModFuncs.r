####################################################################
###########  function 6: check basic numeric properties  ############
####################################################################
# small function
CheckIt = function(vector_in) {
   tmp1 = any(is.na(vector_in) | is.nan(vector_in) | is.infinite(vector_in))
   cat("^^any na or nan:",tmp1,"...")
   
   if(is.list(vector_in)) stop("groupwise data is list","\t")
   
   if (is.vector(vector_in)) {
   cat("Min is:",min(vector_in)," Max is:", max(vector_in),"...")
   cat("length is",length(vector_in),"\t") }
   
   if (is.matrix(vector_in)) {
      nr = nrow(vector_in); nc = ncol(vector_in)
      cat("nrow:",nr,"; ncol:",nc,"\t")
     if (nc == 2) {
          rs1 =  rowSums(vector_in)
         cat("^^min rowSum:",min(rs1),"^^Max rowSum:",max(rs1),"\n") 
         }
      if (nc == 4) {
          rs1 =  rowSums(vector_in[,1:2]) # For fet based on individual marginals
         cat("^^min rowSum:",min(rs1),"^^Max rowSum:",max(rs1),"\n") 
         }    
     }   # third big if
   }  
 

###########################################################################
#### function: updated main grouping with specified minimal group size 
#######################################################################
### function
# FOR PRACTICAL PURPOSED to_merge is always true.

eNetBuilder = function(div_mat = NULL,ngp = NULL, merge_size = NULL,rad = NULL)
  {
  if (!is.matrix(div_mat))
    stop("^^ A matrix with entries as divergences is needed ...","\n")
  
  m_lt = nrow(div_mat)

  ngrptmp = 1; cnt_idx = double(0);  div_mat_in = div_mat
  all_idx = 1:m_lt;  # caution, no longer all_idx = 1:m 
  
  idx_left = all_idx ; itn = length(idx_left)
  centers_id = double(0) ;    grp_idx = vector("list",ngp)

  while(itn & (ngrptmp <= ngp)) {
    # cat("^^Forming group", ngrptmp, "....","\n")
   ball = vector("list",itn);  ball_size = double(itn)

   # compare itn balls
   for (i in 1:itn) {
        ball_id_tmp = which(abs(div_mat_in[i,]) <= rad)
        ball[[i]] = idx_left[ball_id_tmp]    # ids in ball from original 1:nrow(div_mat) indices
        ball_size[i] = length(ball[[i]])    }

    smax_id = which.max(ball_size);  id_mxball = idx_left[smax_id]
    ball_max = ball[[smax_id]];  size_mxball = ball_size[smax_id]

   # cat("^^Ids: e-ball center", id_mxball,". Ball size:", size_mxball, "\n") # "ids in ball:",ball_max,"\n")
    
    centers_id = c(centers_id,id_mxball);     grp_idx[[ngrptmp]] = ball_max
    cnt_idx = c(cnt_idx,ball_max)
    # added check
   # cat("^^# of Id used",length(cnt_idx),"\n")
    
    if (any(cnt_idx <0))     stop("^^negative id in current ids","\n") 
    
    idx_left = all_idx[-cnt_idx]  # take current ids from 1:m

    itn = length(idx_left)
   # cat("^^# of ids left currently",itn,"\n") # ". Idx left:", idx_left,"\n")

    # Follwing is the best partition, case 1:  ngrptmp <= ngp-1
    if (ngrptmp <= ngp-1) {
       if (itn <= merge_size) {  
         # the above condition should be itn <= merge_size instead of itn < merge_size since it is when ngrptmp <= ngp-1
       #  cat("^^Can not form",ngp,"groups each with cardinality no less than",merge_size,"\n")
         cat("^^---- Regroup with smaller e-net size ---- ","\n")
         rad = rad/2  # adjust this
         idx_left = 1:m_lt; itn = length(idx_left)
         ngrptmp = 0;  cnt_idx = double(0);  div_mat_in = div_mat
         }

       if (itn == merge_size & ngrptmp == ngp-1) {
         cat("^^Merging",itn,"ids.",ngp,"e-balls reached as requested","\n")
         grp_idx[[ngp]] =  idx_left
         itn = 0
        }
       
       if (itn > merge_size) {
         # cat("^# of id's left to group is", itn,",continue grouping","\n")
         div_mat_in = div_mat[idx_left,][,idx_left]
         }
    }  # end  if (ngrptmp <= ngp-1)

    # case 2:  ngrptmp = ngp
    if (ngrptmp == ngp) {
      # nothing left 
      if (itn == 0)
        cat("^^",ngp,"e-balls reached as requested","\n")

      if (itn > merge_size) {
     #  cat("^^---Some ids left for more groups. Regroup with larger e-net size --- ","\n")
        rad = 1.5*rad  # adjust this
        idx_left = 1:m_lt; itn = length(idx_left)
         ngrptmp = 0;  cnt_idx = double(0);  div_mat_in = div_mat
       }  
       
       if (itn <= merge_size) {
           grp_idx[[ngp]] = c(grp_idx[[ngp]],idx_left)
          cat("^^Merging",itn,"ids.",ngp,"e-balls reached as requested","\n")
           }
      } # end if (ngrptmp == ngp)

    ngrptmp = ngrptmp + 1

  } # end of while
      return (grp_idx)
} # end of fcuntion

### ###########################
eNetFull = function(metrics_in = NULL, ngrp_in = NULL, merge_size = NULL, rad_in = NULL, mgpsize_in = NULL) {

    if(is.null(mgpsize_in))
     stop("^^Specify minimal group size (MGS)")
    rad_tmp = rad_in;  mgps_tmp = 1  # mininal number of groups and group size

    while (mgps_tmp < mgpsize_in) {

      groups_tmp = eNetBuilder(metrics_in,ngrp_in,merge_size,rad_tmp)
      ngp_tmp = length(groups_tmp)
  #    check minimal group size
      mgps_tmp = min(sapply(groups_tmp,length))
      cat("^^Current achieved minimal group size is", mgps_tmp,"\n")
    } # end of while
    return(groups_tmp)
}   # end of function
 


########################## ########################## ##########################  
################## function:  divergence  matrix   ########################## 
########################## ########################## ##########################

GetDivergenceMatrix = function(pvSpList=NULL)  {   
  
  lgA1 = length(pvSpList)     # caution: no longer m but lgA1
  lgt_div = (lgA1-1)*lgA1/2
  cat("^^Computing", lgt_div, "pairwise divergences. Takes time ...","\n")
  
  infNorm_mat = matrix(0,lgA1,lgA1); div_mat = matrix(0,lgA1,lgA1)
  
  # start of computing pairwise quantities
  for (i in 2:lgA1) { 
       sp_pv1 = pvSpList[[i]];  lg1 = length(sp_pv1)     #pvSpList has been formated
           
     for (j in 1:(i-1)) {       
       sp_pv2 = pvSpList[[j]]; lg2 = length(sp_pv2)    #pvSpList has been formated 
     
       teva = union(sp_pv1,sp_pv2)  # evaluate cdf at union of pvalue supports 
       cdfn1 = pvalueDist(teva,sp_pv1);  cdfn2 = pvalueDist(teva,sp_pv2)
       
       infNorm_mat[i,j] = max(abs(cdfn1-cdfn2)) 
        }
    }   # end of computing pairwise quantities
      
     infNorm_mat_sym = infNorm_mat + t(infNorm_mat)
      
      cat("^^Finished computing matrix of pairwise divergences...","\n")       
      return(infNorm_mat_sym)
      
  }
  
######################### ########################## ##########################   
  ### Identify Approximate Uniform, this function also formats different psupport
  ######################### ########################## ########################## 
 
  Div_Appr_Unif = function(pvSpList=NULL,test_used=NULL,appr_unif_tol = NULL) {
      lg3 = length(pvSpList)
      pv_supp_formatted = vector("list",lg3)
      id_appr_unif = double(0)
       
          for (i in 1:lg3) {             
             sp_pv_tmp = c(0,pvSpList[[i]])
             if (max(abs(diff(sp_pv_tmp))) < appr_unif_tol)      
              id_appr_unif = c(id_appr_unif,i)
          }  
                      
    return(list(pvSpList,id_appr_unif))       
  } # end of func
 
 
       
 ## merged from mod func for generalized estimator
       ##############################################################################
      #### subsection 1: Function to get full tables for  association studies #########
      ###############################################################################
      fulltable <- function(twocoldata)
      {
        # m is the # rows in the data
        m <-nrow(twocoldata)
        # construct 3-3 full table
        tables<- array(0, dim=c(3,3,m))

        controlsum<-sum(twocoldata[,1])
        treatsum<-sum(twocoldata[,2])

         # get 2-by-2 table for each gene
          tables[1,1,]=twocoldata[,1]
          tables[1,2,]=twocoldata[,2]
          tables[2,1,]=controlsum - twocoldata[,1]
          tables[2,2,]= treatsum -twocoldata[,2]

         # get 3rd row and column for each gene
         tables[1,3,]= tables[1,1,]+ tables[1,2,]
         tables[2,3,]= tables[2,1,]+ tables[2,2,]

           tables[3,1,]=rep(controlsum,m)
           tables[3,2,]=rep(treatsum,m)
           tables[3,3,]=rep(treatsum+controlsum,m)

         return(tables)
      }


      ################################################################# ###########################
      #### Subsection 3: Function to Get marginals and cellcounts from generated counts
      #####             for association stuies
      ################################################################# ###########################

      getcellcountsandmarginals <- function(countveccontrol,countvectreat)
      {
        m <- length(countveccontrol)
        twocolsimdata <- cbind(countveccontrol,countvectreat)
        simtables <- fulltable(twocolsimdata)

        simallcellcounts <- vector("list",m)
        for (i in 1:m) {simallcellcounts[[i]] <- simtables[1:2,1:2,i]}

        simallmarginals <- vector("list",m)
        for (i in 1:m) {simallmarginals[[i]] <- c(simtables[1,3,i],simtables[2,3,i],simtables[3,1,i])}

         y2 <- list(simallcellcounts,simallmarginals)
          return(y2)
      }

 ###################################################################################
  #### Subsection 4: Function to two-sided p-value support for FET
  ###################################################################################

    pvalFETSupport <- function(cellmatrix)
      {

          ns = cellmatrix[1,1]; nw = sum(cellmatrix[,1]); nb = sum(cellmatrix[,2]); nd = sum(cellmatrix[1,])
#          x <- x[1, 1]; m <- sum(x[, 1]);    n <- sum(x[, 2]); k <- sum(x[1, ])
#        lo <- max(0, k - n);   hi <- min(k, m);       support <- lo:hi
#        logdc <- dhyper(support, m, n, k, log = TRUE)
          lo <- max(0, nd - nb);   hi <- min(nd, nw)
          support <- lo:hi
          # all probability masses for this distribution
          mass = dhyper(support, nw, nb, nd)
          # realized mass
          realizedmass = dhyper(ns, nw, nb, nd)

           # two-sided pvalue
           pvalue <- sum(mass[which(mass <= realizedmass)])
            # add: sum of probabilities of y such that P(y) < P(X)
           lessProb = sum(mass[which(mass < realizedmass)])
           # eqProb: Psum_y P(y) st P(y) = P(x)
           eqProb = pvalue - lessProb
           
           if (is.na(pvalue)) print("pvalue is NaN")

          ### compute p value support
          lgh2 = length(mass)
          temp <- double(lgh2)
          for (i in 1:lgh2)
          {
            # temp[i] is the two-sided pvalue for i-th table given marginals
            temp[i] <- sum(mass[which(mass <= mass[i])])
          }
           # sort pvalue support
          psupport <- unique(sort(temp,decreasing=FALSE))

          # u either runif or =0.5 or = 1
        randPval = lessProb + mean(runif(50))*eqProb  # average of 50 realizations
        randPval = min(1,randPval)
        
        # mid p-value
        MidPval = lessProb + 0.5*eqProb
        MidPval = min(1,MidPval)
          # return all three p-values
          pvals = c(pvalue,MidPval,randPval)
          
          # save probless to the fist entry
          support<- c(pvals, psupport)
          return(support)
          }

      ################################################################# ######### ######### ######### #########
      #### Subsection 5: Function to compuete deltas (lambda - F(lambda) for adjusted estimator         #########
      ################################################################# ######### #########  ######### #########
      # works for each support
      deviations <- function(lambda,asupport)
      {
            # remove the lst entry of asupport since it is the mean
            bigornot <- lambda >= asupport
            bignum <- sum(bigornot)

             ### if lambda smaller than the smallest pvalue, then delta[i]=lambda
            if ( bignum ==0)    { delta=lambda}

            # if lambda falls into an interval formed by successive points in support
            if ( bignum > 0)
            {      fall= max(which(bigornot))
             delta=lambda- asupport[-1][fall] }
          return(delta)
      }


      #################################################################################
      ####   Subsection 10: Distribution of the p-values  (Generally Applicable)  #########
      #################################################################################

       pvalueDist <- function(peva,support)
          {
            lgh=length(support)
            lgh3 <- length(peva)
            y <- double(lgh3)
            for (i in 1:lgh3)  {
              cpv <- peva[i] >= support
              # the key is to compare where t falls in the support
              if (sum(cpv)==0) {y[i]=0} # because t is strictly less than minimum of support
              else if (sum(cpv)==lgh) {y[i]=1} # because t is no less than maximum of support
              else # t falls in one interval
              {  y[i] <- support[max(which(cpv))]   }
            }
            return(y)
            }


      ################################################################################### #######################
      #### Subsection 14: pvalue supports, under true null  Test Stat ~ Bino for equality of two Poisson   #########
      ###################################################################################   #####################

      pvalueByBinoSupport <- function(marginal)
      {
          #  under null that two poisson rates are equal, the UMP test follows Binomial with p=0.5
          # get all possible masses
          mass <- dbinom(seq(from=0,to=marginal[2]),marginal[2],0.5)

          ## compute pvalue
          realizedmass <- dbinom(marginal[1],marginal[2],0.5)

           # two-sided pvalue
           pvalue <- sum(mass[which(mass <= realizedmass)])
           # add: sum of probabilities of y such that P(y) < P(X)
           lessProb = sum(mass[which(mass < realizedmass)])
           # eqProb: Psum_y P(y) st P(y) = P(x)
           eqProb = pvalue - lessProb
           
           # check p-value
           if (is.na(pvalue)) print("pvalue is NaN")

          ### compute p value support
          lgh2 = length(mass)
          temp <- double(lgh2)
          for (i in 1:lgh2)
          {
            # temp[i] is the two-sided pvalue for i-th table given marginals
            temp[i] <- sum(mass[which(mass <= mass[i])])
          }
          # sort pvalue support
          psupport <- sort(temp,decreasing=FALSE)
          
          # u either runif or =0.5 or = 1
        randPval = lessProb + mean(runif(50))*eqProb  
        randPval = min(1,randPval)
        
        # mid p-value
        MidPval = lessProb + 0.5*eqProb
        MidPval = min(1,MidPval)
          # return all three p-values
          pvals = c(pvalue,MidPval,randPval)
           
          # save mean to the fist entry
          support<- c(pvals, psupport)
          return(support)
          }


     # storey: the new qvalue package stops with error when there is p-value being exactly 1
     # so the following is used
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
########## added new scheme for generalized estimator
GenEstProp <- function(pvector,psupports,tunings=c(0.5,100))
  {
    # m is the coherent length of the argument
    m = length(pvector)

    # define smallest guiding value
    minSupp = unlist(lapply(psupports, min))
    maxmin = max(minSupp)
    cat("^^Maximum of the minimum of each support is: ",maxmin,"\n")
    
    if (maxmin == 1) {
      cat("^^At lest one p-value CDF is a Dirac mass ...","\n")      
      singPsupp = which(minSupp == 1)      
      maxminA = max(minSupp[-singPsupp])
      
      cat("^^Second maximum of the minumum of each non-singleton support is: ",maxminA,"\n")
      if (maxminA < 0.5) {
        tunings = seq(maxminA+tunings[1]*(0.5-maxminA), 0.5,length.out=tunings[2])
      } else{ 
        tunings = c(maxminA)  
      }
      
    } else {
      singPsupp = double(0)
      if (maxmin < 0.5) {
        tunings = seq(maxmin+tunings[1]*(0.5-maxmin), 0.5,length.out=tunings[2])
      } else {  
        tunings = c(maxmin)  
      }
      
    } 
    
    # define which supports are to be searched
    if (length(singPsupp) == 0) {
      schLoc = 1:m
    }  else {
      schLoc = (1:m)[-singPsupp]
    } 
    
    # start search
    Lt = length(tunings)
    est = double(Lt)
    cuts = double(m)
    
    for (j in 1:Lt) {
      lamb = tunings[j]
      
      for (i in schLoc) {
        psupp = psupports[[i]]
        if (sum(psupp<=lamb) >0 ) {
          tmp = unique(max(psupp[which(psupp<=lamb)]))
          loc = max(which(psupp==tmp))
          cuts[i]= psupp[loc]
        } else  {    cuts[i] = NA }
        
      }
      cutstmp = cuts[!is.na(cuts)]
      
      esttmp = sum(as.numeric(pvector[!is.na(cuts)]>cutstmp)/(1-cutstmp))/m + 
        1/((1-lamb)*m)+ sum(is.na(cuts))/m+length(singPsupp)/m
      
      est[j] = min(1,esttmp)
    }
    #genest = min(est[est>0])
    genest = mean(est[!is.nan(est)])
    return(genest)
  }


    ###################################################
    ####### Subsection 19: BH FDR procedure        ########
    ######################################################
    # BH procedure requires input of pvalues and their indices
      # BH procedure requires input of pvalues and their indices
      BHFDRApp <- function(Pvals,FDRlevel)
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
        print ("No rejections by BH procedure")        # when there are no rejections
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


     #### rejection based on adjusted p-values
      FDRAdjustedPval <- function(Pvals,FDRlevel)
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
        cmp3 <- PvalAndIdxOrd[,1] <= FDRlevel
        scmp3 <- sum(cmp3)
    
       # collect rejections if any
        if (scmp3 == 0) {   
        print ("No rejections by BH")        # when there are no rejections
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
               rejAndTresh <- list(BHrej,PvalAndIdxOrd[,1][r])    
              }
            return(rejAndTresh)
      }
  ############################################################# ###########################
      #### Subsection 23: Function to Get marginals and cellcounts from generated counts
      #####             for differential expression
      ################################################################# ###########################

 getcellcountsandmarginals_DE <- function(data_in)
      {
        m <-nrow(data_in)
        # construct 3-3 full table
        tables<- array(0, dim=c(3,3,m))

        tables[1,1,] = data_in[,1]
        tables[2,1,] = data_in[,2]
        tables[1,2,] = data_in[,3] - data_in[,1]
        tables[2,2,] = data_in[,4] - data_in[,2]
        tables[1,3,]= data_in[,3]
        tables[2,3,]= data_in[,4]
        tables[3,3,] = data_in[,3] + data_in[,4]
        tables[3,1,] = data_in[,1] + data_in[,2]
        tables[3,2,] = tables[3,3,] - tables[3,1,]
        
        simallcellcounts <- vector("list",m)
        for (i in 1:m) {simallcellcounts[[i]] <- tables[1:2,1:2,i]}

        simallmarginals <- vector("list",m)
        for (i in 1:m) {simallmarginals[[i]] <- c(tables[1,3,i],tables[2,3,i],tables[3,1,i])}

         y2 <- list(simallcellcounts,simallmarginals)
          return(y2)
      }
      
      
######## one-sided p-value and its support for FET ##########
     # p-value = pbinom(qbinom(.)) + u dbinom(qbinom(.))

     pvalOneSideFETSupport <- function(cellmatrix, Side = "Right")
      {
           ns = cellmatrix[1,1]; nw = sum(cellmatrix[,1]); nb = sum(cellmatrix[,2]); nd = sum(cellmatrix[1,]);
          drx = dhyper(ns, nw, nb, nd)
          lo <- max(0, nd - nb);   hi <- min(nd, nw)
          #          x <- x[1, 1]; m <- sum(x[, 1]);    n <- sum(x[, 2]); k <- sum(x[1, ])
#        lo <- max(0, k - n);   hi <- min(k, m);       support <- lo:hi
#        logdc <- dhyper(support, m, n, k, log = TRUE)

          if (Side == "Right") {
            # get support
            psupport = dhyper(lo:ns, nw, nb, nd) + phyper(lo:ns, nw, nb, nd,lower.tail = FALSE)
            psupport = unique(sort(psupport))
            # tail
            ptail =  phyper(ns, nw, nb, nd,lower.tail = FALSE)
            # type of p-value
               pvalO = ptail + drx
                 pvalMid = ptail + 0.5*drx
                pvalRnd = ptail + mean(runif(50))*drx
           } 
           
           if (Side == "Left") {
             # get support
            psupport = phyper(lo:nd, nw, nb, nd)
            psupport = unique(sort(psupport))
            # tail
             ptail =  phyper(ns, nw, nb, nd)
              # type of p-value
               pvalO = ptail
                 pvalMid = ptail - drx + 0.5*drx
               pvalRnd = ptail - drx + mean(runif(50))*drx
           }

          # save pval and support
          pvals = c(pvalO,pvalMid,pvalRnd)
          support<- c(pvals,psupport)   # this is list now
          return(support)
          }

 ######## one-sided p-value and its support for Binomial ##########
     # p-value = pbinom(qbinom(.)) + u dbinom(qbinom(.))
     pvalOneSideBTSupport <- function(cellvector, Side = "Right")
      {
           ns = cellvector[1]; nt = cellvector[2]
          drx = dbinom(ns, nt,0.5)

          if (Side == "Right") {
            # get support
            psupport = pbinom(0:ns, nt,0.5,lower.tail = FALSE) + dbinom(0:ns, nt,0.5)
            psupport = unique(sort(psupport))
            # tail
            ptail =  pbinom(ns, nt,0.5,lower.tail = FALSE)
            # type of pval
             pvalO = ptail + drx
               pvalMid = ptail + 0.5*drx
                pvalRnd = ptail + mean(runif(50))*drx 
           } 
           
           if (Side == "Left") {
             # get support
            psupport = pbinom(0:nt, nt,0.5)
            psupport = unique(sort(psupport))
            # tail
             ptail =  pbinom(ns, nt,0.5)
              # type of p-value
            pvalO = ptail
               pvalMid = ptail - drx + 0.5*drx
               pvalRnd = ptail - drx + mean(runif(50))*drx   
           }

          # save pvalue and support
          pvals = c(pvalO,pvalMid,pvalRnd)
          support<- c(pvals,psupport)   # this is list now
          return(support)
          }
          

 ### heyse2011 adjusted p-values
    HeyseAdjFunc = function (pvec, pSupps)
      {
        m = length(pvec)
        if (m ==1)
          stop("^^^At least two p-values and their supports are needed.")
        # supports[[i]] is the support for pvec[i]; entries of supports[[i]] are ascendingly ordered wrt index
        
        # order p-values ascendingly
        pvecAscend =  pvec[order(pvec,decreasing = F)]
        
        # create BH sequence of p-values in terms of Heyse2011 rephrase; m*p_(k)/k, where p_(m) is the largest
        pBHAscend = m*pvecAscend/(1:m)
        
        # obtain Q(t) = \sum_{j=1}^m F_j(t)
        # obtain matrix F_j(p_(i))
        stepf = lapply(pSupps, function(x) stepfun(x, c(0, x)))
        p.mat = matrix(NA, m, m)
        for (i in 1:m)  {
          p.mat[i, ] = stepf[[i]](pvecAscend)
        }
        # sum of each column; obtain Q(p_(i)) = \sum_{j=1}^m F_j(p_(i))
        p.sums = apply(p.mat, 2, sum)
        # Q(p_(i))/i
        p.sums = p.sums/(1:m)
        
        # obtain Heller 2012 with orginal p-vec ordering
        pR = double(m)
        pR[m:1] = cummin(p.sums[m:1])
        op = order(pvec)
        pR = pmin(1,pR[order(op)])
        
        
        ##### obtain Heyse2011 with orginal p-vec ordering ####
        pH = double(m)
        pH[m] = pBHAscend[m]
        pH[(m-1):1] = pmin(cummin(pBHAscend[m:2]),p.sums[(m-1):1])
        pH = pmin(1, pH[order(op)])
        
        # choose which to return
        return(pH)
      }
 
 #######################################################################
    #### Subsection 18:  Storey's procedure               #########
    #######################################################################

      StoreyFDREst <- function (simpvalues,FDRlevel,pi0est)
      {
         m <- length(simpvalues)
        # step length in search depends on minimal difference between p values
         sortedpvalvec <- sort(simpvalues)
         diffsortedpvals <- diff(sortedpvalvec)
         # since diff does not compute the differece between minimal p value and zero
         # it should be added, so step length in thresh is the minimum of
         # minimal diff and min p value
         threshstep <- min(min(diffsortedpvals[diffsortedpvals > 0]),min(sortedpvalvec))

        # it is known that FDR threshold is no less than bonferroni, but no larger than
        # FDRlevel. the step length should be sensitive to distances bewteen p values to
        # make practical change in R(t) as t changes
        # sometimes threshstep can be so small that R itself can not handle it

        # Step 1: crude search first   do not start from min(sortedpvalvec) when m is small
        threshveccrude <- seq(from= 0,to=1,by=10^-3)
        ecdfcrude0 <- max(sum(sortedpvalvec <= threshveccrude[1]),1)/m

        istcrude <-1
        stFDPcrude <- pi0est*threshveccrude[1]/ecdfcrude0
         stFDPcrude0 <- stFDPcrude
        while (stFDPcrude <= FDRlevel)
        {
          istcrude <- istcrude+1;
          # check if mat threshol (t=1) is reach, is so, break
          if (istcrude > length(threshveccrude)){print("max threshold (t=1) in crude search for Storey's estimators exhausted");break}
          steCDFcrude <- max(sum(sortedpvalvec <= threshveccrude[istcrude]),1)/m
          stFDPcrude <- pi0est*threshveccrude[istcrude]/steCDFcrude
          if (is.na(stFDPcrude) | is.nan(stFDPcrude)) print("NA or NaN occured as FDP of Storey's estimators")

        }
        if (istcrude == 2) { stsupthreshcrude <- threshveccrude[1]
             stFDPcrudeAtSThresh <- stFDPcrude0    }     else {
              stsupthreshcrude <- threshveccrude[istcrude-1]
             ecdfstsupthreshcrude <- max(sum(sortedpvalvec <= stsupthreshcrude),1)/m
        stFDPcrudeAtSThresh <- pi0est*stsupthreshcrude/ecdfstsupthreshcrude  }
      #  cat("stFDPcrudeAtSThresh is",stFDPcrudeAtSThresh,"at sup thresh",stsupthreshcrude,"\n")

        ############################################
        #### step 2: after a crude threshold has been located, do a finer search from the crude sup threshold
        ### chen's fdr estimator
        # check if mat threshol (t=1) is reach, is so, break

         # use while loop for storey's
         if (stsupthreshcrude == 1) {  print("max threshold (t=1) for Storey's estimators exhausted, no fine search in second round")
            stsupthresh <- stsupthreshcrude
         stFDPAtSThresh <- stFDPcrudeAtSThresh} else {
        threshvecst <- seq(from = stsupthreshcrude,to=threshveccrude[istcrude],
                               by= min((threshveccrude[istcrude]-stsupthreshcrude)/m,FDRlevel/m))
         ecdf0st <- max(sum(sortedpvalvec <= threshvecst[1]),1)/m
         stFDP <- pi0est*threshvecst[1]/ecdf0st
          stFDP0 <- stFDP
         ist <-1
        while (stFDP <= FDRlevel)
        {
          ist <- ist+1;
          if (ist > length(threshvecst)) break
          steCDF <- max(sum(sortedpvalvec <= threshvecst[ist]),1)/m
          stFDP <- pi0est*threshvecst[ist]/steCDF
          if (is.na(stFDP) | is.nan(stFDP)) print("NA or NaN occured as FDP for Storey's estimators in fine search")
        }
        if (ist==2) {stsupthresh <- threshvecst[1]
          stFDPAtSThresh <- stFDP0     } else {
        stsupthresh <- threshvecst[ist-1]
         ecfstsupthresh <- max(sum(sortedpvalvec <= stsupthresh),1)/m
        stFDPAtSThresh <- pi0est*stsupthresh/ecfstsupthresh }
       #  cat("stFDPAtSThresh",stFDPAtSThresh,"at sup thresh",stsupthresh,"\n")
      } # end of esle


        # get discoveries
        stDiscovery <- which(simpvalues <= stsupthresh)
        if (length(stDiscovery)==0) print("No discoveries found by Storey's procedures")

        StDiscoveries <- list(stDiscovery,c(stsupthresh,stFDPAtSThresh))
        return(StDiscoveries)

      }
      
# storey: the new qvalue package stops with error when there is p-value being exactly 1
     # so the following is used
     storeyPi0Est <- function(lambda,pvector)
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
