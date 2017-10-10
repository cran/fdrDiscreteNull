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
  centers_id = double(0) ;    grp_idx = vector("list",0)

  while(itn & (ngrptmp <= ngp)) {
     cat("^^Forming group", ngrptmp, "....","\n")
   ball = vector("list",itn);  ball_size = double(itn)

   # compare itn balls
   for (i in 1:itn) {
        ball_id_tmp = which(abs(div_mat_in[i,]) <= rad)
        ball[[i]] = idx_left[ball_id_tmp]    # ids in ball from original 1:nrow(div_mat) indices
        ball_size[i] = length(ball[[i]])    }

    smax_id = which.max(ball_size);  id_mxball = idx_left[smax_id]
    ball_max = ball[[smax_id]];  size_mxball = ball_size[smax_id]

    cat("^^Ids: e-ball center", id_mxball,". Ball size:", size_mxball, "\n") # "ids in ball:",ball_max,"\n")

    centers_id = c(centers_id,id_mxball);     grp_idx[[ngrptmp]] = ball_max
    cnt_idx = c(cnt_idx,ball_max)
    # added check
    cat("^^# of Id used",length(cnt_idx),"\n")
    
    if (any(cnt_idx <0))     stop("^^negative id in current ids","\n") 
    
    idx_left = all_idx[-cnt_idx]  # take current ids from 1:m

    itn = length(idx_left)
    cat("^^# of ids left currently",itn,"\n") # ". Idx left:", idx_left,"\n")

    # Follwing is the best partition, case 1:  ngrptmp <= ngp-1
    if (ngrptmp <= ngp-1) {
       if (itn <= merge_size) {  
         # the above condition should be itn <= merge_size instead of itn < merge_size since it is when ngrptmp <= ngp-1
         cat("^^Can not form",ngp,"groups each with cardinality no less than",merge_size,"\n")
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
          cat("^# of id's left to group is", itn,",continue grouping","\n")
         div_mat_in = div_mat[idx_left,][,idx_left]
         }
    }  # end  if (ngrptmp <= ngp-1)

    # case 2:  ngrptmp = ngp
    if (ngrptmp == ngp) {
      # nothing left 
      if (itn == 0)
        cat("^^",ngp,"e-balls reached as requested","\n")

      if (itn > merge_size) {
       cat("^^---Some ids left for more groups. Regroup with larger e-net size --- ","\n")
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
  #    cat("^^Current # of groups",ngp_tmp,"\n")
      mgps_tmp = length(groups_tmp[[ngp_tmp]])

    } # end of while
    return(groups_tmp)
}   # end of function
 


########################## ########################## ##########################  
################## function:  divergence  matrix   ########################## 
########################## ########################## ##########################

GetDivergenceMatrix = function(scalfac=NULL,pvSpList=NULL)  {   
  
  lgA1 = length(pvSpList)     # caution: no longer m but lgA1
  lgt_div = (lgA1-1)*lgA1/2
  cat("^^Computing", lgt_div, "pairwise divergences. Takes time ...","\n")
  
  chi_mat = matrix(0,lgA1,lgA1); infNorm_mat = matrix(0,lgA1,lgA1); div_mat = matrix(0,lgA1,lgA1)
  
  # start of computing pairwise quantities
  for (i in 2:lgA1) { 
       sp_pv1 = pvSpList[[i]];  lg1 = length(sp_pv1)     #pvSpList has been formated
           
     for (j in 1:(i-1)) {       
       sp_pv2 = pvSpList[[j]]; lg2 = length(sp_pv2)    #pvSpList has been formated 
       
       chi_mat[i,j] = abs(lg1 - lg2)
     
       teva = union(sp_pv1,sp_pv2)  # evaluate cdf at union of pvalue supports 
       cdfn1 = pvalueDist(teva,sp_pv1);  cdfn2 = pvalueDist(teva,sp_pv2)
       
       infNorm_mat[i,j] = max(abs(cdfn1-cdfn2)) 
        }
    }   # end of computing pairwise quantities
    
      ## compute divergences
      chiWgt = max(chi_mat); infNormWgt = max(infNorm_mat)
      
      chi_mat = chi_mat + t(chi_mat) ;  infNorm_mat = infNorm_mat + t(infNorm_mat)
      
      if (chiWgt ==0 | infNormWgt == 0)   cat("^^ Homogeneous null distributions ...","\n")
  
      if (chiWgt ==0 & infNormWgt != 0)  {
         cat("^^ Homogeneous supports ...","\n") 
         div_mat = scalfac*(infNorm_mat/infNormWgt)
      }
      
      if (chiWgt != 0 & infNormWgt == 0)  {
         cat("^^ Homogeneous masses ...","\n") 
         div_mat = scalfac*(chi_mat/chiWgt)
      }
      
      if (chiWgt != 0 & infNormWgt != 0) {
        cat("^^ Heterogeneous supports and masses ...","\n")
        div_mat = scalfac*(chi_mat + infNorm_mat)/(chiWgt + infNormWgt)
      }
      
      div_mat = div_mat + t(div_mat) # since diagnal is zero, this is correct
      cat("^^Finished computing matrix of pairwise divergences...","\n")
       
      return(list(div_mat,chi_mat/chiWgt,infNorm_mat,chiWgt))
      
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
 
 
########################################################################
### Compute supremum norm to reference Unif
########################################################################
 Div_Ref_Unif = function(pvSpList=NULL,test_used=NULL) {

      lg3 = length(pvSpList)
      d_ref_unif = double(lg3)
      s_supp = double(lg3)

     # caution: psupport has different elements
       for (i in 1:lg3) {

         # using unif as reference distribution
         sp_pv_tmp = c(0,pvSpList[[i]])
         d_ref_unif[i] = max(abs(diff(sp_pv_tmp)))
         s_supp[i] = length(pvSpList[[i]])
       }
    return(list(d_ref_unif,s_supp))
  } # end of func
   
  
############################ #########################
### function compute divergence   as vector
######################################################      
GetDivVec = function(scalfac,pvSpList)  {
  m = length(pvSpList)   
  lg_div_vec = (m-1)*m/2
  
  cat("^^Computing", lg_div_vec, "pairwise divergences. Takes time ...","\n")
  
  chi_vec = double(lg_div_vec);   infNorm_vec = double(lg_div_vec)
  idxPair_div = matrix(0,nrow = lg_div_vec, ncol = 2);   div_vec = double(lg_div_vec) 
  
  # start of computing pairwise quantities
  for (i in 2:m) {
       sp_pv1 = pvSpList[[i]][-1][-1];  lg1 = length(sp_pv1)
        stloc = (i-1)*(i-2)/2  # start loc in the vector
     
     for (j in 1:(i-1)) {
       sp_pv2 = pvSpList[[j]][-1][-1]; lg2 = length(sp_pv2)
     
       teva = union(sp_pv1,sp_pv2)  # evaluate cdf at union of pvalues 
       cdfn1 = pvalueDist(teva,sp_pv1);  cdfn2 = pvalueDist(teva,sp_pv2)
       infNorm = max(abs(cdfn1-cdfn2));  infNorm_vec[stloc+j] = infNorm 
       
       chi = abs(lg1 - lg2);   chi_vec[stloc+j] = chi
       idxPair_div[stloc+j,] =  c(i,j)
        }
    }   # end of computing pairwise quantities
    
      ## compute divergences
      chiWgt = max(chi_vec); infNormWgt = max(infNorm_vec)
      
      if (chiWgt ==0 | infNormWgt == 0) {
        dev_vec = double(lg_div_vec)
        cat("^^ Homogeneous null distributions ...","\n")
      }
      
      if (chiWgt ==0 & infNormWgt != 0)  {
         cat("^^ Homogeneous supports ...","\n") 
         div_vec = scalfac*infNorm_vec/infNormWgt
      }
      
      if (chiWgt != 0 & infNormWgt == 0)  {
         cat("^^ Homogeneous masses ...","\n") 
         div_vec = scalfac*chi_vec/chiWgt
      }
      
      if (chiWgt != 0 & infNormWgt != 0) {
        cat("^^ Heterogeneous supports and masses ...","\n")
      div_vec = scalfac*(chi_vec/chiWgt + infNorm_vec/infNormWgt)
      }
      
      cat("^^Finished computing", lg_div_vec, "pairwise divergences...","\n")
      ## end of computing divergences
      
      # plot extreme statistics
     # chi_L = min(chi_vec); chi_LS = chi_L/chiWgt;  infNorm_L = min(infNorm_vec)
     # div_M = max(div_vec); div_L = min(div_vec)
      
    #  mins = c(chi_L,infNorm_L,div_L);  #  maxs = c(1,infNormWgt,div_M)
      
    #  cat("^^For divergences, (chiMin,chiMax):",chi_L,chiWgt,
     #           "(infNormMin,infNormMax):",infNorm_L,infNormWgt,"(divMin,divMax):",div_L,div_M,"\n")
      
    #  plot(mins,maxs,type="n", main = "Ranges of Differences")
    #  labs = c(expression(chi),expression(infinity),expression(delta))
    #  text(mins,maxs,labels=labs,cex=0.8,col= 1:3)
      par(mfrow = c(1,3))
      boxplot(infNorm_vec, main = "Diffs in Supremum Norms")
      boxplot(chi_vec/chiWgt,main = "Scaled Diffs in Supports",xlab="")    # bwplot requires library(lattice)
      boxplot(div_vec, main ="Divergences")
       
       
      return(list(div_vec,idxPair_div))
      
  }
       
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
      #### Subsection 4: Function to Get supports under true null
      #####                when test statistic is non-central hypergeometric
      ###################################################################################

      ## Function to get support of pvalues,  marginals is a vector as a list component
      # pvalue support is always computed under true null, so set pis=1
      pvalueSupport <- function(marginal)
      {   # if dim(marginals)!=c(3,3),report error
          masstmp <- dnoncenhypergeom(NA,marginal[1],marginal[2],marginal[3],1)
          mass <- masstmp[,2]
          lgh2 = length(mass)
          # change 11/21: matrix into double
           temp <- double(lgh2)
          # temp <- matrix(0,nrow=lgh2,ncol=1)
          for (i in 1:lgh2)
          {
            # temp[i] is the pvalue for i-th table given marginals
            temp[i] <- sum(mass[which(mass <= mass[i])])
          }
          # compute the expectation of pvalue
           meantmp <- sum(temp*mass)
           # sort pvalue support
          psupport <- sort(temp,decreasing=FALSE)

          # save mean to the fist entry
          support<- c(meantmp,psupport)
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

           # temp[i] is the two-sided pvalue for i-th table given marginals
           pvalue <- sum(mass[which(mass <= realizedmass)])
           if (is.na(pvalue)) print("pvalue is NaN")

          ### compute p value support
          lgh2 = length(mass)
          temp <- double(lgh2)
          for (i in 1:lgh2)
          {
            # temp[i] is the two-sided pvalue for i-th table given marginals
            temp[i] <- sum(mass[which(mass <= mass[i])])
          }
          # compute the expectation of pvalue
           meantmp <- sum(temp*mass)
           # sort pvalue support
          psupport <- sort(temp,decreasing=FALSE)

          # save mean to the fist entry
          support<- c(meantmp,pvalue, psupport)
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

     pvalOneSideFETSupport <- function(cellmatrix, lowertail = "F")
      {
           ns = cellmatrix[1,1]; nw = sum(cellmatrix[,1]); nb = sum(cellmatrix[,2]); nd = sum(cellmatrix[1,]);
          drx = dhyper(ns, nw, nb, nd)

          if (lowertail == "F") {
            # get support
            psupport = dhyper(0:ns, nw, nb, nd) + phyper(0:ns, nw, nb, nd,lower.tail = FALSE)
            psupport = unique(sort(psupport))
            # tail
            ptail =  phyper(ns, nw, nb, nd,lower.tail = FALSE)
            # type of p-value
               pvalO = ptail + drx

           } else {
             # get support
            psupport = phyper(0:nd, nw, nb, nd)
            psupport = unique(sort(psupport))
            # tail
             ptail =  phyper(ns, nw, nb, nd)
              # type of p-value
               pvalO = ptail

           }

          # save pval and support
          pvals = c(pvalO)
          support<- c(pvals,psupport)   # this is list now
          return(support)
          }

 ######## one-sided p-value and its support for Binomial ##########
     # p-value = pbinom(qbinom(.)) + u dbinom(qbinom(.))
     pvalOneSideBTSupport <- function(cellvector, lowertail = "F")
      {
           ns = cellvector[1]; nt = cellvector[2]
          drx = dbinom(ns, nt,0.5)

          if (lowertail == "F") {
            # get support
            psupport = pbinom(0:ns, nt,0.5,lower.tail = FALSE) + dbinom(0:ns, nt,0.5)
            psupport = unique(sort(psupport))
            # tail
            ptail =  pbinom(ns, nt,0.5,lower.tail = FALSE)
            # type of pval
             pvalO = ptail + drx
                
           } else {
             # get support
            psupport = pbinom(0:nt, nt,0.5)
            psupport = unique(sort(psupport))
            # tail
             ptail =  pbinom(ns, nt,0.5)
              # type of p-value
            pvalO = ptail
               
           }

          # save pvalue and support
          pvals = c(pvalO)
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
        
        return(pH)
      }