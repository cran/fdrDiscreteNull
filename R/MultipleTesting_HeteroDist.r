################################################################################
#######  Grouped and Weighted GenEst        ############
################################################################################
#load things in the main simulation file

GeneralizedEstimatorsGrouped = function(data_in = NULL,grpby = c("quantileOfRowTotal","kmeans","divergence"),ngrp_in = NULL,
                                         GroupMergeSize = 150, minGroupSize = 50,
                                          test_in = NULL, FET_via_in = NULL,OneSide_in = NULL, FDRlevel_in = NULL,
                                            eNetSize = NULL, unif_tol= 10^-3, Tunings =c(0.5,100)) 
{

   data_in = as.matrix(data_in)
   m = nrow(data_in)
  #################  stat of step 1: UNGrouped   ################
  cat("\n","^^^Implement UNGROUPED FDR procedures ...^^^","\n")
  UngroupedResults = GeneralizedFDREstimators(data=data_in, Test=test_in,
                                  FET_via = FET_via_in, OneSide = OneSide_in, FDRlevel=FDRlevel_in,TuningRange = Tunings) 
                                          
  pval_vec_ungrouped = UngroupedResults$pvalues   # pvalue vector
  pval_Supp_ungrped =  UngroupedResults$pvalSupp   # pvalue supports list 
  
  # pi0 estimated ungrouped by 2 methods; pi0Est_st_ugp is storey Est 
  pi0Est_ugp = UngroupedResults$pi0Est; 
  
  ### end of gathering results   UNGrouped

  #################  start Step 2: GROUPED ONLY  #################
  #################### start: preparations ##################
  cat("\n","^^ ----- Implement GROUPWISE Weighted procedure ------^^^","\n")

   ngrp = ngrp_in

   grpidx_list = vector("list",ngrp);    grpdata_list = vector("list",ngrp)
   grpres_list = vector("list",ngrp);    wgtPval_list = vector("list",ngrp)
   wgtSTPval_list = vector("list",ngrp)

    pi0Est_gped = double(0);    namepiEstgp = double(0);  pi0_ov_gp = 0; 

    grpId_dis = double(0);  grpId_STdis = double(0);  
    
    # caution: if it is FET, the 1st 2 columns are observed counts
    if (test_in == "Binomial Test") {
       totalCounts = rowSums(data_in)
       } else {
     if (test_in == "Fisher's Exact Test" & FET_via_in == "IndividualMarginals")  {
          totalCounts = rowSums(data_in[,1:2])
        }
      if (test_in == "Fisher's Exact Test" & FET_via_in == "PulledMarginals")  {
          totalCounts = rowSums(data_in)
       }
       }
      # caution ids_grps = sort(unique(grp.ids))
      diffCounts = abs(data_in[,1] - data_in[,2])  # if Poisson

     #################### end: preparations ##################

     ######################################################################################
     ############ Section 2: start: actual grouping and groupwise analysis   ################
     ######################################################################################

     # layout grouping strategy first

      if (grpby == "quantileOfRowTotal") {
          cat("^^Grouping by quantiles of row total counts","\n")
          brks = quantile(totalCounts, probs = seq(from=0,to=1,by=1/ngrp), na.rm = TRUE)
          }

     if (grpby == "kmeans") {
         cat("^^Grouping by kmeans on group differences and total counts","\n")
        kmgrps = kmeans(cbind(diffCounts,totalCounts),centers = ngrp)
        }

      ##  reference divergences
     if (grpby == "divergence") {
                 
             cat("^^Grouping by pairwise divergences. Be patient ...","\n")

             # insert:  when  |F_i - id| <= m^-3, take it as uniform   
              pv_supp_And_id_appr_unif = Div_Appr_Unif(pval_Supp_ungrped,test_in, unif_tol)  #(pvSpList=NULL,test_used=NULL,appr_unif_tol = NULL) 
              pv_supp_formatted = pv_supp_And_id_appr_unif[[1]]
              id_appr_unif = pv_supp_And_id_appr_unif[[2]]
              
              lgX = length(id_appr_unif)
              if (lgX > 0)  { 
                  cat("^^",lgX,"of p-value cdf's are very close to and will be identified as Uniform distribution","\n")
                  
                  if (lgX < m) {
                    id_pvDist_far_unif = (1:m)[-id_appr_unif]
                    ng_div =  ngrp-1
                    grpidx_list[[1]] = id_appr_unif   # first group via divergence
                    
                    lgXa =  length(id_pvDist_far_unif)
                    pv_supp_far_unif = vector("list",lgXa)
                     for (iX in 1:lgXa)  
                      pv_supp_far_unif[[iX]] = pv_supp_formatted[[iX]]
                    }  
                    if (lgX == m)
                     stop("^^ All p-value CDF's are within", unif_tol, "infinity-norm from Uniform distribution. Please use Storey's procedure...")
               } else { 
                     cat("^^ No p-value cdf is within", unif_tol, "infinity-norm distance from Uniform distribution","\n") 
                     ng_div =  ngrp 
                     pv_supp_far_unif = pval_Supp_ungrped
                     }
             # end insertion  
             
             ## group those that are far from Unif
             Divs_mat = GetDivergenceMatrix(pv_supp_far_unif)  # change scalfac when needed        
             div_max = max(Divs_mat)     # boxplot may cause memory surge  
             cat("^^ Maximal pairwise divergence is",div_max,"\n")
              if (div_max == 0)   cat("^^ Homogeneous null distributions since their distances are identically 0...","\n")
              
             # two cases for eNet size
             if (is.null(eNetSize))  {
               rad = div_max/(2*ng_div) 
               # make ngrp-1 groups out of those cdf's far from Unif  
               grpidx_list_div = eNetFull(Divs_mat, ng_div, GroupMergeSize, rad, minGroupSize) #(data, ngrp, mergesize, rad, mingrpsize)
               }  else {
                 grpidx_list_div = eNetFull(Divs_mat,ng_div, GroupMergeSize, eNetSize, minGroupSize) 
               }  # end of eNet
         
         } #   end if divergence

     #################### start of grouping and groupwise analysis ####################  
     for (j in 1:ngrp) {
       # start of groupwise analysis by bincounts
        if (grpby == "quantileOfRowTotal") {
           if (j < ngrp)
              grpidx_list[[j]] = which(totalCounts >= brks[j] & totalCounts < brks[j+1])
           if (j == ngrp)
                grpidx_list[[j]] = which(totalCounts >= brks[j] & totalCounts <= brks[j+1])
         }

       # start of groupwise analysis by kmeans
       if (grpby == "kmeans")  grpidx_list[[j]] = which(kmgrps$cluster == j)

       # assign grouping results via divergence
        if (grpby == "divergence") {
         # p-values with cdf's close to Unif are already as the first group grpidx_list[[1]]
         if (lgX > 0 & lgX <m) {
           if (j >1) {
            grpidx_list[[j]] = grpidx_list_div[[j-1]]  
           }
           # if no p-value cdf is close enough to Uniform 
           } else {
               grpidx_list[[j]] = grpidx_list_div[[j]]
             } 
           }
       ############ start of  grouped data, groupwise anlysis  #######
       grpdata_list[[j]] = cbind(data_in[grpidx_list[[j]],],grpidx_list[[j]])

       nc_kep = ncol(data_in)
       grpdata_in = grpdata_list[[j]][,1:nc_kep]
       cat("^^^ Checking number of hypotheses in each group ...", "\n")
       CheckIt(grpdata_in)   # check line

       ### when a group of data is non-empty, analyze them
       if (nrow(grpdata_in) == 0) {
          stop("^^^Group ",j," has no data, please adjust the number of groups to split data into...","\n")
        } else {
         # estimate groupwise pi0
         pi0eg =  GenEstProp(pval_vec_ungrouped[grpidx_list[[j]]],pval_Supp_ungrped[grpidx_list[[j]]],Tunings)
          
          # removed storey's pi0 est for each group
          pi0Est_gped = c(pi0Est_gped, pi0eg)
          namepiEstgp = c(namepiEstgp,paste('pi0Est_gp',j,sep=""))
       }   ##  only analyze non-empty grouped data     
   } ######## end of loop for (j in 1:ngrp) ######
   
   ## display estimated pi0
    cat("\n","^^--Estimated pi0 for each group by generalized estimator",pi0Est_gped,"--","\n") 
    if (any(pi0Est_gped==1))
      cat("^^ --- infinity weight appears --","\n") 
      
       ####################################################################################
       ######### step 3: start of weighting the p-values from the ungrouped procedure  ############
       ####################################################################################
       # grouping by divergence can create intersecting groups, for which HU 2010 JASA weighting does not work

       ### step 3.1: check if parition via divergence is reached
        non_partition = 0
        for (j1 in 2:ngrp) {
             for (j2 in 1:(j1-1)) {
                 non_partition = non_partition + any(grpidx_list[[j1]] %in% grpidx_list[[j2]])  }  }

       union_chk = length(unlist(grpidx_list))== m

        is_partition = 1- non_partition - !union_chk

        if (is_partition) {
               cat("\n","^^^----- Partition by grouping obtained. Ready for weighting ...","\n")
               }  else {
                cat("^^ ------Partition by grouping NOT obtained. Weighting will NOT be applied...","\n") }

    # if partition reached, then compute weights and following steps
    if (grpby != "divergence" | is_partition) {

         ### step 3.2: get overall pi0est and weight ungrouped pvalues
         for (j in 1:ngrp) {
            # pi0eg = grpres_list[[j]][[1]]$pi0Est
              pi0eg = pi0Est_gped[j];    wgt = pi0eg/(1-pi0eg);  

             ###### start of overall pi0est  #############
              pi0_ov_gp = pi0_ov_gp + pi0eg*length(grpidx_list[[j]])/m

             ### weigthing pvalues from ungrouped analysis
             pval_vec_ungrouped_gp = pval_vec_ungrouped[grpidx_list[[j]]]
             wgtPval_list[[j]] =  pval_vec_ungrouped_gp*wgt

             }  # end: overall pi0est and weight ungrouped pvalues


       ######### step 3.4: start: apply grouped and weighted scheme to BH if needed #####
       wgtPval_vec = unlist(wgtPval_list) 
         ### Now on grouped and weighted BH
         if ( pi0_ov_gp >= 1)   {
           cat("^^Overall pi0 est based on all groupwise ests is", pi0_ov_gp,". Setting weighted multiple testing results as zero", "\n")
            TDP_BHwg = 0;   FDP_BHwg = 0;   Dis_BH_gw = double(0)
            wFDR = list("pi0Est" = pi0_ov_gp, "Threshold" = 0, "NumberOfDiscoveries" = 0,"IndicesOfDiscoveries" = double(0))  
            }  else {
                 cat("^^Implement weighted procedure by MEANINGFUL weighting ...","\n")
                 FDRlevel_ada = FDRlevel_in/(1- pi0_ov_gp)
                 Dis_BH_gw <- BHFDRApp(wgtPval_vec,FDRlevel_ada)
              
               wFDR = list("pi0Est" = pi0_ov_gp, "Threshold" = Dis_BH_gw[[2]], "NumberOfDiscoveries" = length(Dis_BH_gw[[1]][,2]),"IndicesOfDiscoveries" = Dis_BH_gw[[1]][,2])   

            }  ## end weighted grouped BH by Gen est


    } # end if partition, then compute weights

  #################### step 4: start: collection all results #####################
  # if partition not reached, overall pi0 est is undefined
  if (grpby == "divergence" & non_partition != 0) {
      pi0_ov_gp = NA }

  nmpi0es = c("pi0E_GE","pi0E_gGE",namepiEstgp)
  pi0EstTmp = matrix(0,1,length(nmpi0es))
  pi0EstTmp[1,] = c(pi0Est_ugp,pi0_ov_gp,pi0Est_gped)
  colnames(pi0EstTmp) =  nmpi0es

 # get results
 results_grp = list("aBHH"=UngroupedResults$aBHH,"BHH"=UngroupedResults$BHH,"BH"=UngroupedResults$BH,
 "aBH"=UngroupedResults$aBH,
 "wFDR"=wFDR, "pi0estAll"=pi0EstTmp,"pvalues" =UngroupedResults$pvalues,"pvalSupp" =UngroupedResults$pvalSupp,
 "adjpval"=UngroupedResults$adjpval)   
} # end of main function
