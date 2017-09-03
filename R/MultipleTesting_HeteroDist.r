################################################################################
#######  Grouped and Weighted GenEst        ############
################################################################################
#load things in the main simulation file

GeneralizedEstimatorsGrouped = function(data_in = NULL,grpby = c("quantileOfRowTotal","kmeans","divergence"),ngrp_in = NULL,
                                          test_in = NULL,FET_via_in = NULL,lowerTail_in = NULL, FDRlevel_in = NULL,
                                            RefDivergence = NULL,eNetSize = NULL, unif_tol= 10^-3, Tunings =c(0.5,100)) 
{

   data_in = as.matrix(data_in)
   m = nrow(data_in)
  #################  stat of step 1: UNGrouped   ################
  cat("\n","^^^Implement UNGROUPED FDR procedures ...^^^","\n")
  UngroupedResults = GeneralizedFDREstimators(data=data_in, Test=test_in,
                                  FET_via = FET_via_in,lowerTail = lowerTail_in, FDRlevel=FDRlevel_in,TuningRange = Tunings) 
                                          
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
        
        ## allow reference divergence
        if (RefDivergence == "Yes") { 
             
             cat("^^Grouping by divergence: compute supremum norm w.r.t. Unif ...","\n")
             d_supp = Div_Ref_Unif(pval_Supp_ungrped,test_in)
             du_vec = d_supp[[1]]; ss_vec = d_supp[[2]]
             
             du_max = max(du_vec); du_min = min(du_vec)  #minimal is always zero
             ss_max = max(ss_vec); ss_min = min(ss_vec)
      
             cat("Max and Min RefNorm:", du_max,",", du_min, "; Max and Min supp size:",ss_max,",",ss_min,"\n")
             
            # div_ref = du_vec/du_max + ss_vec/ss_max  # option A # relaxed triangular inequality
             
             div_ref = du_vec/du_max   # option B #relaxed triangular inequality
             
             div_ref_max = max(div_ref) ;  div_ref_min = min(div_ref); diff_div = div_ref_max - div_ref_min
             
             div_brks = seq(from=div_ref_min, to=div_ref_max, by=diff_div/ngrp)
           } # end of reference divergence
         
         ## if not reference divergence  
         if (RefDivergence == "No") {
             
             cat("^^Grouping by pairwises divergences. Be patient ...","\n")

             # insert:  when  |F_i - id| <= m^-3, take it as uniform   
              pv_supp_And_id_appr_unif = Div_Appr_Unif(pval_Supp_ungrped,test_in, unif_tol)  #(pvSpList=NULL,test_used=NULL,appr_unif_tol = NULL) 
              pv_supp_formatted = pv_supp_And_id_appr_unif[[1]]
              id_appr_unif = pv_supp_And_id_appr_unif[[2]]
              
              lgX = length(id_appr_unif)
              if (lgX > 0)  { 
                  cat("^^",lgX,"of p-value cdf's are very close to and will be identified as Unif","\n")
                  
                  if (lgX < m) {
                    id_pvDist_far_unif = (1:m)[-id_appr_unif]
                    ng_div =  ngrp-1
                   # grpidx_list[[1]] = id_appr_unif   # first group via divergence
                    
                    lgXa =  length(id_pvDist_far_unif)
                    pv_supp_far_unif = vector("list",lgXa)
                     for (iX in 1:lgXa)  
                      pv_supp_far_unif[[iX]] = pv_supp_formatted[[iX]]
                    }  
               } else { 
                     cat("^^All p-value cdf's are NOT close enough to Unif","\n") 
                     ng_div =  ngrp 
                     pv_supp_far_unif = pval_Supp_ungrped
                     }
             # end insertion  
             
             ## group those that are far from Unif
             scalfac_in = 1
             Diff_mat = GetDivergenceMatrix(scalfac_in,pv_supp_far_unif)  # change scalfac when needed
             Divs_mat = Diff_mat[[1]]; chi_mat = Diff_mat[[2]]; infNorm_mat = Diff_mat[[3]]; chi_max = Diff_mat[[4]]
    
             chi_vec = as.vector(chi_mat[lower.tri(chi_mat,diag=FALSE)])
             infNorm_vec = as.vector(infNorm_mat[lower.tri(infNorm_mat,diag=FALSE)])
             Divs_vec = as.vector(Divs_mat[lower.tri(Divs_mat,diag=FALSE)])
    
             div_max = max(Divs_vec); div_min = 0 #min(Divs_vec)  #minimal is always zero
    
             cat("Max infNorm is", max(infNorm_vec), "; max ChiSD is",chi_max,"; max Div is",div_max,"\n")
    
             # boxplot may cause memory surge
    
             # two cases for eNet size
             if (is.null(eNetSize))  {
               rad = (div_max - div_min)/(2*ng_div) 
               # make ngrp-1 groups out of those cdf's far from Unif  
               grpidx_list_div = eNetFull(Divs_mat, ng_div, 30, rad, 30) #(data, ngrp, merge, rad, grpsize) 
               }  else {
                 grpidx_list_div = eNetFull(Divs_mat,ng_div, 30, eNetSize, 30)  # last arg, minimal gp size, 2nd von last merge size
               }  # end of eNet
    
            # ngp_div = length(grpidx_list_div)
             cat("^^Number of groups for cdf's far from Unif via divergence is:", ng_div,"\n")
           }  # end if not reference divergence  
         
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


         
       ############ start of  grouped data, groupwise anlysis  #######
       grpdata_list[[j]] = cbind(data_in[grpidx_list[[j]],],grpidx_list[[j]])

       # in order to group, original data is added with a col of total counts
       nc_kep = ncol(data_in)

       grpdata_in = grpdata_list[[j]][,1:nc_kep]

       CheckIt(grpdata_in)   # check line

       ### when a group of data is non-empty, analyze them
       if (nrow(grpdata_in) == 0) {
          stop("^^^Group ",j," has no data, please adjust the number of group to split data into...","\n")
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
