pTReCASE_multapp<-function(Y, Y1, Y2, Z, X, rhos, F_out, geno_pos, 
                           gene_start, gene_end, min_ASE_Total=8, 
                           min_nASE=10, eps=1e-5, convg=1e-5, 
                           maxIter = 200, maxIter2 = 200, 
                           Perform_CT_co = 1e-3, CT_pco = 1e-3, 
                           Obs_Score = 1, useASE = 1, 
                           cis_window = 100000, gene_add=0){
  #-------------------------------------------------------------------#
  # Specifying Placeholder Values for Older options                   #
  #-------------------------------------------------------------------#
  # Don't allow users to change these deprecated options but still 
  # passes them to the function where there are disregarded.
  lik_convg   = 1e-20
  Testing_Ind = 0
  SNP_Names   = NULL
  Gene_Names  = NULL
  
  #-------------------------------------------------------------------#
  # check input                                                    #
  #-------------------------------------------------------------------#
  
  if(any(is.na(Y))){
    stop("Y has missing values")
  }
  
  if(any(is.na(Y1))){
    stop("Y1 has missing values")
  }
  
  if(any(is.na(Y2))){
    stop("Y2 has missing values")
  }
  
  if(any(is.na(Z))){
    stop("Z has missing values")
  }
  
  if(any(is.na(X))){
    stop("X has missing values")
  }
  
  if(any(is.na(rhos))){
    stop("rhos has missing values")
  }
  
  if(nrow(Y)!=nrow(Z)){
    stop("y and z have different lengths")
  }
  
  
  if(nrow(Y)!=nrow(Z)){
    stop("y and z have different lengths")
  }
  
  if(nrow(Y)!=length(rhos)){
    stop("y and rhos have different lengths.\n")
  }
  
  if(nrow(Y)!=nrow(Y1)||ncol(Y)!=ncol(Y1)){
    stop("y and y1 have different dimensions.\n")
  }
  
  if(nrow(Y1)!=nrow(Y2)||ncol(Y1)!=ncol(Y2)){
    stop("y1 and y2 have different dimensions.\n")
  }
  
  if(length(geno_pos)!=ncol(Z)){
    stop("geno_pos and Z have different dimensions")
  }
  
  if(length(gene_start)!=length(gene_end)){
    stop("Gene Start and Gene Ends should reflect same number of genes")
  }
  
  if(length(gene_start)!=ncol(Y)){
    stop("Gene Endpoints and No. of Genes should be Identical!")
  }
  
  if(!all(Z %in% c(0,1,3,4))){
    stop("z must have values of 0, 1, 3, or 4 \n")
  }
  
  if(!all(Y>=0)){
    stop("Y must be non-negative")
  }
  
  if(!all(Y>=0)){
    stop("Y1 must be non-negative")
  }
  
  if(!all(Y2>=0)){
    stop("Y2 must be non-negative")
  }
  
  if(!all((Y1+Y2)<=Y)){
    stop("Y1+Y2 must be less than or equal to Y!")
  }
  
  if(Testing_Ind == 2 && is.null(F_out)){
    F_out  = c("TReCASE_Compout.txt", "Place_Holder")
  }
  
  if(is.null(SNP_Names)){
    SNP_Names = as.character(c(1:ncol(Z)))
  }
  
  if(is.null(Gene_Names)){
    Gene_Names = as.character(c(1:ncol(Y)))
  }
  
  #-------------------------------------------------------------------#
  # APPLY THE FUNCTION                                                #
  #-------------------------------------------------------------------#
  maxAS   = max(Y1+Y2)
  ctvec   = c(0:maxAS)
  rctvec  = c(maxAS:0)
  lctvec  = lgamma(ctvec+1)
  lrctvec = lgamma(rctvec+1)

  outList = .Call("TReCASE_mtest_only", Y, Y1, Y2, Z, X, rhos, 
                  min_ASE_Total, min_nASE,eps, convg, lik_convg,
                  maxIter, maxIter2, F_out, SNP_Names, Gene_Names,
                  ctvec, lctvec, rctvec, lrctvec, maxAS, 
                  Perform_CT_co, CT_pco, Obs_Score, useASE, 
                  geno_pos, gene_start, gene_end, cis_window, 
                  gene_add)
  
  return(outList)
}
