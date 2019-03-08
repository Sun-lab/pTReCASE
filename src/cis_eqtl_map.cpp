#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <time.h>
#include <string.h>
#include <R_ext/Applic.h>
#include <RcppEigen.h>

#define MIN(x, y) (((x) > (y)) ? (y) : (x))
#define ABS(x) (x>0 ? (x) : (-x))

RcppExport SEXP cis_eqtl_map(SEXP SNP_Loc_, SEXP genes_Start_, SEXP genes_End_,
							 SEXP cis_window_, SEXP n_genes_, SEXP n_snps_, SEXP MASK_out_){
	// RcppEigen Declarations
	Rcpp::IntegerVector SNP_Loc(SNP_Loc_);
	Rcpp::IntegerVector genes_Start(genes_Start_);
	Rcpp::IntegerVector genes_End(genes_End_);
	Rcpp::IntegerMatrix MASK_out(MASK_out_);
	
	int cis_window = Rcpp::as<int>(cis_window_);
	int n_genes    = Rcpp::as<int>(n_genes_);
	int n_snps     = Rcpp::as<int>(n_snps_);
	
	int endl = 0, endr = 0;
	int n_tests = 0;
	
	int i=0, j = 0;
	for(j=0;j<n_snps;j++){
		for(i=0;i<n_genes;i++){
      endl = genes_Start[i]-cis_window;
		  endr = genes_End[i]+cis_window;
		  if((SNP_Loc[j]<=endr)&&(SNP_Loc[j]>=endl)){
				MASK_out(i,j) = 1;
				n_tests += 1;
			} else {
				MASK_out(i,j) = 0;
			}
		}
	}
	Rprintf("No. Tests: %d\n",n_tests);
	return Rcpp::wrap(MASK_out);
}

RcppExport SEXP closest_snp(SEXP SNP_Loc_, SEXP genes_Start_, SEXP genes_End_,SEXP n_genes_,
                            SEXP n_snps_){
  //RcppEigen Declarations
  Rcpp::IntegerVector SNP_Loc(SNP_Loc_);
  Rcpp::IntegerVector genes_Start(genes_Start_);
  Rcpp::IntegerVector genes_End(genes_End_);
  
  int n_genes    = Rcpp::as<int>(n_genes_);
  int n_snps     = Rcpp::as<int>(n_snps_);
  
  // New Vector:
  Rcpp::IntegerVector snp_ind(n_genes);
  
  // Determining the closest SNP
  int i=0, j=0, distl = 0, distr = 0, dist=0, dist_min = 0;
  for(i=0;i<n_genes;i++){
    snp_ind[i] = 0;
    for(j=0;j<n_snps;j++){
      distl = ((genes_Start[i]-SNP_Loc[j])>=0) ? (genes_Start[i]-SNP_Loc[j]):(SNP_Loc[j]-genes_Start[i]);
      distr = ((genes_End[i]-SNP_Loc[j])>=0) ? (genes_End[i]-SNP_Loc[j]):(SNP_Loc[j]-genes_End[i]);
      dist  = (distl>distr) ? (distr):(distl);
      if((j==0)||(dist<dist_min)){
        snp_ind[i] = j+1;
        dist_min = dist;
      } else {
        
      }
    }
  }
  
  return Rcpp::wrap(snp_ind);
}
