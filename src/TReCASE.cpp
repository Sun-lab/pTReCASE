#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <time.h>
#include <string.h>
#include <RcppEigen.h>
#include "TReCASE.h"
#include "TReC.h"
#include "ASE.h"
#include "glm.h"
#include "sort_func.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define ABS(x) (x>0 ? (x) : (-x))
#define IDX_i 2
#define IDX_j 1
//#define DEBUG_drop 1

// TReC only declaration
int TReC_sfit(double* para0, double* ex1, double* ex2, double* X, int M,
              int maxIter2, double* rhos_ASall,
              int maxIter, double eps, double convg, double lik_convg,
              double* lower_psi, double* upper_psi, int* nbd_psi,
              double* curr_parap, double* lower_keg, double* upper_keg,
              int* nbd_keg, double* curr_parak, double* init_parak,
              double* Beta, double* Beta_old, double* offset, int* Yint,
              double* wa, int* iwa, double* g1, SEXP x1,
              double* phi_out, double* psi_out, double* LLik_out,
              int* fail_out, int* iter_out, double* final_duration,
              double* KEG_duration, double* NB_duration, double * psi_duration,
              double* gv, double* tv, double* Xv, double* cv, double** Bv, int* mask_v);

void CisTrans_Score(double* score_val, double* pval, double* ScoreVec_,
                    double* mu_, double* musq_,double* Delta1_,double* Delta2_,double* Delta3_,
                    double* Ibb_, double* Ibe_, double* Iee_, double* Iep_, double* Iea_,
                    double* Ipp_, double* Ipa_, double* Iaa_, double* M1_, double* M2_,  
                    double* ex1, double* ex2, double* rhosAS, int M, int n, double* Xmat_,
                    double* ctvec, double* lgct, double* rct_vec, double* lgrct, double* pvec, double* tmpctvec, double* tmprctvec,
                    double* Dmu_, double* deps_, double* offset_, double* Beta_,
                    double* outvec, double psi, double phi, double kappa, double eta, double gamma, int maxAS);
 
 void CisTrans_ObsScore(double* score_val, double* pval, double* ScoreVec_,
                        double* mu_, double* musq_,double* Delta1_,double* Delta2_,double* Delta3_,
                        double* Ibb_, double* Ibe_, double* Iee_, double* Iep_, double* Iea_,
                        double* Ipp_, double* Ipa_, double* Iaa_, double* M1_, double* M2_,  
                        double* ex1, double* ex2, double* rhosAS, int M, int n, double* Xmat_,
                        double* ctvec, double* lgct, double* rct_vec, double* lgrct, double* pvec, double* tmpctvec, double* tmprctvec,
                        double* Dmu_, double* deps_, double* offset_, double* Beta_,
                        double* outvec, double psi, double phi, double kappa, double eta, double gamma, int maxAS,
                        double* Ibt_, double* Iet_, double* Itt_); 
  
/*** Necessary for the Inclusion of the VMMIN backup routine ***/
static double ** Lmatrix(int n)
{
  int   i;
  double **m;
  
  m = (double **) R_alloc(n, sizeof(double *));
  for (i = 0; i < n; i++)
    m[i] = (double *) R_alloc((i + 1), sizeof(double));
  return m;
}

/*****************************************************************
 * TReC Likelihoods and Gradients
 *
 * All the code necessary to compute the likelihood and gradient 
 * of the TReC component of the pur.TReCASE model. Three functions
 * are contained:
 *     (1)
 *     (2)
 *     (3)
 *
 * INPUTS:
 *  (1)
 *       -
 *  (2)
 *       -
 *  (3)
 *       -
 *****************************************************************/
// para 0 : initial estimates of kappa, eta and gamma in that order
// ex1: Extra parameters for optimization of KEG (stored as in dj_sort.cpp / n_subj, n_AS_HetZ, n0,n1,n2)
// ex2: Extra parameters for optimization of psi (stored nAS/n0/n1/n2/vector nT/vector nB/vector Pis - the probability of an AS read mapping to an allele )
// X: Double pointer to matrix of covariates for NB regression for each subject
// M: The number of covariates utilized in the NB regression
// maxIter: Number of iterations the TReCASE iterative loop runs
// maxIter2: Number of iterations Negative Binomial/ Poisson Fits run
// rhos_ASall: The tumor purities for all subjects with sufficient AS reads
// eps: G.O.F. value for NB regression (if weights are less than a value then there are fit issues)
// Convg: Value that specifies whether or trecase model has converged
// lik_convg: How little does the likelihood have to change between updates to be converged (flatness)
// lower_: Lower bounds for paramters
// upper_: upper bounds for parameters
// nbd_:  The number of bounds for parameters (0 - no bounds, 1- lower bound only, 2- lower + upper bound, 3- upper only)
// curr_para: Current parameter values to use in optimization
// init_para: initial parameters for use in optimization
// Beta: Regression coefficients for NB
// Beta old: Vector for storing old coefficients
// offset: vector storing regression offset for each individual
// Yint: Integer vector used in computing the negative binomial likelihood at the MLE
// wa: vector created for lbfgsb1 to avoid garbage collection issues in R
// iwa: same purpose as above
// g1: same purpose as above
// x1: same purpose as above

// remaining variables are output variables

int TReCASE_sfit(double* para0, double* ex1, double* ex2, double* X, int M,
            		 int maxIter2, double* rhos_ASall,
						     int maxIter, double eps, double convg, double lik_convg,
						     double* lower_psi, double* upper_psi, int* nbd_psi,
						     double* curr_parap, double* lower_keg, double* upper_keg,
						     int* nbd_keg, double* curr_parak, double* init_parak,
						     double* Beta, double* Beta_old, double* offset, int* Yint,
						     double* wa, int* iwa, double* g1, SEXP x1,
						     double* phi_out, double* psi_out, double* LLik_out,
						     int* fail_out, int* iter_out, double* final_duration,
						     double* KEG_duration, double* NB_duration, double* Psi_duration,
						     double* gv, double* tv, double* Xv, double* cv, double** Bv, int* mask_v){
  // Starting the Timer 
  clock_t tstart = clock();
    
  // Setting up Data-specific Parameters:
  int n_subj = 0, nAS_Hz = 0, H0 = 0;
  
  n_subj  = ceil(ex1[0]-0.5);
  nAS_Hz  = ceil(ex2[1]-0.5)+ceil(ex2[3]-0.5);
  //n0      = ceil(ex1[2]-0.5);
  //n1      = ceil(ex1[3]-0.5);
  
  H0 = ceil(ex1[5]-0.5);
  
  // L-BFGS-B declarations
  /*** Non-Unique Parameters ***/
  int lmm, fail, fncount, grcount, maxit1, nREPORT, fail_stat=-1;
  double Fmin = 0.0, factr = 0.0, pgtol = 0.0;
  
  lmm     = 5; // Number of terms involved in hessian approximation
  fail    = 0;
  fncount = 0;
  grcount = 0;
  maxit1  = 200;
  nREPORT = 5;
  
  factr   = 1e7;
  pgtol   = 0.0;
  
  char msg[1023];
  
  /*** Psi Optimization Parameters ***/
  int npara_p  = 1;
  double psi_curr = 0.5,
         psi_old  = 0.5;
  
  lower_psi[0] = -50;
  upper_psi[0] = 100;
  
  lower_psi[1] = -50;
  upper_psi[1] = 100;
  
  nbd_psi[0] = 0;
  nbd_psi[1] = 0;
  
  curr_parap[0] = log(0.5);
  curr_parap[1] = log(0.5);
  
  //Rprintf("Initialized Curr_Parap\n");
  
  /*** KEG Optimization Parameters ***/
  int npara_keg = 3;
  double kappa = 1.0, eta = 1.0, gamma = 1.0;
  
  lower_keg[0] = -50;
  lower_keg[1] = -50;
  lower_keg[2] = -50;
  
  upper_keg[0] = 50;
  upper_keg[1] = 50;
  upper_keg[2] = 50;
  
  nbd_keg[0] = 0;
  nbd_keg[1] = 0;
  nbd_keg[2] = 0;
  
  //Rprintf("Initialized L,U,NBD for KEG\n");
  
  curr_parak[0] = log(para0[0]);
  kappa = para0[0];
  if(H0==0){
    eta   = para0[1];
    gamma = para0[2];
    
    curr_parak[1] = log(eta);
    curr_parak[2] = log(gamma);
  } else if(H0==1){
    npara_keg = 2;
    eta       = 1.0;
    gamma     = para0[1];
    
    curr_parak[1] = log(gamma);
  } else if(H0==2){
    npara_keg = 2;
    eta       = para0[1];
    gamma     = 1.0;
    
    curr_parak[1] = log(eta);
  }
  
  double abstol = 0.0, reltol = 0.0;
  abstol = 1e-8;
  reltol = 1e-8;  
  
  //Rprintf("TReCASE_Sfit (KEG): %d, %d, %d\n",kappa,eta,gamma);
  //Rprintf("H0(%d) // Kappa // EoG // GoN: %f, %f, %f\n",H0,curr_parak[0],curr_parak[1],curr_parak[2]);
  
  //clock_t init_e  = clock();
  //double init_dur = ((double)(init_e - tstart))/CLOCKS_PER_SEC; 
  
  /**************************************************************************************
  * Initialize Beta and Phi
  **************************************************************************************/
  clock_t NB_s, NB_e;
  
  NB_s = clock();
  
  int i=0;
  double phi_curr = 0.5, phi_old = 0.5, logLNB = 0.0;
  
  for(i=0;i<n_subj;i++){
    Yint[i] = ceil(ex1[(9+i)]-0.5);
  }
  
  //Rprintf("Yint[last] = %d\n",Yint[(n_subj-1)]);
  
  /*** Compute the Offset ***/
  comp_offset(kappa,eta,gamma,ex1,offset);
  
  /*** Fit Poisson as Starter ***/
  IRLS_Pois_fit((ex1+9),X,n_subj,M,offset,
                maxIter2,eps,convg,Beta);
  
  //Rprintf("AIP // Beta Vals: %f, %f, %f\n",Beta[0],Beta[1],Beta[2]);
  
  phi_initialize((ex1+9),X,n_subj,M,offset,
                 Beta,&phi_curr);
  
  //Rprintf("Initial Phi: %f\n",phi_curr);
  
  /*** Fit The Negative Binomial Model ***/
  IRLS_NB_fit((ex1+9),Yint,X,n_subj,M,offset,
              maxIter2,eps,convg,0,&phi_curr,&logLNB,Beta);
  
  //Rprintf("AINB // Beta Vals: %f, %f, %f, %f\n",Beta[0],Beta[1],Beta[2]);
  
  // Will initialize nu0 in the loop to minimize the number of times its computed
  
  ex1[7]  = phi_curr;
  phi_old = phi_curr;
  
  NB_e = clock();
  double NB_dur = ((double)(NB_e-NB_s))/CLOCKS_PER_SEC;
  
  //Rprintf("Initial Phi: %.7f\n",phi_curr);
  
  /**************************************************************************************
  *  Initialize Psi
  **************************************************************************************/
  clock_t Psi_s, Psi_e;
  Psi_s = clock();
  
  /*** Compute Pis ***/
  compute_pis(kappa,eta,gamma,rhos_ASall,ex2);
  
  //Rprintf("Past Compute Pis!\n");
  
  /*** Update optimization parameters ***/
  fail = 0;
  
  psi_optroutine(npara_p,lmm,curr_parap,lower_psi,upper_psi,
                 nbd_psi,&Fmin,&fail,ex2,factr,pgtol,&fncount,&grcount,
                 maxit1,msg,nREPORT,wa,iwa,g1,x1);
  
  Psi_e = clock();
  double Psi_dur = ((double)(Psi_e-Psi_s))/CLOCKS_PER_SEC;
  
  //Rprintf("After psi_optroutine: %f // Fail Stat: %d\n",curr_parap[0],fail);
  
  if(fail==100){
    fail_stat=100;
    *psi_out   = psi_curr;
    *phi_out = phi_curr;
    *LLik_out   = -1e10;
    *fail_out = fail_stat;
    *iter_out      = 0;
    *final_duration  = 0;
    *KEG_duration   = 0;
    *NB_duration    = 0;
    *Psi_duration  = 0;
    return(fail_stat);
  }
  
  /*** Update the ex1 parameter ***/
  ex1[8] = exp(curr_parap[0]);
  
  psi_curr = ex1[8];
  
  /**************************************************************************************
  * Fit Loop
  **************************************************************************************/	
  int g=0, h=0;
  double logLik0 = -1e10, logLik1 = -1e10, pdTmp = 0.0, parDiff = 0.0;
  double diff_calc = 0.0, KEG_dur = 0.0, abs_0 = 0.01;
  clock_t KEG_s, KEG_e;
  for(g=0;g<maxIter;g++){
    //if(g%50==0){
    //  Rprintf("Inside Loop!\n");
    //}
    parDiff = 0.0;
    //ologLik = logLik0;
    /********************************************************************************* 
    * Fit KEG
    *********************************************************************************/
    /*** Time the KEG Procedure ***/
    KEG_s = clock();
    
    /*** Prep the Necessities ***/
    fail    = 0;
    fncount = 0;
    grcount = 0;
    
    for(h=0;h<npara_keg;h++){
      init_parak[h] = curr_parak[h];
    }
    
    compute_nu0(X,n_subj,M,Beta,ex1);
    
    /*** Fit using L-BFGS-B ***/
    #ifdef DEBUG_drop
      Rprintf("Iter (%d): Loglik0 %.7f\n",g,logLik0);
    #endif
      
    KEG_optroutine(npara_keg,lmm,curr_parak,lower_keg,upper_keg,
                   nbd_keg,&Fmin,&fail,ex1,factr,pgtol,
                   &fncount,&grcount,maxit1,msg,nREPORT,wa,iwa,g1,x1);
    
    if(fail>1){
      // KEG are the horses of this routine, we should attempt to improve
      // fit in case of failure. Move me to a position where I won't fail!
      
      for(h=0;h<npara_keg;h++){
        curr_parak[h] = init_parak[h];
      }
      
      fail = 0;
      KEG_boptroutine(npara_keg, curr_parak, &Fmin, maxit1,mask_v,abstol,reltol,
                      nREPORT,ex1,&fncount,&grcount,&fail,gv,tv,Xv,cv,Bv,x1);
    }
    
    /*** Update Timer ***/
    KEG_e = clock();
    KEG_dur += ((double)(KEG_e-KEG_s))/CLOCKS_PER_SEC;
    
    /*** Update and Compare Likelihood ***/
    // Re-Establish Kappa, Eta, and Gamma:
    if(fail>1){
      //Rprintf("Failed to Fit (%d) K//E//G!\n",g);
      // Do not update the parameters in the event of failure
      
      if(fail==100){
        fail_stat=100;
        *psi_out   = psi_curr;
        *phi_out   = phi_curr;
        *LLik_out   = logLik1;
        *fail_out = fail_stat;
        *iter_out      = g;
        *final_duration  = 0;
        *KEG_duration   = 0;
        *NB_duration    = 0;
        *Psi_duration  = 0;
        return(fail_stat);
      }
      
      for(h=0;h<npara_keg;h++){
        curr_parak[h] = init_parak[h];
      }
      
      // Keep the likelihood the same
      logLik1 = logLik0;
      
      pdTmp = 0.0;
      
      fail_stat = 3;
      
    } else {
      if(H0 == 0){
        kappa = exp(curr_parak[0]);
        eta   = exp(curr_parak[1]);
        gamma = exp(curr_parak[2]);
      } else if(H0 == 1){
        kappa = exp(curr_parak[0]);
        gamma = exp(curr_parak[1]);
      } else if(H0 == 2){
        kappa = exp(curr_parak[0]);
        eta   = exp(curr_parak[1]);
      }
      
      //if(g % 50 == 0){
      //Rprintf("Iter (%d) -- K//E//G: %.7f, %.7f, %.7f\n",g,kappa,eta,gamma);
      //}
      
      // Update the Likelihood
      logLik1 = -Fmin;
      logLik1 += compute_BBHz(nAS_Hz,ex2,psi_curr);
      
      pdTmp = 0.0;
      
      for(h=0;h<npara_keg;h++){
        diff_calc = exp(curr_parak[h])-exp(init_parak[h]);
        diff_calc = (diff_calc>=0) ? (diff_calc):(-diff_calc);
        pdTmp = (diff_calc>pdTmp) ? diff_calc:pdTmp;
      }
    }
    
    abs_0 = (logLik0>=0) ? (logLik0):(-logLik0);
    
    if(g>1){
      if(((logLik1 - logLik0)/abs_0)<-1e-5){
        Rprintf("Likelihood Decreases for K.E.G.\n");
        fail_stat = 1;
        break;
      }
    }
    
  #ifdef DEBUG_drop
    Rprintf("Loglik1 (After KEG):%.7f\n",logLik1);
  #endif
    
    logLik0 = logLik1;
    
    if(parDiff<pdTmp){
      parDiff = pdTmp;
    }
    
    /********************************************************************************* 
    * Re-fit Beta and Phi
    *********************************************************************************/
    /*** Time the NB Event  ***/
    NB_s = clock();
    
    /*** Compute the Offset ***/
    comp_offset(kappa,eta,gamma,ex1,offset);
    
    /*** Fit The Negative Binomial Model ***/
    phi_old = phi_curr;
    IRLS_NB_fit((ex1+9),Yint,X,n_subj,M,offset,
                maxIter2,eps,convg,0,&phi_curr,&logLNB,Beta);
    
    /*** Time the NB event ***/
    NB_e = clock();
    NB_dur += ((double)(NB_e-NB_s))/CLOCKS_PER_SEC;
    
    /*** Updating the convergence Parameters ***/
    logLik1 = logLNB;
    pdTmp = ((phi_curr-phi_old)>=0)? (phi_curr-phi_old):(phi_old-phi_curr);
    
    #ifdef DEBUG_drop
      /*** Compute Pis ***/
      compute_pis(kappa,eta,gamma,rhos_ASall,ex2);
      Rprintf("Loglik1 (After NB):%.7f\n",logLik1-logLASE_psi(1,&psi_curr,ex2));
    #endif
    
    if(parDiff<pdTmp){
      parDiff = pdTmp;
    }
    
    //Rprintf("Phi parDiff: %.7f\n",parDiff);
    
    ex1[7] = phi_curr;
    
    //	Rprintf("Current Phi (%d): %.7f\n",g,phi_curr);
    
    /********************************************************************************* 
    * Re-fit psi
    *********************************************************************************/
    /*** Time the Psi ***/
    Psi_s = clock();
    
    /*** Compute Pis ***/
    compute_pis(kappa,eta,gamma,rhos_ASall,ex2);
    
    /*** Update optimization parameters ***/
    fail = 0;
    
    psi_old = psi_curr;
    
    /*** Fit Routine ***/
    psi_optroutine(npara_p,lmm,curr_parap,lower_psi,upper_psi,
                   nbd_psi,&Fmin,&fail,ex2,factr,pgtol,&fncount,&grcount,
                   maxit1,msg,nREPORT,wa,iwa,g1,x1);
    
    if(fail>1){
      // Psi may be the cause of some failures. Offer a backup routine in the event of failure.
      fail = 0;
      psi_boptroutine(npara_p, curr_parap, &Fmin, maxit1,mask_v,abstol,reltol,
                      nREPORT,ex2,&fncount,&grcount,&fail,gv,tv,Xv,cv,Bv,x1);
    }
    
    /*** Time the Psi ***/
    Psi_e = clock();
    Psi_dur += ((double)(Psi_e-Psi_s))/CLOCKS_PER_SEC;
    
    /*** Checking Fail Status ***/
    
    if(fail>1){
      //Rprintf("Failed to Fit (%d) Psi!\n",g);
      // Do not allow update if fit fails
      
      if(fail==100){
        fail_stat=100;
        *psi_out   = psi_curr;
        *phi_out = phi_curr;
        *LLik_out   = logLik1;
        *fail_out = fail_stat;
        *iter_out      = g;
        *final_duration  = 0;
        *KEG_duration   = 0;
        *NB_duration    = 0;
        *Psi_duration  = 0;
        return(fail_stat);
      }
      
      Rprintf("Psi Failed!\n");
      
      curr_parap[0] = log(psi_old);
      
      // Keep the likelihood the same 
      logLik1 += -logLASE_psi(1,curr_parap,ex2,x1);
    } else {
      psi_curr = exp(curr_parap[0]);
      ex1[8]   = psi_curr;
      
      // Update the likelihood:
      logLik1 += -Fmin;
    }
    
    pdTmp = ((psi_curr-psi_old)>=0) ? (psi_curr-psi_old):(psi_old-psi_curr);
    
    #ifdef DEBUG_drop
      Rprintf("Loglik1 (After Psi):%.7f\n",logLik1);
    #endif
    
    if(parDiff<pdTmp){
      parDiff = pdTmp;
    }
    
    //Rprintf("Psi parDiff: %.7f\n",parDiff);

    /********************************************************************************* 
    * Check Likelihood and Convergence
    *********************************************************************************/
    abs_0 = (logLik0>=0) ? (logLik0):(-logLik0);
    
    if(((logLik1 - logLik0)/abs_0)<-1e-5){
      Rprintf("Likelihood Decreases for Estimating Beta, Phi, and Psi.\n");
      fail_stat = 1;
      break;
    }
    
    logLik0 = logLik1;
    
    //likDiff = ((logLik1-ologLik)>=0) ? (logLik1-ologLik):(ologLik-logLik1);
    
    //Rprintf("Iter (%d):\n LikDiff: %f \n K: %f\n E: %f\n G: %f\n Phi: %f\n Psi: %f\n parDiff: %f\n", g, likDiff, kappa, eta, gamma, phi_curr, psi_curr, parDiff);
    
    if((parDiff<convg)){
      fail_stat = (fail_stat==3) ? (3):(0);
      break;
    } 
  }
  
  //Rprintf("Outside Loop! K // E // G: %f // %f // %f \n",kappa,eta,gamma);
  //Rprintf("Curr_Parak: %f, %f, %f\n", curr_parak[0],curr_parak[1],curr_parak[2]);
  
  /********************************************************************************* 
  * STORE THE OUTPUT
  *********************************************************************************/
  // Failure Status
  if(fail_stat==-1){
    fail_stat = 2;
  }
  
  // End Timer
  clock_t tot_dur = clock() - tstart;
  double fin_dur = ((double)tot_dur)/CLOCKS_PER_SEC;
  
  // Return Output:
  *psi_out   = psi_curr;
  *phi_out  = phi_curr;
  *LLik_out  = logLik1;
  *fail_out = fail_stat;
  *iter_out = g;
  *final_duration  = fin_dur;
  *KEG_duration   = KEG_dur;
  *NB_duration   = NB_dur;
  *Psi_duration = Psi_dur;
  
  //Rprintf("Made it to the return statement\n");
  
  return(fail_stat);
}  

RcppExport SEXP TReCASE_mtest_only(SEXP Y_, SEXP Y1_, SEXP Y2_, SEXP Z_, SEXP X_, SEXP rhos_,
								   SEXP min_ASE_Total_, SEXP min_nASE_,
								   SEXP eps_, SEXP convg_, SEXP lik_convg_, SEXP maxIter_, SEXP maxIter2_,
								   SEXP F_out_, SEXP SNP_Names_, SEXP Gene_Names_, SEXP ctvec_,
								   SEXP lctvec_, SEXP rctvec_, SEXP lrctvec_, SEXP maxAS_, SEXP Perform_CT_co_, SEXP CT_pco_,
								   SEXP Obs_Score_, SEXP useASE_, SEXP genot_Info_, SEXP gene_start_, SEXP gene_end_, SEXP cis_window_,
								   SEXP gene_add_){
	// RcppEigen Declarations
	Rcpp::IntegerMatrix Y(Y_);
	Rcpp::IntegerMatrix Y1(Y1_);
	Rcpp::IntegerMatrix Y2(Y2_);
	Rcpp::IntegerMatrix Z(Z_);
	Rcpp::NumericMatrix X(X_);
	Rcpp::NumericVector rhos(rhos_);
	
	Rcpp::CharacterVector F_out(F_out_);
	
	Rcpp::CharacterVector SNP_Names(SNP_Names_);
	Rcpp::CharacterVector Gene_Names(Gene_Names_);
	
	Rcpp::IntegerVector genot_Info(genot_Info_);
	Rcpp::IntegerVector gene_start(gene_start_);
	Rcpp::IntegerVector gene_end(gene_end_);
	
	int useASE     = Rcpp::as<int>(useASE_);
	int cis_window = Rcpp::as<int>(cis_window_);
	int gene_add   = Rcpp::as<int>(gene_add_);
	
	// Output File creation:
	int out = 0, k =0;
	
  FILE* fcomp;
  fcomp = fopen(F_out[0],"w");
  
  if(fcomp == NULL){
    Rprintf("Out File unable to be opened. Ensure directory exists!\n");
    out = 1;
    return Rcpp::wrap(out);
  }
	
	// Extract sample data
	int n_subj = X.nrow();
	int M      = X.ncol();
	
	fprintf(fcomp,"SNP\tGene\tTotalDuration\tKEG_Duration\tNB_Duration\tPsi_Duration\tLRT_Eta\tP_Eta\tLRT_Gamma\tP_Gamma\tLik\tLikEta\tLikGamma\tKappa\tEta\tGamma\tPhi\tPsi\tFail_Full\tFail_Full_eta\tFail_full_gamma\tFail_Eta\tFail_Gamma");
	for(k=0;k<M;k++){
		fprintf(fcomp,"\tCoeff_%d",k);
	}
	fprintf(fcomp,"\tCisTrans_Score\tCisTrans_Pval\tCisTrans_Fail\n");
	
	// Iterative Application Controls
	int min_ASE_Total = Rcpp::as<int>(min_ASE_Total_);
	int min_nASE      = Rcpp::as<int>(min_nASE_);
	int n_genes       = Y.ncol();
	int n_snps        = Z.ncol();
	int maxIter       = Rcpp::as<int>(maxIter_);
	int maxIter2      = Rcpp::as<int>(maxIter2_);
	double eps        = Rcpp::as<double>(eps_);
	double convg      = Rcpp::as<double>(convg_);
	double lik_convg  = Rcpp::as<double>(lik_convg_);
	int ctest, endl, endr;
	
	// Allocation for temporary matrices
	int n0 = 0, n1 = 0, n3 = 0, n4 = 0,
		nAS_all = 0, nAS_h = 0,
		len_ex1 = 9+6*n_subj,
		len_ex2 = 4+3*n_subj;
		
	Rcpp::NumericVector Betatmp(M);	
	Rcpp::NumericMatrix Xt(n_subj,M);
	Rcpp::IntegerVector Yt(n_subj);
	Rcpp::IntegerVector Y1t(n_subj);
	Rcpp::IntegerVector Y2t(n_subj);
	Rcpp::NumericVector rhot(n_subj);
	Rcpp::NumericVector rhos_ASall(n_subj);
	
	Rcpp::IntegerVector idx0(n_subj);
	Rcpp::IntegerVector idx1(n_subj);
	Rcpp::IntegerVector idx3(n_subj);
	Rcpp::IntegerVector idx4(n_subj);
	
	Rcpp::NumericVector ex1(len_ex1);
	Rcpp::NumericVector ex2(len_ex2);
	
	Rcpp::NumericVector para0(3);
		para0[0] = 1.0;
		para0[1] = 1.0;
		para0[2] = 1.0;
	
	/*** CisTrans Score Test necessities: ***/ 
	// Vector Creation //
	Rcpp::NumericVector ctvec(ctvec_);
	Rcpp::NumericVector lctvec(lctvec_);
	Rcpp::NumericVector rctvec(rctvec_);
	Rcpp::NumericVector lrctvec(lrctvec_);
	
	double CT_score = 0.0,
	       CT_pval  = 0.0;
	int maxAS = Rcpp::as<int>(maxAS_);
	int Obs_Score = Rcpp::as<int>(Obs_Score_);
	int CT_Fail = 0;
	double Perform_CT_co = Rcpp::as<double>(Perform_CT_co_);
	double CT_pco = Rcpp::as<double>(CT_pco_);
	
	Rcpp::NumericVector outvec(4);
	Rcpp::NumericVector tmpctvec((maxAS+1));
	Rcpp::NumericVector tmprctvec((maxAS+1));
	Rcpp::NumericVector p_vec((maxAS+1));
	Rcpp::NumericVector mu_vec(n_subj);
	Rcpp::NumericVector musq(n_subj);
	
	// Matrix Placeholders //
	// According to RCPP Quick_Ref by Eddelbuettel, these are initialized
	// to 0, so it should not cause a problem for our CT routine.
	Rcpp::NumericVector CT_ScoreVec(2);
	
	Rcpp::NumericMatrix Ibb(M,M);
	Rcpp::NumericMatrix Ibe(M,3);
	Rcpp::NumericMatrix Iee(3,3);
	Rcpp::NumericMatrix Ibt(M,1);
	Rcpp::NumericMatrix Iet(3,1);
	double Itt = 0;
	Rcpp::NumericMatrix Iep(3,1);
	Rcpp::NumericMatrix Iea(3,2);
	double Ipp = 0;
	Rcpp::NumericMatrix Ipa(1,2);
	Rcpp::NumericMatrix Iaa(2,2);
	Rcpp::NumericMatrix M1((M+5),(M+5));
	Rcpp::NumericMatrix M2((M+5),2);
	Rcpp::NumericMatrix Delta1(n_subj,n_subj);
	Rcpp::NumericMatrix Delta2(n_subj,n_subj);
	Rcpp::NumericMatrix Delta3(n_subj,n_subj);
	
	Rcpp::NumericMatrix Dmu(n_subj,3);
	Rcpp::NumericMatrix deps(n_subj,3);
	
	/*** L-BFGS-B 1 Parameters ***/
	int lmm = 5, npara = 3;
	double *wa, *g1;
	int *iwa;
	SEXP x1;
	Rcpp::NumericVector Rcpp_x1(3);
	PROTECT(x1 = Rcpp::wrap(Rcpp_x1));
	//consider replacing with simple Calloc
	wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
	iwa = (int*) R_alloc((3*npara),sizeof(int));
	g1 = (double *)R_alloc(npara, sizeof(double));
	
	/*** Psi Optimization Parameters ***/
	Rcpp::NumericVector lower_psi(2);
	Rcpp::NumericVector upper_psi(2);
	
	Rcpp::IntegerVector nbd_psi(2);
	
	Rcpp::NumericVector curr_parap(2);
	
	/*** KEG Optimization Parameters ***/
	Rcpp::NumericVector lower_keg(3);
	Rcpp::NumericVector upper_keg(3);
	
	Rcpp::IntegerVector nbd_keg(3);
	
	Rcpp::NumericVector curr_parak(3);
	Rcpp::NumericVector curr_parak1a(3);
	Rcpp::NumericVector curr_parak1b(3);
	Rcpp::NumericVector init_parak(3);
	
	// Storage for easier call in Printing
  Rcpp::NumericVector KEG(3);
  
  /*** KEG Backup Optimization Parameters ***/
  Rcpp::IntegerVector mask_v(3);
  mask_v[0] = 1;
  mask_v[1] = 1;
  mask_v[2] = 1;
  
  Rcpp::NumericVector gv(3);
  Rcpp::NumericVector tv(3);
  Rcpp::NumericVector Xv(3);
  Rcpp::NumericVector cv(3);
  
  double** Bv = Lmatrix(3);
  
	/*** Regression Parameters ***/
	Rcpp::NumericVector Beta(M);
	Rcpp::NumericVector Beta1a(M);
	Rcpp::NumericVector Beta1b(M);
	Rcpp::NumericVector Beta_old(M);
	Rcpp::NumericVector offset(n_subj);
	
	Rcpp::IntegerVector Yint(n_subj);
	
	/*** OUTPUT STORAGE	***/
	double phi_out, psi_out, LLik_out, final_duration, KEG_duration,
	       NB_duration, Psi_duration;
	int    fail_out, iter_out;
	
	// Loop
	int i=0,j=0,test_ct=0;
	double lrt_eta = 0.0, lrt_gamma = 0.0, p_eta = 0.0, p_gamma = 0.0;
	Rcpp::List outList(n_snps);
	Rcpp::List tmpList(n_genes);
	Rcpp::List List1;
	Rcpp::List List1a;
	Rcpp::List List1b;
	Rcpp::List List2;
	Rcpp::List List3;
	for(j=0;j<n_snps;j++){
	  Rprintf("SNP  %d Start! --- %d Tests Completed ---\n",(j+1),test_ct);
		for(i=0;i<n_genes;i++){
		  //Rprintf("SNP %d // Gene %d: Start!\n",j,i);
		  // Are gene and SNP close enough:
		  endl = gene_start[i]-cis_window;
		  endr = gene_end[i]+cis_window;
		  
		  if((genot_Info[j]<=endr)&&(genot_Info[j]>=endl)){
		    ctest = 1;
		  } else {
		    ctest = 0;
		  }
		  
			if(ctest==0){
			    // Don't Allow to run
			} else {
			  // Sort the Data
			  DJ_Rcpp_sort(X.begin(),rhos.begin(),(Y.begin()+i*n_subj),(Y1.begin()+i*n_subj),(Y2.begin()+i*n_subj),(Z.begin()+j*n_subj),
                  Xt.begin(),rhot.begin(),Yt.begin(),Y1t.begin(),Y2t.begin(),
                  &n0, &n1, &n3, &n4, &nAS_all, &nAS_h,
                  n_subj,M,min_ASE_Total,
                  idx0.begin(),idx1.begin(),idx3.begin(),idx4.begin());
			  
			  if((nAS_all>=min_nASE)&&(useASE==1)){
			    test_ct += 1;
			    // Set the Extras
			    set_extras(Yt.begin(),Y1t.begin(),Y2t.begin(),rhot.begin(),n0,(n1+n3),
                  n_subj,nAS_all,nAS_h, min_ASE_Total, min_nASE,
                  ex1.begin(), ex2.begin(), rhos_ASall.begin());
			    
			    //Rprintf("Past Set Extras!\n");
			    
			    try{
			      // Fit the Models
			      ex1[5] = 1;
			      TReCASE_sfit(para0.begin(),ex1.begin(),ex2.begin(),Xt.begin(),M,maxIter2,rhos_ASall.begin(),maxIter,eps,convg,lik_convg,
                      lower_psi.begin(),upper_psi.begin(),nbd_psi.begin(),curr_parap.begin(),lower_keg.begin(),upper_keg.begin(),
                      nbd_keg.begin(),curr_parak.begin(),init_parak.begin(),Beta.begin(),Beta_old.begin(),offset.begin(),Yint.begin(),
                      wa,iwa,g1,x1,&phi_out,&psi_out,&LLik_out,&fail_out,&iter_out,&final_duration,
                      &KEG_duration,&NB_duration,&Psi_duration,gv.begin(),tv.begin(),Xv.begin(),cv.begin(),Bv, mask_v.begin());
			      
			      /* Store this output */
			      // Don't want to save curr_parak or betas from this
			      //List2["Param"] =  Rcpp::clone<Rcpp::NumericVector>(curr_parak);
			      //List2["bs"] =  Rcpp::clone<Rcpp::NumericVector>(Beta);
			      List2["Psi"] = psi_out;
			      List2["Phi"] = phi_out;
			      List2["Lik"] = LLik_out;
			      List2["Fail_Stat"] = fail_out;
			      List2["Iter"] = iter_out;
			      List2["Duration"] =  final_duration;
			      List2["KEG_Dur"] =  KEG_duration;
			      List2["NB_Dur"] =  NB_duration;
			      List2["Psi_Dur"] = Psi_duration;
			      
			      if(Rcpp::as<int>(List2["Fail_Stat"])==100){
			        throw std::invalid_argument("Non-finite Optimization in Eta Hypothesis");
			      }
			      
			      ex1[5] = 0;
			      para0[0] = exp(curr_parak[0]);
			      para0[1] = 1.0;
			      para0[2] = exp(curr_parak[1]);
			      TReCASE_sfit(para0.begin(),ex1.begin(),ex2.begin(),Xt.begin(),M,maxIter2,rhos_ASall.begin(),maxIter,eps,convg,lik_convg,
                      lower_psi.begin(),upper_psi.begin(),nbd_psi.begin(),curr_parap.begin(),lower_keg.begin(),upper_keg.begin(),
                      nbd_keg.begin(),curr_parak1a.begin(),init_parak.begin(),Beta1a.begin(),Beta_old.begin(),offset.begin(),Yint.begin(),
                      wa,iwa,g1,x1,&phi_out,&psi_out,&LLik_out,&fail_out,&iter_out,&final_duration,
                      &KEG_duration,&NB_duration,&Psi_duration,gv.begin(),tv.begin(),Xv.begin(),cv.begin(),Bv, mask_v.begin());
			      
			      List1a["Param"] =  curr_parak1a;
			      List1a["bs"] =  Beta1a;
			      List1a["Psi"] = psi_out;
			      List1a["Phi"] = phi_out;
			      List1a["Lik"] = LLik_out;
			      List1a["Fail_Stat"] = fail_out;
			      List1a["Iter"] = iter_out;
			      List1a["Duration"] =  final_duration;
			      List1a["KEG_Dur"] =  KEG_duration;
			      List1a["NB_Dur"] =  NB_duration;
			      List1a["Psi_Dur"] = Psi_duration;
			      
			      //Rprintf("Past Eta \n");
			      
			      ex1[5] = 2;
			      para0[0] = 1.0;
			      para0[1] = 1.0;
			      TReCASE_sfit(para0.begin(),ex1.begin(),ex2.begin(),Xt.begin(),M,maxIter2,rhos_ASall.begin(),maxIter,eps,convg,lik_convg,
                      lower_psi.begin(),upper_psi.begin(),nbd_psi.begin(),curr_parap.begin(),lower_keg.begin(),upper_keg.begin(),
                      nbd_keg.begin(),curr_parak.begin(),init_parak.begin(),Beta.begin(),Beta_old.begin(),offset.begin(),Yint.begin(),
                      wa,iwa,g1,x1,&phi_out,&psi_out,&LLik_out,&fail_out,&iter_out,&final_duration,
                      &KEG_duration,&NB_duration,&Psi_duration,gv.begin(),tv.begin(),Xv.begin(),cv.begin(),Bv, mask_v.begin());
			      
			      //List3["Param"] =  Rcpp::clone<Rcpp::NumericVector>(curr_parak);
			      //List3["bs"] =  Rcpp::clone<Rcpp::NumericVector>(Beta);
			      List3["Psi"] = psi_out;
			      List3["Phi"] = phi_out;
			      List3["Lik"] = LLik_out;
			      List3["Fail_Stat"] = fail_out;
			      List3["Iter"] = iter_out;
			      List3["Duration"] =  final_duration;
			      List3["KEG_Dur"] =  KEG_duration;
			      List3["NB_Dur"] =  NB_duration;
			      List3["Psi_Dur"] = Psi_duration;
			      
			      if(Rcpp::as<int>(List3["Fail_Stat"])==100){
			        throw std::invalid_argument("Non-finite Optimization in Gamma");
			      }
			      
			      ex1[5] = 0;
			      para0[0] = exp(curr_parak[0]);
			      para0[1] = exp(curr_parak[1]);
			      para0[2] = 1.0;
			      TReCASE_sfit(para0.begin(),ex1.begin(),ex2.begin(),Xt.begin(),M,maxIter2,rhos_ASall.begin(),maxIter,eps,convg,lik_convg,
                      lower_psi.begin(),upper_psi.begin(),nbd_psi.begin(),curr_parap.begin(),lower_keg.begin(),upper_keg.begin(),
                      nbd_keg.begin(),curr_parak1b.begin(),init_parak.begin(),Beta1b.begin(),Beta_old.begin(),offset.begin(),Yint.begin(),
                      wa,iwa,g1,x1,&phi_out,&psi_out,&LLik_out,&fail_out,&iter_out,&final_duration,
                      &KEG_duration,&NB_duration,&Psi_duration,gv.begin(),tv.begin(),Xv.begin(),cv.begin(),Bv, mask_v.begin());
			      List1b["Param"] =  curr_parak1b;
			      List1b["bs"] =  Beta1b;
			      List1b["Psi"] = psi_out;
			      List1b["Phi"] = phi_out;
			      List1b["Lik"] = LLik_out;
			      List1b["Fail_Stat"] = fail_out;
			      List1b["Iter"] = iter_out;
			      List1b["Duration"] =  final_duration;
			      List1b["KEG_Dur"] =  KEG_duration;
			      List1b["NB_Dur"] =  NB_duration;
			      List1b["Psi_Dur"] = Psi_duration;
			      
			      if((Rcpp::as<int>(List1a["Fail_Stat"])==100)&&(Rcpp::as<int>(List1b["Fail_Stat"])==100)){
			        throw std::invalid_argument("Non-finite Optimization in Full Model");
			      } else if((Rcpp::as<int>(List1a["Fail_Stat"])==100)&&(Rcpp::as<int>(List1b["Fail_Stat"])!=100)){
			        List1 = List1b;
			      } else if((Rcpp::as<int>(List1a["Fail_Stat"])!=100)&&(Rcpp::as<int>(List1b["Fail_Stat"])==100)){
			        List1 = List1a;
			      } else {
			        if((Rcpp::as<double>(List1a["Lik"])>=Rcpp::as<double>(List1b["Lik"]))){
			          List1 = List1a;
			        } else {
			          List1 = List1b;
			        }
			      }
			      
			      KEG = Rcpp::as<Rcpp::NumericVector>(List1["Param"]);
			      
			      // Conduct the Tests
			      if(Rcpp::as<int>(List1["Fail_Stat"])!=1&&Rcpp::as<int>(List2["Fail_Stat"])!=1){
			        lrt_eta = -2*(Rcpp::as<double>(List2["Lik"])-Rcpp::as<double>(List1["Lik"]));
			        p_eta   = R::pchisq(lrt_eta,1,0,0);
			      } else {
			        p_eta   = -99.00;
			        lrt_eta = -1.0e8;
			      }
			      
			      if(Rcpp::as<int>(List1["Fail_Stat"])!=1&&Rcpp::as<int>(List3["Fail_Stat"])!=1){
			        lrt_gamma = -2*(Rcpp::as<double>(List3["Lik"])-Rcpp::as<double>(List1["Lik"]));
			        p_gamma   = R::pchisq(lrt_gamma,1,0,0);
			      } else {
			        p_gamma   = -99.00;
			        lrt_gamma = -1.0e8;
			      }
			      
			      //Rprintf("To Timing Calc\n");
			      
			      // Compile results of timers
			      double Tot_dur = Rcpp::as<double>(List1a["Duration"])+Rcpp::as<double>(List1b["Duration"])+Rcpp::as<double>(List2["Duration"])+Rcpp::as<double>(List3["Duration"]);
			      double KEG_dur = Rcpp::as<double>(List1a["KEG_Dur"])+Rcpp::as<double>(List1b["KEG_Dur"])+Rcpp::as<double>(List2["KEG_Dur"])+Rcpp::as<double>(List3["KEG_Dur"]);
			      double NB_dur = Rcpp::as<double>(List1a["NB_Dur"])+Rcpp::as<double>(List1b["NB_Dur"])+Rcpp::as<double>(List2["NB_Dur"])+Rcpp::as<double>(List3["NB_Dur"]);
			      double Psi_dur = Rcpp::as<double>(List1a["Psi_Dur"])+Rcpp::as<double>(List1b["Psi_Dur"])+Rcpp::as<double>(List2["Psi_Dur"])+Rcpp::as<double>(List3["Psi_Dur"]);
			      
			      //Rprintf("Past Timing Calc\n");
			      
			      // Store the Output
			      Betatmp = Rcpp::as<Rcpp::NumericVector>(List1["bs"]);
			      fprintf(fcomp,"%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.18f\t%.10f\t%.18f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\t%d\t%d\t%d\t%d",j,(i+gene_add),Tot_dur,KEG_dur,NB_dur,Psi_dur,lrt_eta,p_eta,lrt_gamma,p_gamma,Rcpp::as<double>(List1["Lik"]),Rcpp::as<double>(List2["Lik"]),Rcpp::as<double>(List3["Lik"]),exp(KEG[0]),exp(KEG[1]),exp(KEG[2]),Rcpp::as<double>(List1["Phi"]),Rcpp::as<double>(List1["Psi"]),Rcpp::as<int>(List1["Fail_Stat"]),Rcpp::as<int>(List1a["Fail_Stat"]),Rcpp::as<int>(List1b["Fail_Stat"]),Rcpp::as<int>(List2["Fail_Stat"]),Rcpp::as<int>(List3["Fail_Stat"]));
			      for(k=0;k<M;k++){
			        fprintf(fcomp,"\t%.10f",Betatmp[k]);
			      }
			      
			      // Conduct CisTrans Score Test //
			      /* TWo Dimensional Score test for differential action in 
			      * ASE and TReC Components of our model which could lead
			      * to spurious results in eQTL testing. */
			      if(p_gamma>=0&&p_gamma<=Perform_CT_co){
			        if(Obs_Score==1){
			          CT_Fail = 0;
			          CisTrans_ObsScore(&CT_score,&CT_pval,CT_ScoreVec.begin(),
                               mu_vec.begin(),musq.begin(),Delta1.begin(),Delta2.begin(),Delta3.begin(),
                               Ibb.begin(),Ibe.begin(),Iee.begin(),Iep.begin(),Iea.begin(),
                               &Ipp,Ipa.begin(),Iaa.begin(),M1.begin(),M2.begin(),
                               ex1.begin(),ex2.begin(),rhos_ASall.begin(),M,n_subj,Xt.begin(),
                               ctvec.begin(),lctvec.begin(),rctvec.begin(),lrctvec.begin(),p_vec.begin(),tmpctvec.begin(),tmprctvec.begin(),
                               Dmu.begin(), deps.begin(),offset.begin(),Betatmp.begin(),
                               outvec.begin(),Rcpp::as<double>(List1["Psi"]),Rcpp::as<double>(List1["Phi"]),
                               exp(KEG[0]),exp(KEG[1]),exp(KEG[2]),maxAS,
                               Ibt.begin(),Iet.begin(),&Itt);
			        }
			        
			        if(Obs_Score==0||(Obs_Score==1&&CT_score<0)){
			          CT_Fail = 1;
			          CisTrans_Score(&CT_score,&CT_pval,CT_ScoreVec.begin(),
                            mu_vec.begin(),musq.begin(),Delta1.begin(),Delta2.begin(),Delta3.begin(),
                            Ibb.begin(),Ibe.begin(),Iee.begin(),Iep.begin(),Iea.begin(),
                            &Ipp,Ipa.begin(),Iaa.begin(),M1.begin(),M2.begin(),
                            ex1.begin(),ex2.begin(),rhos_ASall.begin(),M,n_subj,Xt.begin(),
                            ctvec.begin(),lctvec.begin(),rctvec.begin(),lrctvec.begin(),p_vec.begin(),tmpctvec.begin(),tmprctvec.begin(),
                            Dmu.begin(), deps.begin(),offset.begin(),Betatmp.begin(),
                            outvec.begin(),Rcpp::as<double>(List1["Psi"]),Rcpp::as<double>(List1["Phi"]),
                            exp(KEG[0]),exp(KEG[1]),exp(KEG[2]),maxAS);
			        }
			        
			        if(CT_score<0){
			          CT_Fail = 2;
			        }
			      } else {
			        CT_score = -1e8;
			        CT_pval = 1.0;
			        CT_Fail = 5;
			      }
			      
			      // Storing CisTrans Score Test Values //
			      fprintf(fcomp,"\t%.10f\t%.10f\t%d\n",CT_score,CT_pval,CT_Fail);
			      
			      //Rprintf("Past Storing Output\n");
			    } catch(std::invalid_argument& e){
			      Rprintf("Error Caught!\n");
			      fprintf(fcomp,"%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\t%d\t%d\t%d\t%d",j,(i+gene_add),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,4,4,4,4,4);
			      for(k=0;k<M;k++){
			        fprintf(fcomp,"\t%.10f",0.0);
			      }
			      fprintf(fcomp,"\t%.10f\t%.10f\t%d\n",-1e8,1.0,5);
			    }
			    
			    // Parameter reset
			    para0[0] = 1.0;
			    para0[1] = 1.0;
			    para0[2] = 1.0;
			    
			    //Rprintf("Past Parameter Reset\n");
			  }
			  
  			if((nAS_all<min_nASE)||useASE==0||(CT_pval<CT_pco)||(CT_Fail==2)){
  			  //Rprintf("We are inside the TReC only piece! \n");
  			  // TReC only
  			  test_ct += 1;
  			  // Set the Extras
  			  set_extras(Yt.begin(),Y1t.begin(),Y2t.begin(),rhot.begin(),n0,(n1+n3),
                  n_subj,nAS_all,nAS_h, min_ASE_Total, min_nASE,
                  ex1.begin(), ex2.begin(), rhos_ASall.begin());
  			  
  			  //Rprintf("Past Set Extras!\n");
  			  
  			  try{
  			    // Fit the Models
  			    ex1[5] = 1;
  			    TReC_sfit(para0.begin(),ex1.begin(),ex2.begin(),Xt.begin(),M,maxIter2,rhos_ASall.begin(),maxIter,eps,convg,lik_convg,
                      lower_psi.begin(),upper_psi.begin(),nbd_psi.begin(),curr_parap.begin(),lower_keg.begin(),upper_keg.begin(),
                      nbd_keg.begin(),curr_parak.begin(),init_parak.begin(),Beta.begin(),Beta_old.begin(),offset.begin(),Yint.begin(),
                      wa,iwa,g1,x1,&phi_out,&psi_out,&LLik_out,&fail_out,&iter_out,&final_duration,
                      &KEG_duration,&NB_duration,&Psi_duration,gv.begin(),tv.begin(),Xv.begin(),cv.begin(),Bv, mask_v.begin());
  			    
  			    /* Store this output */
  			    // Don't want to save curr_parak or betas from this
  			    //List2["Param"] =  Rcpp::clone<Rcpp::NumericVector>(curr_parak);
  			    //List2["bs"] =  Rcpp::clone<Rcpp::NumericVector>(Beta);
  			    List2["Psi"] = psi_out;
  			    List2["Phi"] = phi_out;
  			    List2["Lik"] = LLik_out;
  			    List2["Fail_Stat"] = fail_out;
  			    List2["Iter"] = iter_out;
  			    List2["Duration"] =  final_duration;
  			    List2["KEG_Dur"] =  KEG_duration;
  			    List2["NB_Dur"] =  NB_duration;
  			    List2["Psi_Dur"] = Psi_duration;
  			    
  			    if(Rcpp::as<int>(List2["Fail_Stat"])==100){
  			      throw std::invalid_argument("Non-finite Optimization in Eta Hypothesis");
  			    }
  			    
  			    ex1[5] = 0;
  			    para0[0] = exp(curr_parak[0]);
  			    para0[1] = 1.0;
  			    para0[2] = exp(curr_parak[1]);
  			    TReC_sfit(para0.begin(),ex1.begin(),ex2.begin(),Xt.begin(),M,maxIter2,rhos_ASall.begin(),maxIter,eps,convg,lik_convg,
                      lower_psi.begin(),upper_psi.begin(),nbd_psi.begin(),curr_parap.begin(),lower_keg.begin(),upper_keg.begin(),
                      nbd_keg.begin(),curr_parak1a.begin(),init_parak.begin(),Beta1a.begin(),Beta_old.begin(),offset.begin(),Yint.begin(),
                      wa,iwa,g1,x1,&phi_out,&psi_out,&LLik_out,&fail_out,&iter_out,&final_duration,
                      &KEG_duration,&NB_duration,&Psi_duration,gv.begin(),tv.begin(),Xv.begin(),cv.begin(),Bv, mask_v.begin());
  			    
  			    List1a["Param"] =  curr_parak1a;
  			    List1a["bs"] =  Beta1a;
  			    List1a["Psi"] = psi_out;
  			    List1a["Phi"] = phi_out;
  			    List1a["Lik"] = LLik_out;
  			      List1a["Fail_Stat"] = fail_out;
  			    List1a["Iter"] = iter_out;
  			    List1a["Duration"] =  final_duration;
  			    List1a["KEG_Dur"] =  KEG_duration;
  			    List1a["NB_Dur"] =  NB_duration;
  			    List1a["Psi_Dur"] = Psi_duration;
  			    
  			    //Rprintf("Past Eta \n");
  			    
  			    ex1[5] = 2;
  			    para0[0] = 1.0;
  			    para0[1] = 1.0;
  			    TReC_sfit(para0.begin(),ex1.begin(),ex2.begin(),Xt.begin(),M,maxIter2,rhos_ASall.begin(),maxIter,eps,convg,lik_convg,
                      lower_psi.begin(),upper_psi.begin(),nbd_psi.begin(),curr_parap.begin(),lower_keg.begin(),upper_keg.begin(),
                      nbd_keg.begin(),curr_parak.begin(),init_parak.begin(),Beta.begin(),Beta_old.begin(),offset.begin(),Yint.begin(),
                      wa,iwa,g1,x1,&phi_out,&psi_out,&LLik_out,&fail_out,&iter_out,&final_duration,
                      &KEG_duration,&NB_duration,&Psi_duration,gv.begin(),tv.begin(),Xv.begin(),cv.begin(),Bv, mask_v.begin());
  			    
  			    //List3["Param"] =  Rcpp::clone<Rcpp::NumericVector>(curr_parak);
  			    //List3["bs"] =  Rcpp::clone<Rcpp::NumericVector>(Beta);
  			    List3["Psi"] = psi_out;
  			    List3["Phi"] = phi_out;
  			    List3["Lik"] = LLik_out;
  			    List3["Fail_Stat"] = fail_out;
  			    List3["Iter"] = iter_out;
  			    List3["Duration"] =  final_duration;
  			    List3["KEG_Dur"] =  KEG_duration;
  			    List3["NB_Dur"] =  NB_duration;
  			    List3["Psi_Dur"] = Psi_duration;
  			    
  			    if(Rcpp::as<int>(List3["Fail_Stat"])==100){
  			      throw std::invalid_argument("Non-finite Optimization in Gamma");
  			    }
  			    
  			    ex1[5] = 0;
  			    para0[0] = exp(curr_parak[0]);
  			    para0[1] = exp(curr_parak[1]);
  			    para0[2] = 1.0;
  			    TReC_sfit(para0.begin(),ex1.begin(),ex2.begin(),Xt.begin(),M,maxIter2,rhos_ASall.begin(),maxIter,eps,convg,lik_convg,
                      lower_psi.begin(),upper_psi.begin(),nbd_psi.begin(),curr_parap.begin(),lower_keg.begin(),upper_keg.begin(),
                      nbd_keg.begin(),curr_parak1b.begin(),init_parak.begin(),Beta1b.begin(),Beta_old.begin(),offset.begin(),Yint.begin(),
                      wa,iwa,g1,x1,&phi_out,&psi_out,&LLik_out,&fail_out,&iter_out,&final_duration,
                      &KEG_duration,&NB_duration,&Psi_duration,gv.begin(),tv.begin(),Xv.begin(),cv.begin(),Bv, mask_v.begin());
  			    List1b["Param"] =  curr_parak1b;
  			    List1b["bs"] =  Beta1b;
  			    List1b["Psi"] = psi_out;
  			    List1b["Phi"] = phi_out;
  			    List1b["Lik"] = LLik_out;
  			    List1b["Fail_Stat"] = fail_out;
  			    List1b["Iter"] = iter_out;
  			    List1b["Duration"] =  final_duration;
  			    List1b["KEG_Dur"] =  KEG_duration;
  			    List1b["NB_Dur"] =  NB_duration;
  			    List1b["Psi_Dur"] = Psi_duration;
  			    
  			    if((Rcpp::as<int>(List1a["Fail_Stat"])==100)&&(Rcpp::as<int>(List1b["Fail_Stat"])==100)){
  			      throw std::invalid_argument("Non-finite Optimization in Full Model");
  			    } else if((Rcpp::as<int>(List1a["Fail_Stat"])==100)&&(Rcpp::as<int>(List1b["Fail_Stat"])!=100)){
  			      List1 = List1b;
  			    } else if((Rcpp::as<int>(List1a["Fail_Stat"])!=100)&&(Rcpp::as<int>(List1b["Fail_Stat"])==100)){
  			      List1 = List1a;
  			    } else {
  			      if((Rcpp::as<double>(List1a["Lik"])>=Rcpp::as<double>(List1b["Lik"]))){
  			        List1 = List1a;
  			      } else {
  			        List1 = List1b;
  			      }
  			    }
  			    
  			    KEG = Rcpp::as<Rcpp::NumericVector>(List1["Param"]);
  			    
  			    // Conduct the Tests
  			    if(Rcpp::as<int>(List1["Fail_Stat"])!=1&&Rcpp::as<int>(List2["Fail_Stat"])!=1){
  			      lrt_eta = -2*(Rcpp::as<double>(List2["Lik"])-Rcpp::as<double>(List1["Lik"]));
  			      p_eta   = R::pchisq(lrt_eta,1,0,0);
  			    } else {
  			      p_eta   = -99.00;
  			      lrt_eta = -1.0e8;
  			    }
  			    
  			    if(Rcpp::as<int>(List1["Fail_Stat"])!=1&&Rcpp::as<int>(List3["Fail_Stat"])!=1){
  			      lrt_gamma = -2*(Rcpp::as<double>(List3["Lik"])-Rcpp::as<double>(List1["Lik"]));
  			      p_gamma   = R::pchisq(lrt_gamma,1,0,0);
  			    } else {
  			      p_gamma   = -99.00;
  			      lrt_gamma = -1.0e8;
  			    }
  			    
  			    //Rprintf("To Timing Calc\n");
  			    
  			    // Compile results of timers
  			    double Tot_dur = Rcpp::as<double>(List1a["Duration"])+Rcpp::as<double>(List1b["Duration"])+Rcpp::as<double>(List2["Duration"])+Rcpp::as<double>(List3["Duration"]);
  			    double KEG_dur = Rcpp::as<double>(List1a["KEG_Dur"])+Rcpp::as<double>(List1b["KEG_Dur"])+Rcpp::as<double>(List2["KEG_Dur"])+Rcpp::as<double>(List3["KEG_Dur"]);
  			    double NB_dur = Rcpp::as<double>(List1a["NB_Dur"])+Rcpp::as<double>(List1b["NB_Dur"])+Rcpp::as<double>(List2["NB_Dur"])+Rcpp::as<double>(List3["NB_Dur"]);
  			    double Psi_dur = Rcpp::as<double>(List1a["Psi_Dur"])+Rcpp::as<double>(List1b["Psi_Dur"])+Rcpp::as<double>(List2["Psi_Dur"])+Rcpp::as<double>(List3["Psi_Dur"]);
  			    
  			    //Rprintf("Past Timing Calc\n");
  			    
  			    // Store the Output
  			    Betatmp = Rcpp::as<Rcpp::NumericVector>(List1["bs"]);
  			    fprintf(fcomp,"%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.18f\t%.10f\t%.18f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\t%d\t%d\t%d\t%d",j,(i+gene_add),Tot_dur,KEG_dur,NB_dur,Psi_dur,lrt_eta,p_eta,lrt_gamma,p_gamma,Rcpp::as<double>(List1["Lik"]),Rcpp::as<double>(List2["Lik"]),Rcpp::as<double>(List3["Lik"]),exp(KEG[0]),exp(KEG[1]),exp(KEG[2]),Rcpp::as<double>(List1["Phi"]),Rcpp::as<double>(List1["Psi"]),Rcpp::as<int>(List1["Fail_Stat"]),Rcpp::as<int>(List1a["Fail_Stat"]),Rcpp::as<int>(List1b["Fail_Stat"]),Rcpp::as<int>(List2["Fail_Stat"]),Rcpp::as<int>(List3["Fail_Stat"]));
  			    for(k=0;k<M;k++){
  			      fprintf(fcomp,"\t%.10f",Betatmp[k]);
  			    }
  			    fprintf(fcomp,"\t%.10f\t%.10f\t%d\n",-1e8,1.0,5);
  			    
  			    //Rprintf("Past Storing Output\n");
  			  } catch(std::invalid_argument& e){
  			    Rprintf("Error Caught!\n");
  			    fprintf(fcomp,"%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\t%d\t%d\t%d\t%d",j,(i+gene_add),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,4,4,4,4,4);
  			    for(k=0;k<M;k++){
  			      fprintf(fcomp,"\t%.10f",0.0);
  			    }
  			    fprintf(fcomp,"\t%.10f\t%.10f\t%d\n",-1e8,1.0,5);
  			  }
  			  
  			  // Parameter reset
  			  para0[0] = 1.0;
  			  para0[1] = 1.0;
  			  para0[2] = 1.0;
  			  
  			  //Rprintf("Past Parameter Reset\n");
  			  
  			} 
  			/* End a single Gene-SNP Pair test. Reset the the value of hte
  			 * Score Test P-Value */
  			CT_pval = 1.0;
			}
		}
	}
	//Rprintf("To Memory Free Portion!\n");
	fclose(fcomp);
	UNPROTECT(1);
	return Rcpp::wrap(out);
}

// Failure Status Codes:
// 0   - Converged
// 1   - Optimization routine failure (likelihood decrease)
// 2   - Iteration Limit Reached
// 3   - KEG optimization failure causing a cessation in parameter updates
// 4   - Results from a call of 100 (see below) and prints out no information for test
// 100 - Optimization routine failure in LBFGS (infinite likelihood) 

// CT Failure Code:
// 0 - Observed Score Works
// 1 - Used CT_Expected and it Worked
// 2 - Invertibility issues with Hessian
// 5 - Did not Conduct
