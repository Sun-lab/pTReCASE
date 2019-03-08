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
//#define DEBUG_err 1

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
                 double* KEG_duration, double* NB_duration, double* Psi_duration,
                 double* gv, double* tv, double* Xv, double* cv, double** Bv, int* mask_v){
  // Starting the Timer 
  clock_t tstart = clock();
  
  // Setting up Data-specific Parameters:
  int n_subj = 0, H0 = 0;
  
  n_subj  = ceil(ex1[0]-0.5);
  
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
    
    KEG_optroutine2(npara_keg,lmm,curr_parak,lower_keg,upper_keg,
                   nbd_keg,&Fmin,&fail,ex1,factr,pgtol,
                   &fncount,&grcount,maxit1,msg,nREPORT,wa,iwa,g1,x1);
    
    if(fail>1){
      // KEG are the horses of this routine, we should attempt to improve
      // fit in case of failure. Move me to a position where I won't fail!
      
      for(h=0;h<npara_keg;h++){
        curr_parak[h] = init_parak[h];
      }
      
      fail = 0;
      KEG_boptroutine2(npara_keg, curr_parak, &Fmin, maxit1,mask_v,abstol,reltol,
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
        *psi_out   = -1.0;
        *phi_out   = phi_curr;
        *LLik_out   = logLik1;
        *fail_out = fail_stat;
        *iter_out      = g;
        *final_duration  = 0;
        *KEG_duration   = 0;
        *NB_duration    = 0;
        *Psi_duration  = 0.0;
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
    
    if(parDiff<pdTmp){
      parDiff = pdTmp;
    }
    
    //Rprintf("Phi parDiff: %.7f\n",parDiff);
    
    ex1[7] = phi_curr;
    
    //	Rprintf("Current Phi (%d): %.7f\n",g,phi_curr);
    
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
  *psi_out   = -1.00;
  *phi_out  = phi_curr;
  *LLik_out  = logLik1;
  *fail_out = fail_stat;
  *iter_out = g;
  *final_duration  = fin_dur;
  *KEG_duration   = KEG_dur;
  *NB_duration   = NB_dur;
  *Psi_duration = 0.00;
  
  //Rprintf("Made it to the return statement\n");
  
  return(fail_stat);
}  

