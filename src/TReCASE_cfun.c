#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "Defn1.h"
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <time.h>
#include "TReC.h"
#include "ASE.h"

//#define DEBUG 1

/*****************************************************************
 * (0) Specifying the lbfgsb1                                    *
 *****************************************************************/
void lbfgsb1(int n, int m, double *x, double *l, double *u, int *nbd,
             double *Fmin, optimfn1 fminfn1, optimgr1 fmingr1, int *fail,
             void *ex, double factr, double pgtol,
             int *fncount, int *grcount, int maxit, char *msg,
             int trace, int nREPORT, double *wa, int *iwa, double *g, SEXP x1);

/*****************************************************************
 * (1) Combined TReC+ASE likelihoods                             *
 *****************************************************************/
double neg_logLik_TS(int n, double* para, void* ex, SEXP x1){
	/* Main Function Body */
	#ifdef DEBUG
		Rprintf("Current Parameters:\n log(Kappa): %5.5f\n log(Eta): %5.5f\n log(Gamma): %5.5f\n",
				 para[0],para[1],para[2]);
	#endif
	
	double logL = neg_logLTReC_TS(n,para,ex,x1)+neg_logLASE_TS(n,para,ex,x1);
	
	#ifdef DEBUG
		Rprintf("Resulting LogL = %f\n",logL);
	#endif
	
	/* Return Results */
	return logL;
}

/*****************************************************************
 * (2) Combined TReC+ASE gradients                               *
 *****************************************************************/
 // Acts as a wrapper for the function to call from lbfgsb1. Creates
 // a memory hold for the hessian that it deletes while doing nothing
 // with it.
 
void neg_gradLik_TSnh(int n, double* para, double* gr, void* ex, SEXP x1){
	/* Declare and Initialize */
	// hh is deprecated. Do not use it.
	double hh;
	
	/* Function Body */
	// To remove calls to R's memory management system, I create three variables
	// which C will destroy at the conclusion of the function.
	double g1 = 0.0, g2=0.0, g3 = 0.0;
	
	// Trec gradient
	neg_gradTReC_TS(n,para,gr,&hh,ex,x1);
	  g1 = gr[0];
	  g2 = gr[1];
	  if(n==3){
	    g3 = gr[2];
	  }

	
	// ASE gradient
	neg_gradASE_TS(n,para,gr,&hh,ex,x1);
	
	// Combine
  gr[0]+=g1;
  gr[1]+=g2;
  if(n==3){
    gr[2]+=g3;
  }
}

void neg_gradTReC_TSnh(int n, double* para, double* gr, void* ex, SEXP x1){
  /* Memory Allocation using R */
  double hh = 0.0;
  
  /* Function Body */
  // hh is deprecated. Do not use it!
  neg_gradTReC_TS(n,para,gr,&hh,ex,x1);
}


/*****************************************************************
 * (3) Combined TReC+ASE gradients                               *
 *****************************************************************/
/* Hessian Not Needed! Hence deprecation of hh.
void neg_gradLik_TSh(int n, double* para, double* gr, double* Hessian, void* ex){
	// Declare and Initialize 
	int i=0, j=0;
	
	double* trec_gr = new double[ n ];
	double* trec_h  = new double[ (n*n) ];
	double* ase_gr  = new double[ n ];
	double* ase_h   = new double[ (n*n) ];
	
	// Function Body 
	neg_gradTReC_TS(n,para,trec_gr,trec_h,ex);
	neg_gradASE_TS(n,para,ase_gr,ase_h,ex);
	
	// Allocate Gradient and Initial Hessian
	for(i=0;i<n;i++){
		gr[i] = trec_gr[i]+ase_gr[i];
	}
	
	for(j=0;j<(n*n);j++){
		Hessian[i] = trec_h[j]+ase_h[j];
	}
	
	// Deallocate Heap Memory 
	delete[] trec_gr;
	delete[] trec_h;
	delete[] ase_gr;
	delete[] ase_h;
}
*/

/*****************************************************************
 * (4) vmmin_dj                                                  *
 *****************************************************************/
/*  BFGS variable-metric method, based on Pascal code
in J.C. Nash, `Compact Numerical Methods for Computers', 2nd edition,
converted by p2c then re-crafted by B.D. Ripley */

#define stepredn	0.2
#define acctol		0.0001
#define reltest		10.0

void
  vmmin_dj(int n0, double *b, double *Fmin, optimfn1 fminfn, optimgr1 fmingr,
           int maxit, int trace, int *mask,
           double abstol, double reltol, int nREPORT, void *ex,
           int *fncount, int *grcount, int *fail,
           double* g, double* t, double* X, double* c, double** B, SEXP x1)
  {
    Rboolean accpoint, enough;
    //double *g, *t, *X, *c, **B;
    int   count, funcount, gradcount;
    double f, gradproj;
    int   i, j, ilast, iter = 0;
    double s, steplength;
    double D1, D2;
    int   n, *l;
    
    if (maxit <= 0) {
      *fail = 0;
      *Fmin = fminfn(n0, b, ex,x1);
      *fncount = *grcount = 0;
      return;
    }
    
    if (nREPORT <= 0)
      error(_("REPORT must be > 0 (method = \"BFGS\")"));
    l = (int *) R_alloc(n0, sizeof(int));
    n = 0;
    for (i = 0; i < n0; i++) if (mask[i]) l[n++] = i;
    //g = vect(n0);
    //t = vect(n);
    //X = vect(n);
    //c = vect(n);
    //B = Lmatrix(n);
    f = fminfn(n0, b, ex,x1);
    if (!R_FINITE(f)){
      *fail = 100;
      return;
    }
    if (trace) Rprintf("initial  value %f \n", f);
    *Fmin = f;
    funcount = gradcount = 1;
    fmingr(n0, b, g, ex,x1);
    iter++;
    ilast = gradcount;
    
    do {
      if (ilast == gradcount) {
        for (i = 0; i < n; i++) {
          for (j = 0; j < i; j++) B[i][j] = 0.0;
          B[i][i] = 1.0;
        }
      }
      for (i = 0; i < n; i++) {
        X[i] = b[l[i]];
        c[i] = g[l[i]];
      }
      gradproj = 0.0;
      for (i = 0; i < n; i++) {
        s = 0.0;
        for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
        for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
        t[i] = s;
        gradproj += s * g[l[i]];
      }
      
      if (gradproj < 0.0) {	/* search direction is downhill */
steplength = 1.0;
        accpoint = FALSE;
        do {
          count = 0;
          for (i = 0; i < n; i++) {
            b[l[i]] = X[i] + steplength * t[i];
            if (reltest + X[i] == reltest + b[l[i]]) /* no change */
count++;
          }
          if (count < n) {
            f = fminfn(n0, b, ex,x1);
            funcount++;
            accpoint = R_FINITE(f) &&
              (f <= *Fmin + gradproj * steplength * acctol);
            if (!accpoint) {
              steplength *= stepredn;
            }
          }
        } while (!(count == n || accpoint));
        enough = (f > abstol) &&
          fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);
        /* stop if value if small or if relative change is low */
        if (!enough) {
          count = n;
          *Fmin = f;
        }
        if (count < n) {/* making progress */
        *Fmin = f;
          fmingr(n0, b, g, ex,x1);
          gradcount++;
          iter++;
          D1 = 0.0;
          for (i = 0; i < n; i++) {
            t[i] = steplength * t[i];
            c[i] = g[l[i]] - c[i];
            D1 += t[i] * c[i];
          }
          if (D1 > 0) {
            D2 = 0.0;
            for (i = 0; i < n; i++) {
              s = 0.0;
              for (j = 0; j <= i; j++)
                s += B[i][j] * c[j];
              for (j = i + 1; j < n; j++)
                s += B[j][i] * c[j];
              X[i] = s;
              D2 += s * c[i];
            }
            D2 = 1.0 + D2 / D1;
            for (i = 0; i < n; i++) {
              for (j = 0; j <= i; j++)
                B[i][j] += (D2 * t[i] * t[j]
                              - X[i] * t[j] - t[i] * X[j]) / D1;
            }
          } else {	/* D1 < 0 */
        ilast = gradcount;
          }
        } else {	/* no progress */
        if (ilast < gradcount) {
          count = 0;
          ilast = gradcount;
        }
        }
      } else {		/* uphill search */
        count = 0;
        if (ilast == gradcount) count = n;
        else ilast = gradcount;
        /* Resets unless has just been reset */
      }
      if (trace && (iter % nREPORT == 0))
        Rprintf("iter%4d value %f\n", iter, f);
      if (iter >= maxit) break;
      if (gradcount - ilast > 2 * n)
        ilast = gradcount;	/* periodic restart */
    } while (count != n || ilast != gradcount);
    if (trace) {
      Rprintf("final  value %f \n", *Fmin);
      if (iter < maxit) Rprintf("converged\n");
      else Rprintf("stopped after %i iterations\n", iter);
    }
    *fail = (iter < maxit) ? 0 : 1;
    *fncount = funcount;
    *grcount = gradcount;
  }



/*****************************************************************
 * (5) lbfgsb_wrappers                                            *
 *****************************************************************/
void KEG_optroutine(int npara,int lmm,double* initPara,double* lower,
				double* upper, int* nbd, double* Fmin, int* fail, double* ex,
				double factr, double pgtol, int* fncount, int* grcount, int maxit1,
				char* msg, int nREPORT, double* wa, int* iwa, double* g1, SEXP x1){
	// Simple wrapper for lbfgsb
	lbfgsb1(npara,lmm,initPara,lower,upper,nbd,Fmin,neg_logLik_TS,
		        neg_gradLik_TSnh, fail, (void *) ex, factr, pgtol,
		        fncount, grcount, maxit1, msg, 0, nREPORT,wa,iwa,g1,x1);
}


void KEG_boptroutine(int npara,double* initPara,double* Fmin, int maxit1,
                     int* mask, double abstol, double reltol, int nReport,
                     double* ex, int* fncount, int*grcount, int* fail,
                     double* gv, double* tv, double* Xv, double* cv, double** Bv,
                     SEXP x1){
  vmmin_dj(npara,initPara,Fmin,neg_logLik_TS,neg_gradLik_TSnh,maxit1,
        0,mask,abstol,reltol,nReport,(void *)ex, fncount,grcount,fail,
        gv,tv,Xv,cv,Bv,x1);
                     }

void psi_optroutine(int npara,int lmm,double* initPara,double* lower,
				double* upper, int* nbd, double* Fmin, int* fail, double* ex,
				double factr, double pgtol, int* fncount, int* grcount, int maxit1,
				char* msg, int nREPORT, double* wa, int* iwa, double* g1, SEXP x1){
	// Simple wrapper for lbfgsb
	lbfgsb1(npara,lmm,initPara,lower,upper,nbd,Fmin,logLASE_psi,
		        gradASE_psi, fail, (void *) ex, factr, pgtol,
		        fncount, grcount, maxit1, msg, 0, nREPORT,wa,iwa,g1,x1);
}


void psi_boptroutine(int npara,double* initPara,double* Fmin, int maxit1,
                     int* mask, double abstol, double reltol, int nReport,
                     double* ex, int* fncount, int*grcount, int* fail,
                     double* gv, double* tv, double* Xv, double* cv, double** Bv,
                     SEXP x1){
  vmmin_dj(npara,initPara,Fmin,logLASE_psi,gradASE_psi,maxit1,
           0,mask,abstol,reltol,nReport,(void *)ex, fncount,grcount,fail,gv,tv,Xv,cv,Bv,x1);
                     }

void KEG_optroutine2(int npara,int lmm,double* initPara,double* lower,
                    double* upper, int* nbd, double* Fmin, int* fail, double* ex,
                    double factr, double pgtol, int* fncount, int* grcount, int maxit1,
                    char* msg, int nREPORT, double* wa, int* iwa, double* g1, SEXP x1){
  // Simple wrapper for lbfgsb
  lbfgsb1(npara,lmm,initPara,lower,upper,nbd,Fmin,neg_logLTReC_TS,
            neg_gradTReC_TSnh, fail, (void *) ex, factr, pgtol,
            fncount, grcount, maxit1, msg, 0, nREPORT,wa,iwa,g1,x1);
}


void KEG_boptroutine2(int npara,double* initPara,double* Fmin, int maxit1,
                     int* mask, double abstol, double reltol, int nReport,
                     double* ex, int* fncount, int*grcount, int* fail,
                     double* gv, double* tv, double* Xv, double* cv, double** Bv,
                     SEXP x1){
  vmmin_dj(npara,initPara,Fmin,neg_logLTReC_TS,neg_gradTReC_TSnh,maxit1,
           0,mask,abstol,reltol,nReport,(void *)ex, fncount,grcount,fail,
           gv,tv,Xv,cv,Bv,x1);
                     }
