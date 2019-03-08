#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <time.h>

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
 
 /****************************************************************
  * (1) loglikNB_pmy
  ****************************************************************/
  
double loglikNB_pmy (double phi, int n, double* mu, int* y){
	double logL = 0.0, vphi, lik0 = 0.0, mui=0.0;
	int i=0, yi=0;
	
	vphi = 1/phi;
	lik0 = vphi*log(vphi)-lgamma(vphi);
	
	for(i=0;i<n;i++){
		yi = y[i];
		mui = mu[i];
		if(yi==0){
			logL = logL + vphi*log(vphi) - vphi*log(vphi+mui);
		} else {
			logL = logL + lgamma(yi + vphi) - lgamma(yi + 1.0) + yi*log(mui);
			logL = logL - (vphi + yi)*log(vphi + mui) + lik0;
		}
	}
	
	return(logL);
}

/****************************************************************
 * (2) neg.logLTReC.TS
 ****************************************************************/
  
 double neg_logLTReC_TS(int n, double* para, void* ex, SEXP x1){
	 /* Initialize the variables */
	 double kappa = 0.0, eta = 0.0, gamma = 0.0, phi, vphi, delta, cs, zeta, logL = 0.0,
	        lik0 = 0.0, mui=0.0;
	 double* exPara;
	 int i=0,n_subj=0,H0,idxc=0,idxr=0,idxn=0,
		 n0 = 0, n1 = 0, yi = 0;
	 
	 /* Deal with Extra Parameters */
	 exPara = (double *) ex;
	 n_subj = ceil(exPara[0]-0.5);
	 n0     = ceil(exPara[2]-0.5);
	 n1     = ceil(exPara[3]-0.5);
	 H0     = ceil(exPara[5]-0.5);
	 phi    = exPara[7];
	 
	 kappa = exp(para[0]);
	 
	 if(H0 == 0){
		eta = exp(para[1]);
		gamma = exp(para[2]);
	} else if(H0 == 1){
		eta = 1.0;
		gamma = exp(para[1]);
	} else if(H0 == 2){
		eta = exp(para[1]);
		gamma = 1.0;
	}

	vphi = 1.0/phi;
	lik0 = vphi*log(vphi)-lgamma(vphi);
	
	 for(i=0;i<n0;i++){
		 idxc = 9+i;
		 idxr = idxc+n_subj;
		 idxn = idxr+n_subj;
		 
		 yi = ceil(exPara[idxc]-0.5);
		 
		 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
		 cs    = (double)((1.0 - exPara[idxr])/delta);
		 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
		 
		 mui = exPara[idxn]*delta;
		 
		 if(yi==0){
		   logL = logL + vphi*log(vphi) - vphi*log(vphi+mui);
		 } else {
		   logL = logL + lgamma(yi + vphi) - lgamma(yi + 1.0) + yi*log(mui);
		   logL = logL - (vphi + yi)*log(vphi + mui) + lik0;
		 }
	 }
	 
	 for(i=n0;i<(n0+n1);i++){
		 idxc = 9+i;
		 idxr = idxc+n_subj;
		 idxn = idxr+n_subj;
		 
		 yi = ceil(exPara[idxc]-0.5);
		 
		 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
		 cs    = (double)((1.0 - exPara[idxr])/delta);
		 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
		 
		 mui = exPara[idxn]*delta;
		 mui *= (double)((1.0 + zeta)/2.0);
		 
		 if(yi==0){
		   logL = logL + vphi*log(vphi) - vphi*log(vphi+mui);
		 } else {
		   logL = logL + lgamma(yi + vphi) - lgamma(yi + 1.0) + yi*log(mui);
		   logL = logL - (vphi + yi)*log(vphi + mui) + lik0;
		 }
	 }
	 
	 for(i=(n0+n1);i<n_subj;i++){
		 idxc = 9+i;
		 idxr = idxc+n_subj;
		 idxn = idxr+n_subj;
		 
		 yi = ceil(exPara[idxc]-0.5);
		 
		 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
		 cs    = (double)((1.0 - exPara[idxr])/delta);
		 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
		 
		 mui = exPara[idxn]*delta;
		 mui *= (double) zeta;
		 
		 if(yi==0){
		   logL = logL + vphi*log(vphi) - vphi*log(vphi+mui);
		 } else {
		   logL = logL + lgamma(yi + vphi) - lgamma(yi + 1.0) + yi*log(mui);
		   logL = logL - (vphi + yi)*log(vphi + mui) + lik0;
		 }
	 }
	 
	 return(-logL);
 }

/*****************************************************************
 * (3) neg.gradTReC.TS
 *****************************************************************/
 
 /* Code is long and somewhat unwieldy. This is due to the "unpacking"
	of the if-else statements to ensure that loops are not constantly
	checking conditions across all subjects.*/
 
 void neg_gradTReC_TS(int n, double* para, double* gr, 
					  double* Hessian, void* ex, SEXP x1){
	 /* Initialize the variables */
	 double kappa = 0.0, eta = 0.0, gamma = 0.0, phi, delta, cs, zeta, nu1,
			sum1= 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0,
			sum5 = 0.0, sum6 = 0.0, sum7 = 0.0,dl_dmu,
			dl_dkappa = 0.0, dl_dlogeta = 0.0, dl_dloggamma = 0.0;
	 double dmu_dkappa = 0.0, dmu_deta = 0.0, dmu_dgamma_p = 0.0,
			d2l_dmu2 = 0.0;
	 double* exPara;
	 int i=0,n_subj=0,H0,Hess,idxc=0,idxr=0,idxn=0,yi,
		 n0 = 0, n1 = 0;
	 
	 /* Deal with Extra Parameters */
	 exPara = (double *) ex;
	 n_subj = ceil(exPara[0]-0.5);
	 n0     = ceil(exPara[2]-0.5);
	 n1     = ceil(exPara[3]-0.5);
	 H0     = ceil(exPara[5]-0.5);
	 Hess   = ceil(exPara[6]-0.5);
	 phi    = exPara[7];
	 
	 kappa = exp(para[0]);
	 
	 if(Hess == 1){
		 if(H0 == 0){
			/*****************************************************
			 * Parameters 
			 *****************************************************/
			 
			 eta = exp(para[1]);
			 gamma = exp(para[2]);
			 
			/****************************************************
			 * Gradient Computation
			 ****************************************************/
			 
			 for(i=0;i<n0;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 nu1 = exPara[idxn]*delta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 dmu_dkappa = exPara[idxn]*exPara[idxr];
				 dl_dkappa += dl_dmu*dmu_dkappa;

				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/
				  
				 dmu_deta = 0.0;
				 dl_dlogeta += 0.0;
				 
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = 0.0;
				 dl_dloggamma += 0.0;
				 
				/*********************************************************
				 * Hessian Computation   
				 *********************************************************/
				 d2l_dmu2 = (double)(-yi)/(nu1*nu1); 
				 d2l_dmu2 += phi*(1.0+phi*yi)/((1.0+phi*nu1)*(1.0+phi*nu1));
			
				 sum1 += d2l_dmu2*dmu_dkappa*dmu_dkappa;
				 sum2 += d2l_dmu2*dmu_deta*dmu_deta;
				 sum3 += d2l_dmu2*dmu_deta*dmu_dkappa;
				 sum4 += d2l_dmu2*dmu_dgamma_p*dmu_dgamma_p;
				 sum5 += d2l_dmu2*dmu_dgamma_p*dmu_dkappa;
				 sum6 += 0.0;
				 sum7 += d2l_dmu2*dmu_dgamma_p*dmu_deta;
			 }
			 
			 for(i=n0;i<(n0+n1);i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double)((1.0+zeta)/2.0);
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = (0.5*exPara[idxn]*exPara[idxr]*(1.0+gamma));
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				 /*********************************************************
				  * dl.dlogeta  
				  *********************************************************/
				  
				 dmu_deta = 0.5*exPara[idxn]*(1.0 - exPara[idxr]);
				 dl_dlogeta += dl_dmu*dmu_deta; 
				 
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = 0.5*exPara[idxn]*exPara[idxr];
				 dl_dloggamma += dl_dmu*dmu_dgamma_p;
		
				/*********************************************************
				 * Hessian Computation   
				 *********************************************************/
				 d2l_dmu2 = (double)(-yi)/(nu1*nu1); 
				 d2l_dmu2 += phi*(1.0+phi*yi)/((1.0+phi*nu1)*(1.0+phi*nu1));
			
				 sum1 += d2l_dmu2*dmu_dkappa*dmu_dkappa;
				 sum2 += d2l_dmu2*dmu_deta*dmu_deta;
				 sum3 += d2l_dmu2*dmu_deta*dmu_dkappa;
				 sum4 += d2l_dmu2*dmu_dgamma_p*dmu_dgamma_p;
				 sum5 += d2l_dmu2*dmu_dgamma_p*dmu_dkappa;
				 sum6 += dl_dmu*0.5*exPara[idxn]*exPara[idxr];
				 sum7 += d2l_dmu2*dmu_dgamma_p*dmu_deta;
			 }
			 
			 for(i=(n0+n1);i<n_subj;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double) zeta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = exPara[idxn]*exPara[idxr]*gamma;
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				 /*********************************************************
				  * dl.dlogeta  
				  *********************************************************/
				  
				 dmu_deta = exPara[idxn]*(1.0 - exPara[idxr]);
				 dl_dlogeta += dl_dmu*dmu_deta;
				 
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = exPara[idxn]*exPara[idxr];
				 dl_dloggamma += dl_dmu*dmu_dgamma_p;
				 
				/*********************************************************
				 * Hessian Computation   
				 *********************************************************/
				 d2l_dmu2 = (double)(-yi)/(nu1*nu1); 
				 d2l_dmu2 += phi*(1.0+phi*yi)/((1.0+phi*nu1)*(1.0+phi*nu1));
			
				 sum1 += d2l_dmu2*dmu_dkappa*dmu_dkappa;
				 sum2 += d2l_dmu2*dmu_deta*dmu_deta;
				 sum3 += d2l_dmu2*dmu_deta*dmu_dkappa;
				 sum4 += d2l_dmu2*dmu_dgamma_p*dmu_dgamma_p;
				 sum5 += d2l_dmu2*dmu_dgamma_p*dmu_dkappa;
				 sum6 += dl_dmu*exPara[idxn]*exPara[idxr];
				 sum7 += d2l_dmu2*dmu_dgamma_p*dmu_deta;
			 } 
		 } else if(H0 == 1){
			/*****************************************************
			 * Parameters 
			 *****************************************************/
			 
			 eta = 1.0;
			 gamma = exp(para[1]);
			 
			/****************************************************
			 * Gradient Computation
			 ****************************************************/
			 
			 for(i=0;i<n0;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 nu1 = exPara[idxn]*delta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 dmu_dkappa = exPara[idxn]*exPara[idxr];
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				 /*********************************************************
				  * dl.dlogeta  
				  *********************************************************/
				  
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = 0.0;
				 dl_dloggamma += 0.0;

				/*********************************************************
				 * Hessian Computation   
				 *********************************************************/
				 d2l_dmu2 = (double)(-yi)/(nu1*nu1); 
				 d2l_dmu2 += phi*(1.0+phi*yi)/((1.0+phi*nu1)*(1.0+phi*nu1));
			
				 sum1 += d2l_dmu2*dmu_dkappa*dmu_dkappa;
				 sum4 += d2l_dmu2*dmu_dgamma_p*dmu_dgamma_p;
				 sum5 += d2l_dmu2*dmu_dgamma_p*dmu_dkappa;	
				 sum6 += 0.0;
			 }
			 
			 for(i=n0;i<(n0+n1);i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double)((1.0+zeta)/2.0); 
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = (0.5*exPara[idxn]*exPara[idxr]*(1.0+gamma));
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				 /*********************************************************
				  * dl.dlogeta  
				  *********************************************************/
				  
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = 0.5*exPara[idxn]*exPara[idxr];
				 dl_dloggamma += dl_dmu*dmu_dgamma_p;

				/*********************************************************
				 * Hessian Computation   
				 *********************************************************/
				 d2l_dmu2 = (double)(-yi)/(nu1*nu1); 
				 d2l_dmu2 += phi*(1.0+phi*yi)/((1.0+phi*nu1)*(1.0+phi*nu1));
			
				 sum1 += d2l_dmu2*dmu_dkappa*dmu_dkappa;
				 sum4 += d2l_dmu2*dmu_dgamma_p*dmu_dgamma_p;
				 sum5 += d2l_dmu2*dmu_dgamma_p*dmu_dkappa;	
				 sum6 += dl_dmu*0.5*exPara[idxn]*exPara[idxr];
			 }
			 
			 for(i=(n0+n1);i<n_subj;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double) zeta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = exPara[idxn]*exPara[idxr]*gamma;
				 dl_dkappa += dl_dmu*dmu_dkappa;

				/*********************************************************
				  * dl.dlogeta  
				  *********************************************************/
				  
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = exPara[idxn]*exPara[idxr];
				 dl_dloggamma += dl_dmu*dmu_dgamma_p;
	
				/*********************************************************
				 * Hessian Computation   
				 *********************************************************/
				 d2l_dmu2 = (double)(-yi)/(nu1*nu1); 
				 d2l_dmu2 += phi*(1.0+phi*yi)/((1.0+phi*nu1)*(1.0+phi*nu1));
			
				 sum1 += d2l_dmu2*dmu_dkappa*dmu_dkappa;
				 sum4 += d2l_dmu2*dmu_dgamma_p*dmu_dgamma_p;
				 sum5 += d2l_dmu2*dmu_dgamma_p*dmu_dkappa;	
				 sum6 += dl_dmu*exPara[idxn]*exPara[idxr];
			 }

		 } else if(H0 == 2){
			/*****************************************************
			 * Parameters 
			 *****************************************************/
			 
			 eta = exp(para[1]);
			 gamma = 1.0;
			 
			/****************************************************
			 * Gradient Computation
			 ****************************************************/
			 
			 for(i=0;i<n0;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 nu1 = exPara[idxn]*delta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 dmu_dkappa = exPara[idxn]*exPara[idxr];
				 dl_dkappa += dl_dmu*dmu_dkappa;

				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/
				  
				 dmu_deta = 0.0;
				 dl_dlogeta += 0.0;
				 
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 
				/*********************************************************
				 * Hessian Computation   
				 *********************************************************/
				 d2l_dmu2 = (double)(-yi)/(nu1*nu1); 
				 d2l_dmu2 += phi*(1.0+phi*yi)/((1.0+phi*nu1)*(1.0+phi*nu1));
			
				 sum1 += d2l_dmu2*dmu_dkappa*dmu_dkappa;
				 sum2 += d2l_dmu2*dmu_deta*dmu_deta;
				 sum3 += d2l_dmu2*dmu_deta*dmu_dkappa;
			 }
			 
			 for(i=n0;i<(n0+n1);i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double)((1.0+zeta)/2.0); 
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = (0.5*exPara[idxn]*exPara[idxr]*(1.0+gamma));
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				 /*********************************************************
				  * dl.dlogeta  
				  *********************************************************/
				  
				 dmu_deta = 0.5*exPara[idxn]*(1.0 - exPara[idxr]);
				 dl_dlogeta += dl_dmu*dmu_deta; 
				 
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/

				/*********************************************************
				 * Hessian Computation   
				 *********************************************************/
				 d2l_dmu2 = (double)(-yi)/(nu1*nu1); 
				 d2l_dmu2 += phi*(1.0+phi*yi)/((1.0+phi*nu1)*(1.0+phi*nu1));
			
				 sum1 += d2l_dmu2*dmu_dkappa*dmu_dkappa;
				 sum2 += d2l_dmu2*dmu_deta*dmu_deta;
				 sum3 += d2l_dmu2*dmu_deta*dmu_dkappa;
			 }
			 
			 for(i=(n0+n1);i<n_subj;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double) zeta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = exPara[idxn]*exPara[idxr]*gamma;
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				 /*********************************************************
				  * dl.dlogeta  
				  *********************************************************/
				  
				 dmu_deta = exPara[idxn]*(1.0 - exPara[idxr]);
				 dl_dlogeta += dl_dmu*dmu_deta;
				 
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/

				/*********************************************************
				 * Hessian Computation   
				 *********************************************************/
				 d2l_dmu2 = (double)(-yi)/(nu1*nu1); 
				 d2l_dmu2 += phi*(1.0+phi*yi)/((1.0+phi*nu1)*(1.0+phi*nu1));
			
				 sum1 += d2l_dmu2*dmu_dkappa*dmu_dkappa;
				 sum2 += d2l_dmu2*dmu_deta*dmu_deta;
				 sum3 += d2l_dmu2*dmu_deta*dmu_dkappa;
			 }

		 }
		 
		/************************************************************
		 * FINISH GRADIENT and HESSIAN
		 ************************************************************/
		/* Gradient */
		 dl_dkappa    *= kappa;
		 dl_dlogeta   *= eta;
	     dl_dloggamma *= kappa*gamma;
		 
		 gr[0] = -dl_dkappa;
		 
		 if(H0 == 0){
			/* Gradient */
			 gr[1] = -dl_dlogeta;
			 gr[2] = -dl_dloggamma;
			 
			/* Hessian */
			sum4 *= kappa*kappa;
			sum5 *= kappa;
			sum7 *= kappa; 
		 
			Hessian[0] = -(dl_dkappa+kappa*kappa*sum1);
			Hessian[1] = -(eta*kappa*sum3);
			Hessian[2] = -(gamma*kappa*(sum5+sum6));
			Hessian[3] = Hessian[1];
			Hessian[4] = -(dl_dlogeta+eta*eta*sum2);
			Hessian[5] = -(eta*gamma*sum7);
			Hessian[6] = Hessian[2];
			Hessian[7] = Hessian[5];
			Hessian[8] = -(dl_dloggamma+gamma*gamma*sum4);
		} else if(H0 == 1){
		   /* Gradient */
		    gr[1] = -dl_dloggamma;
			
		   /* Hessian */
			sum4 *= kappa*kappa;
			sum5 *= kappa;
			sum7 *= kappa;
			
			Hessian[0] = -(dl_dkappa+kappa*kappa*sum1);
			Hessian[1] = -(gamma*kappa*(sum5+sum6));
			Hessian[2] = Hessian[1];
			Hessian[3] = -(dl_dloggamma+gamma*gamma*sum4);
		} else if(H0 == 2){
		   /* Gradient */
		    gr[1] = -dl_dlogeta;
			
		   /* Hessian */
			Hessian[0] = -(dl_dkappa+kappa*kappa*sum1);
			Hessian[1] = -(eta*kappa*sum3);
			Hessian[2] = Hessian[1];
			Hessian[3] = -(dl_dlogeta+eta*eta*sum2);
		}
	 } else {
		 if(H0 == 0){
			/*****************************************************
			 * Parameters 
			 *****************************************************/
			 
			 eta = exp(para[1]);
			 gamma = exp(para[2]); 
			 
			/****************************************************
			 * Gradient Computation
			 ****************************************************/
			 
			 for(i=0;i<n0;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 nu1 = exPara[idxn]*delta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 dmu_dkappa = exPara[idxn]*exPara[idxr];
				 dl_dkappa += dl_dmu*dmu_dkappa;

				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/
				 
				 dmu_deta = 0.0;
				 dl_dlogeta += 0.0;
				 
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = 0.0;
				 dl_dloggamma += 0.0;
			 }
			 
			 for(i=n0;i<(n0+n1);i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double)((1.0+zeta)/2.0); 
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = (0.5*exPara[idxn]*exPara[idxr]*(1.0+gamma));
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/
					
				 dmu_deta = 0.5*exPara[idxn]*(1.0-exPara[idxr]);
				 dl_dlogeta += dl_dmu*dmu_deta;
				 
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = 0.5*exPara[idxn]*exPara[idxr];
				 dl_dloggamma += dl_dmu*dmu_dgamma_p;
			 }
			 
			 for(i=(n0+n1);i<n_subj;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double) zeta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = exPara[idxn]*exPara[idxr]*gamma;
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/

				 dmu_deta = exPara[idxn]*(1.0 - exPara[idxr]);
				 dl_dlogeta += dl_dmu*dmu_deta;
				 
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = exPara[idxn]*exPara[idxr];
				 dl_dloggamma += dl_dmu*dmu_dgamma_p;
			 }
		 } else if(H0 == 1){
			/*****************************************************
			 * Parameters 
			 *****************************************************/
			 
			 eta = 1.0;
			 gamma = exp(para[1]);
			 
			/****************************************************
			 * Gradient Computation
			 ****************************************************/
			 
			 for(i=0;i<n0;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 nu1 = exPara[idxn]*delta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 dmu_dkappa = exPara[idxn]*exPara[idxr];
				 dl_dkappa += dl_dmu*dmu_dkappa;

				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/

				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = 0.0;
				 dl_dloggamma += 0.0;
			 }
			 
			 for(i=n0;i<(n0+n1);i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double)((1.0+zeta)/2.0); 
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = (0.5*exPara[idxn]*exPara[idxr]*(1.0+gamma));
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/

				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = 0.5*exPara[idxn]*exPara[idxr];
				 dl_dloggamma += dl_dmu*dmu_dgamma_p;
			 }
			 
			 for(i=(n0+n1);i<n_subj;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double) zeta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = exPara[idxn]*exPara[idxr]*gamma;
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/

				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/
				 dmu_dgamma_p = exPara[idxn]*exPara[idxr];
				 dl_dloggamma += dl_dmu*dmu_dgamma_p;
			 }
		 } else if(H0 == 2){
		    /*****************************************************
			 * Parameters 
			 *****************************************************/
			 
			 eta = exp(para[1]);
			 gamma = 1.0;
	
			/****************************************************
			 * Gradient Computation
			 ****************************************************/
			 
			 for(i=0;i<n0;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 nu1 = exPara[idxn]*delta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 dmu_dkappa = exPara[idxn]*exPara[idxr];
				 dl_dkappa += dl_dmu*dmu_dkappa;

				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/
				 
				 dmu_deta = 0.0;
				 dl_dlogeta += 0.0;

				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/

			 }
			 
			 for(i=n0;i<(n0+n1);i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double)((1.0+zeta)/2.0); 
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = (0.5*exPara[idxn]*exPara[idxr]*(1.0+gamma));
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/
				
				 dmu_deta = 0.5*exPara[idxn]*(1.0 - exPara[idxr]);
				 dl_dlogeta += dl_dmu*dmu_deta;

				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/

			 }
			 
			 for(i=(n0+n1);i<n_subj;i++){
				 idxc = 9+i;
				 idxr = idxc+n_subj;
				 idxn = idxr+n_subj;
				 
				 yi = ceil(exPara[idxc]-0.5);
				 
				 delta = (double)(1.0 - exPara[idxr] + exPara[idxr]*kappa);
				 cs    = (double)((1.0 - exPara[idxr])/delta);
				 zeta  = (double)(cs*eta + (1.0 - cs)*gamma);
				 
				 nu1 = exPara[idxn]*delta;
				 nu1 *= (double) zeta;
				 
				 /*********************************************************
				  * dl.dmu   
				  *********************************************************/
				 
				 dl_dmu = (double) yi/nu1;
				 dl_dmu -= (double)(1.0+phi*yi)/(1.0+phi*nu1);
				 
				 /*********************************************************
				  * dl.dkappa
				  *********************************************************/
				 
				 dmu_dkappa = exPara[idxn]*exPara[idxr]*gamma;
				 dl_dkappa += dl_dmu*dmu_dkappa;
				 
				/*********************************************************
				 * dl.dlogeta  
				 *********************************************************/

				 dmu_deta = exPara[idxn]*(1.0 - exPara[idxr]);
				 dl_dlogeta += dl_dmu*dmu_deta; 
				 
				/*********************************************************
				 * dl.dloggamma   
				 *********************************************************/

			 }
		 }
		 
		/************************************************************
		 * FINISH GRADIENT
		 ************************************************************/
		 dl_dkappa    *= kappa;
		 dl_dlogeta   *= eta;
	     dl_dloggamma *= kappa*gamma;
		 
		 gr[0] = -dl_dkappa;
		 
		 if(H0 == 0){
			 gr[1] = -dl_dlogeta;
			 gr[2] = -dl_dloggamma;
		} else if(H0 == 1){
		    gr[1] = -dl_dloggamma;
		} else if(H0 == 2){
		    gr[1] = -dl_dlogeta;
		}
	 }
 }
 
 /*****************************************************************
  * (4) grad.Hess.wrapper
  *****************************************************************/
 /* void neg_gradTReC_TS_nhess(int n,double* para,double* gr,void* ex){
	 double* Hessian;
	 neg_gradTReC_TS(n,para,gr,Hessian,ex);
 }
 
  void neg_gradTReC_TS_hess(int n,double* para,double* gr, 
						    double* Hessian,void* ex){
	 neg_gradTReC_TS(n,para,gr,Hessian,ex);
 } */
