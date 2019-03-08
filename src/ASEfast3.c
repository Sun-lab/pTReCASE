#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <time.h>

/*****************************************************************
 * ASE Likelihoods and Gradients
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
 * (1) loglikBB
 ****************************************************************/
  
 double loglikBB(int n, int* nB, int* nTotal, double* pis, double psi){
  /* Declare Variables */
  double logL=0.0, pi1, logLi = 0.0, ai_pi1;
  int i, ni, nBi;
  
  /* Function Body */
  for(i=0;i<n;i++){
	  ni  = nTotal[i];
	  nBi = nB[i];
	  pi1 = pis[i];
	  ai_pi1 = 1.0 - pi1;
	  
	  logLi = lchoose(ni,nBi);
	  
	  logLi += lgammafn(nBi+pi1/psi)+lgammafn(ni-nBi+ai_pi1/psi)+lgammafn(1.0/psi)-
			   lgammafn(ni+1.0/psi)-lgammafn(pi1/psi)-lgammafn(ai_pi1/psi);
	
	  logL += logLi;
  }
  
  return(logL);
 }
  
/****************************************************************
 * (2) neg.logLASE.TS
 ****************************************************************/
  
 double neg_logLASE_TS(int n, double* para, void* ex, SEXP x1){
	/* Declare Variables */
	double kappa = 0.0, eta = 0.0, gamma = 0.0, psi, delta, cs, zeta,
	       logL = 0.0, logLi = 0.0, pi1 = 0.0, ai_pi1 = 0.0;
	double* exPara;
	int n_subj, nAS, H0, idxnt, idxnb, idxra,i,
	    nBi, ni;
	
	/* Extra Parameters */
	exPara = (double *) ex;
	n_subj = ceil(exPara[0]-0.5);
	nAS    = ceil(exPara[1]-0.5);
	H0     = ceil(exPara[5]-0.5);
	psi    = exPara[8];
	
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
	
	for(i=0;i<nAS;i++){
		idxnt = 9+3*n_subj+i;
		idxnb = idxnt+nAS;
		idxra = idxnb+nAS;
		
		delta  = (double)((1.0-exPara[idxra])+exPara[idxra]*kappa);
		cs     = (double)((1.0-exPara[idxra])/delta);
		zeta   = (double)(cs*eta+(1.0-cs)*gamma);
		pi1 = (double)(zeta/(zeta+1.0));
		
		ni  = ceil(exPara[idxnt]-0.5);
		nBi = ceil(exPara[idxnb]-0.5); 
		
		// Likelihood Computations
		ai_pi1 = 1.0 - pi1;
		
		logLi = lchoose(ni,nBi);
		
		logLi += lgammafn(nBi+pi1/psi)+lgammafn(ni-nBi+ai_pi1/psi)+lgammafn(1.0/psi)-
		  lgammafn(ni+1.0/psi)-lgammafn(pi1/psi)-lgammafn(ai_pi1/psi);
		
		logL += logLi;
	}
	
	return(-logL);
 }

/****************************************************************
 * (3) neg.logLASE.TS
 ****************************************************************/
 void neg_gradASE_TS(int n, double* para, double* gr, double* Hessian,
					 void* ex, SEXP x1){
   /* Initialize the Variables */
	double kappa = 0.0, eta = 0.0, gamma = 0.0, psi, delta, cs, zeta, pi_i,
		   dpi_dzeta, dl_dpi, tmp_val, dzeta_dkappa, dl_dkappa = 0.0,
		   dl_deta = 0.0, dl_dgamma = 0.0, d2l_dpi2, d2pi_dzeta2,
		   d2zeta_dkappa2, dcs_dkappa, dpi_dzeta_sq, d2lxdpi_sq,
		   dlxd2pi, dlxdpi, vec_d2l_dl, sum1 = 0.0, sum2 = 0.0,
		   sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, 
		   sum7 = 0.0, sum8 = 0.0, sum9 = 0.0, sum10 = 0.0,
		   dzeta_dkappa_sq = 0.0,d2l_dkappa2 = 0.0,
		   pi_div_psi = 0.0, aipi_div_psi = 0.0;
	double* exPara;
	int n_subj, nAS, Hess, H0, idxnt, idxnb, idxra, ni, nBi,
		i, j, k, l, m;
		
   /* Extra Parameters */
	exPara = (double *) ex;
	n_subj = ceil(exPara[0]-0.5);
	nAS    = ceil(exPara[1]-0.5);
	H0     = ceil(exPara[5]-0.5);
	psi    = exPara[8];
	Hess   = ceil(exPara[6]-0.5);
	
   /* Speed Driven Rewrite */
   /* Evaluate conditions only once */
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
	
	if(Hess == 1){
		// I don't use the Hessian, so the method of calculation here 
		// will not be updated.
		if(H0 == 0){
			for(i=0;i<nAS;i++){
				idxnt = 9+3*n_subj+i;
				idxnb = idxnt+nAS;
				idxra = idxnb+nAS;
				
				delta = (double)((1.0-exPara[idxra])+exPara[idxra]*kappa);
				cs    = (double)((1.0-exPara[idxra])/delta);
				zeta  = (double)(cs*eta+(1.0-cs)*gamma);
				pi_i  = (double)(zeta/(1.0+zeta));
				
				dpi_dzeta = (double)(1.0/pow((1.0+zeta),2));
				
			   /*********************************************************
				* dl.dpi
				*********************************************************/
				
				dl_dpi = 0.0;
				
				ni  = ceil(exPara[idxnt]-0.5);
				nBi = ceil(exPara[idxnb]-0.5); 
				
				if(nBi>0){
					for(j=0;j<nBi;j++){
						tmp_val = (double)j*psi;
						dl_dpi += (double)(1.0/(pi_i+tmp_val));
					}
				}
				
				if(nBi<ni){
					for(k=0;k<(ni-nBi);k++){
						tmp_val = (double)k*psi;
						dl_dpi -= (double)(1.0/(1.0-pi_i+tmp_val));
					}
				}
				
			   /*********************************************************
				* dl_d_logkappa
				*********************************************************/
				
				dzeta_dkappa = (1.0-exPara[idxra])*exPara[idxra]/pow(delta,2);
				dl_dkappa += dl_dpi*dpi_dzeta*dzeta_dkappa;
				
			   /*********************************************************
				* dl_d_logeta
				*********************************************************/
				
				dl_deta += dl_dpi*dpi_dzeta*cs;
				
			   /*********************************************************
				* dl_d_loggamma
				*********************************************************/
				
				dl_dgamma += dl_dpi*dpi_dzeta*(1.0-cs);
				
			   /*********************************************************
			    * Hessian Computations
				*********************************************************/
			   /****************************************************
				* d2l_dpi2   
				****************************************************/
				
				d2l_dpi2 = 0.0;
				
				if(nBi>0){
					for(l=0;l<nBi;l++){
						tmp_val = (double)l*psi;
						d2l_dpi2 -= (double)(1.0/pow((pi_i+tmp_val),2));
					}
				}
				
				if(nBi<ni){
					for(m=0;m<(ni-nBi);m++){
						tmp_val = (double)m*psi;
						d2l_dpi2 -= (double)(1.0/pow((1.0-pi_i+tmp_val),2));
					}
				}

			   /****************************************************
				* Other Terms   
				****************************************************/
				
				d2pi_dzeta2     = (double)(-2.0/pow((1.0+zeta),3));
				d2zeta_dkappa2  = (double)(-2.0*(gamma-eta)*dzeta_dkappa*exPara[idxra]/delta);;
				dcs_dkappa      = (double)(-(1.0-exPara[idxra])*exPara[idxra]/pow(delta,2));

				dpi_dzeta_sq    = dpi_dzeta*dpi_dzeta;
				dzeta_dkappa_sq = (gamma-eta)*(gamma-eta)*dzeta_dkappa*dzeta_dkappa;
				
				d2lxdpi_sq = d2l_dpi2*dpi_dzeta_sq;
				dlxd2pi    = dl_dpi*d2pi_dzeta2;
				dlxdpi     = dl_dpi*dpi_dzeta;
				vec_d2l_dl = d2lxdpi_sq+dlxd2pi;
				
				sum1  += d2lxdpi_sq*dzeta_dkappa_sq;
				sum2  += dlxd2pi*dzeta_dkappa_sq;
				sum3  += dlxdpi*d2zeta_dkappa2; 
				
				sum4  += cs*cs*vec_d2l_dl;
				sum5  += vec_d2l_dl*dzeta_dkappa*cs;
				sum6  += dlxdpi*dcs_dkappa;
							
				sum7  += (1.0-cs)*(1.0-cs)*vec_d2l_dl;
				sum8  += vec_d2l_dl*dzeta_dkappa*(1.0-cs);
				sum9  += dlxdpi*dcs_dkappa;
				
				sum10 += cs*(1.0-cs)*vec_d2l_dl;
			}
		} else if(H0 == 1){
			for(i=0;i<nAS;i++){
				idxnt = 9+3*n_subj+i;
				idxnb = idxnt+nAS;
				idxra = idxnb+nAS;
				
				delta = (double)((1.0-exPara[idxra])+exPara[idxra]*kappa);
				cs    = (double)((1.0-exPara[idxra])/delta);
				zeta  = (double)(cs*eta+(1.0-cs)*gamma);
				pi_i  = (double)(zeta/(1.0+zeta));
				
				dpi_dzeta = (double)(1.0/pow((1.0+zeta),2));
				
			   /*********************************************************
				* dl.dpi
				*********************************************************/
				
				dl_dpi = 0.0;
				
				ni  = ceil(exPara[idxnt]-0.5);
				nBi = ceil(exPara[idxnb]-0.5); 
				
				if(nBi>0){
					for(j=0;j<nBi;j++){
						tmp_val = (double)j*psi;
						dl_dpi += (double)(1.0/(pi_i+tmp_val));
					}
				}
				
				if(nBi<ni){
					for(k=0;k<(ni-nBi);k++){
						tmp_val = (double)k*psi;
						dl_dpi -= (double)(1.0/(1.0-pi_i+tmp_val));
					}
				}
				
			   /*********************************************************
				* dl_d_logkappa
				*********************************************************/
				
				dzeta_dkappa = (1.0-exPara[idxra])*exPara[idxra]/pow(delta,2);
				dl_dkappa += dl_dpi*dpi_dzeta*dzeta_dkappa;
				
			   /*********************************************************
				* dl_d_logeta
				*********************************************************/

			   /*********************************************************
				* dl_d_loggamma
				*********************************************************/
				
				dl_dgamma += dl_dpi*dpi_dzeta*(1.0-cs);
				
			   /*********************************************************
			    * Hessian Computations
				*********************************************************/
			   /****************************************************
				* d2l_dpi2   
				****************************************************/
				
				d2l_dpi2 = 0.0;
				
				if(nBi>0){
					for(l=0;l<nBi;l++){
						tmp_val = (double)l*psi;
						d2l_dpi2 -= (double)(1.0/pow((pi_i+tmp_val),2));
					}
				}
				
				if(nBi<ni){
					for(m=0;m<(ni-nBi);m++){
						tmp_val = (double)m*psi;
						d2l_dpi2 -= (double)(1.0/pow((1.0-pi_i+tmp_val),2));
					}
				}

			   /****************************************************
				* Other Terms   
				****************************************************/
				
				d2pi_dzeta2     = (double)(-2.0/pow((1.0+zeta),3));
				d2zeta_dkappa2  = (double)(-2.0*(gamma-eta)*dzeta_dkappa*exPara[idxra]/delta);;
				dcs_dkappa      = (double)(-(1.0-exPara[idxra])*exPara[idxra]/pow(delta,2));

				dpi_dzeta_sq    = dpi_dzeta*dpi_dzeta;
				dzeta_dkappa_sq = (gamma-eta)*(gamma-eta)*dzeta_dkappa*dzeta_dkappa;
				
				d2lxdpi_sq = d2l_dpi2*dpi_dzeta_sq;
				dlxd2pi    = dl_dpi*d2pi_dzeta2;
				dlxdpi     = dl_dpi*dpi_dzeta;
				vec_d2l_dl = d2lxdpi_sq+dlxd2pi;
				
				sum1  += d2lxdpi_sq*dzeta_dkappa_sq;
				sum2  += dlxd2pi*dzeta_dkappa_sq;
				sum3  += dlxdpi*d2zeta_dkappa2; 
							
				sum7  += (1.0-cs)*(1.0-cs)*vec_d2l_dl;
				sum8  += vec_d2l_dl*dzeta_dkappa*(1.0-cs);
				sum9  += dlxdpi*dcs_dkappa;
			}
		} else if(H0 == 2){
			for(i=0;i<nAS;i++){
				idxnt = 9+3*n_subj+i;
				idxnb = idxnt+nAS;
				idxra = idxnb+nAS;		
				
				delta = (double)((1.0-exPara[idxra])+exPara[idxra]*kappa);
				cs    = (double)((1.0-exPara[idxra])/delta);
				zeta  = (double)(cs*eta+(1.0-cs)*gamma);
				pi_i  = (double)(zeta/(1.0+zeta));
				
				dpi_dzeta = (double)(1.0/pow((1.0+zeta),2));
				
			   /*********************************************************
				* dl.dpi
				*********************************************************/
				
				dl_dpi = 0.0;
				
				ni  = ceil(exPara[idxnt]-0.5);
				nBi = ceil(exPara[idxnb]-0.5); 
				
				if(nBi>0){
					for(j=0;j<nBi;j++){
						tmp_val = (double)j*psi;
						dl_dpi += (double)(1.0/(pi_i+tmp_val));
					}
				}
				
				if(nBi<ni){
					for(k=0;k<(ni-nBi);k++){
						tmp_val = (double)k*psi;
						dl_dpi -= (double)(1.0/(1.0-pi_i+tmp_val));
					}
				}
				
			   /*********************************************************
				* dl_d_logkappa
				*********************************************************/
				
				dzeta_dkappa = (1.0-exPara[idxra])*exPara[idxra]/pow(delta,2);
				dl_dkappa += dl_dpi*dpi_dzeta*dzeta_dkappa;
				
			   /*********************************************************
				* dl_d_logeta
				*********************************************************/
				
				dl_deta += dl_dpi*dpi_dzeta*cs;
				
			   /*********************************************************
				* dl_d_loggamma
				*********************************************************/
				
			   /*********************************************************
			    * Hessian Computations
				*********************************************************/
			   /****************************************************
				* d2l_dpi2   
				****************************************************/
				
				d2l_dpi2 = 0.0;
				
				if(nBi>0){
					for(l=0;l<nBi;l++){
						tmp_val = (double)l*psi;
						d2l_dpi2 -= (double)(1.0/pow((pi_i+tmp_val),2));
					}
				}
				
				if(nBi<ni){
					for(m=0;m<(ni-nBi);m++){
						tmp_val = (double)m*psi;
						d2l_dpi2 -= (double)(1.0/pow((1.0-pi_i+tmp_val),2));
					}
				}

			   /****************************************************
				* Other Terms   
				****************************************************/
				
				d2pi_dzeta2     = (double)(-2.0/pow((1.0+zeta),3));
				d2zeta_dkappa2  = (double)(-2.0*(gamma-eta)*dzeta_dkappa*exPara[idxra]/delta);;
				dcs_dkappa      = (double)(-(1.0-exPara[idxra])*exPara[idxra]/pow(delta,2));

				dpi_dzeta_sq    = dpi_dzeta*dpi_dzeta;
				dzeta_dkappa_sq = (gamma-eta)*(gamma-eta)*dzeta_dkappa*dzeta_dkappa;
				
				d2lxdpi_sq = d2l_dpi2*dpi_dzeta_sq;
				dlxd2pi    = dl_dpi*d2pi_dzeta2;
				dlxdpi     = dl_dpi*dpi_dzeta;
				vec_d2l_dl = d2lxdpi_sq+dlxd2pi;
				
				sum1  += d2lxdpi_sq*dzeta_dkappa_sq;
				sum2  += dlxd2pi*dzeta_dkappa_sq;
				sum3  += dlxdpi*d2zeta_dkappa2; 
				
				sum4  += cs*cs*vec_d2l_dl;
				sum5  += vec_d2l_dl*dzeta_dkappa*cs;
				sum6  += dlxdpi*dcs_dkappa;
			}
		}
		
	   /*************************************************************
		* Finish the Hessian
		*************************************************************/
		
		d2l_dkappa2 = kappa*(gamma-eta)*dl_dkappa+kappa*kappa*(sum1+sum2+sum3);
		
		if(H0 == 0){
			Hessian[0] = -d2l_dkappa2;
			Hessian[1] = -kappa*eta*((gamma-eta)*sum5+sum6);
			Hessian[2] = -kappa*gamma*((gamma-eta)*sum8-sum9);
			Hessian[3] = Hessian[1];
			Hessian[4] = -(eta*dl_deta+eta*eta*sum4);
			Hessian[5] = -eta*gamma*sum10;
			Hessian[6] = Hessian[2];
			Hessian[7] = Hessian[5];
			Hessian[8] = -(gamma*dl_dgamma+gamma*gamma*sum7);
		} else if(H0 == 1){
			Hessian[0] = -d2l_dkappa2;
			Hessian[1] = -kappa*gamma*((gamma-eta)*sum8-sum9);
			Hessian[2] = Hessian[1];
			Hessian[3] = -(gamma*dl_dgamma+gamma*gamma*sum7);
		} else if(H0 == 2){
			Hessian[0] = -d2l_dkappa2;
			Hessian[1] = -kappa*eta*((gamma-eta)*sum5+sum6);
			Hessian[2] = Hessian[1];
			Hessian[3] = -(eta*dl_deta+eta*eta*sum4);
		}
	} else{
		if(H0 == 0){		
			for(i=0;i<nAS;i++){
				idxnt = 9+3*n_subj+i;
				idxnb = idxnt+nAS;
				idxra = idxnb+nAS;
				
				delta = (double)((1.0-exPara[idxra])+exPara[idxra]*kappa);
				cs    = (double)((1.0-exPara[idxra])/delta);
				zeta  = (double)(cs*eta+(1.0-cs)*gamma);
				pi_i  = (double)(zeta/(1.0+zeta));
				
				dpi_dzeta = (double)(1.0/pow((1.0+zeta),2));
				
			   /*********************************************************
				* dl.dpi
				*********************************************************/
				
				dl_dpi = 0.0;
				
				ni  = ceil(exPara[idxnt]-0.5);
				nBi = ceil(exPara[idxnb]-0.5); 
				
				pi_div_psi   = pi_i/psi;
				aipi_div_psi = (1.0-pi_i)/psi; 
				
				dl_dpi = (1.0/psi)*(digamma(nBi+pi_div_psi)-digamma(ni-nBi+aipi_div_psi)-
									digamma(pi_div_psi)+digamma(aipi_div_psi));
				
			   /*********************************************************
				* dl_d_logkappa
				*********************************************************/
				
				dzeta_dkappa = (1.0-exPara[idxra])*exPara[idxra]/pow(delta,2);
				dl_dkappa += dl_dpi*dpi_dzeta*dzeta_dkappa;
				
			   /*********************************************************
				* dl_d_logeta
				*********************************************************/
				
				dl_deta += dl_dpi*dpi_dzeta*cs;
				
			   /*********************************************************
				* dl_d_loggamma
				*********************************************************/
				
				dl_dgamma += dl_dpi*dpi_dzeta*(1.0-cs);
			}
		} else if(H0 == 1){
			for(i=0;i<nAS;i++){
				idxnt = 9+3*n_subj+i;
				idxnb = idxnt+nAS;
				idxra = idxnb+nAS;
				
				delta = (double)((1.0-exPara[idxra])+exPara[idxra]*kappa);
				cs    = (double)((1.0-exPara[idxra])/delta);
				zeta  = (double)(cs*eta+(1.0-cs)*gamma);
				pi_i  = (double)(zeta/(1.0+zeta));
				
				dpi_dzeta = (double)(1.0/pow((1.0+zeta),2));
				
			   /*********************************************************
				* dl.dpi
				*********************************************************/
				
				dl_dpi = 0.0;
				
				ni  = ceil(exPara[idxnt]-0.5);
				nBi = ceil(exPara[idxnb]-0.5); 
				
				pi_div_psi   = pi_i/psi;
				aipi_div_psi = (1.0-pi_i)/psi; 
				
				dl_dpi = (1.0/psi)*(digamma(nBi+pi_div_psi)-digamma(ni-nBi+aipi_div_psi)-
									digamma(pi_div_psi)+digamma(aipi_div_psi));
				
			   /*********************************************************
				* dl_d_logkappa
				*********************************************************/
				
				dzeta_dkappa = (1.0-exPara[idxra])*exPara[idxra]/pow(delta,2);
				dl_dkappa += dl_dpi*dpi_dzeta*dzeta_dkappa;
				
			   /*********************************************************
				* dl_d_logeta
				*********************************************************/

			   /*********************************************************
				* dl_d_loggamma
				*********************************************************/
				
				dl_dgamma += dl_dpi*dpi_dzeta*(1.0-cs);
			}
		} else if(H0 == 2){
			for(i=0;i<nAS;i++){
				idxnt = 9+3*n_subj+i;
				idxnb = idxnt+nAS;
				idxra = idxnb+nAS;			
				
				delta = (double)((1.0-exPara[idxra])+exPara[idxra]*kappa);
				cs    = (double)((1.0-exPara[idxra])/delta);
				zeta  = (double)(cs*eta+(1.0-cs)*gamma);
				pi_i  = (double)(zeta/(1.0+zeta));
				
				dpi_dzeta = (double)(1.0/pow((1.0+zeta),2));
				
			   /*********************************************************
				* dl.dpi
				*********************************************************/
				
				dl_dpi = 0.0;
				
				ni  = ceil(exPara[idxnt]-0.5);
				nBi = ceil(exPara[idxnb]-0.5); 
				
				pi_div_psi   = pi_i/psi;
				aipi_div_psi = (1.0-pi_i)/psi; 
				
				dl_dpi = (1.0/psi)*(digamma(nBi+pi_div_psi)-digamma(ni-nBi+aipi_div_psi)-
									digamma(pi_div_psi)+digamma(aipi_div_psi));
				
			   /*********************************************************
				* dl_d_logkappa
				*********************************************************/
				
				dzeta_dkappa = (1.0-exPara[idxra])*exPara[idxra]/pow(delta,2);
				dl_dkappa += dl_dpi*dpi_dzeta*dzeta_dkappa;
				
			   /*********************************************************
				* dl_d_logeta
				*********************************************************/
				
				dl_deta += dl_dpi*dpi_dzeta*cs;
				
				
			}
		}
	}
	
   /*************************************************************
    * Finish the Gradient
	*************************************************************/
	
	dl_dkappa *= kappa*(gamma-eta);
	dl_deta   *= eta;
	dl_dgamma *= gamma;
	
	if(H0 == 0){
		gr[0] = -dl_dkappa;
		gr[1] = -dl_deta;
		gr[2] = -dl_dgamma;
	} else if(H0 == 1){
		gr[0] = -dl_dkappa;
		gr[1] = -dl_dgamma;
	} else if(H0 == 2){
		gr[0] = -dl_dkappa;
		gr[1] = -dl_deta;
	}
 }
 
/****************************************************************
 * (4) Compute the Pis of logLASE.psi
 ****************************************************************/  
void compute_pis(double kappa, double eta, double gamma, double* rhos,
				 double* ex){
	// Declare Variables 
	int nAS = 0, n0 = 0, n1 = 0, i=0;
	double delta = 0.0, cs = 0.0, zeta = 0.0;
	double* exPara;
	
	/* Extra Parameters */
	exPara = (double *) ex;
	nAS    = ceil(exPara[0]-0.5);
	n0     = ceil(exPara[1]-0.5);
	n1     = ceil(exPara[2]-0.5);
	
	// Compute and Store Pis
	int dbl_nASp4 = 4+2*nAS;
	for(i=0;i<n0;i++){
		ex[dbl_nASp4+i] = 0.5;
	}
	for(i=n0;i<(n0+n1);i++){
		delta = (double)((1.0 - rhos[i])+rhos[i]*kappa);
		cs    = (double)((1.0 - rhos[i])/delta);
		zeta  = (double)(cs*eta+(1.0-cs)*gamma);
		ex[dbl_nASp4+i]  = (double)(zeta/(zeta+1.0));
	}
	for(i=(n0+n1);i<nAS;i++){
		ex[dbl_nASp4+i] = 0.5;
	}
}

/****************************************************************
 * (5) logLASE.psi
 ****************************************************************/ 
double logLASE_psi(int n, double* para, void* ex, SEXP x1){
	// Declare Value
	int nAS_all = 0, i=0, ni=0, nBi = 0;
	double logL = 0.0, logLi=0.0, psi = para[0], pi1 = 0.0, ai_pi1 = 0.0;
	double* exPara;
	
	//Setting current psi value
	psi = exp(para[0]);
	
	// Handling Extras
	exPara = (double *) ex;
	nAS_all = ceil(exPara[0]-0.5);
	
	// Setting indices
	int dbl_nASp4 = 4+2*nAS_all;
	int nASp4 = nAS_all+4;
	
	for(i=0;i<nAS_all;i++){
		ni = ceil(exPara[4+i]-0.5);
		nBi     = ceil(exPara[nASp4+i]-0.5);
		pi1    = exPara[dbl_nASp4+i];
		
	  // Compute the likelihood for the individual:
	  // Gamma digamma used for speed issues in large counts.
	  ai_pi1 = 1.0 - pi1;
	  logLi = lchoose(ni,nBi);
	  
	  logLi += lgammafn(nBi+pi1/psi)+lgammafn(ni-nBi+ai_pi1/psi)+lgammafn(1.0/psi)-
	    lgammafn(ni+1.0/psi)-lgammafn(pi1/psi)-lgammafn(ai_pi1/psi);
	  
	  logL += logLi;
	}
	
	// Return Results
	return(-logL);
} 

/****************************************************************
 * (6) gradASE.psi
 ****************************************************************/
void gradASE_psi(int n, double* para, double* gr, void* ex, SEXP x1){
	/* Declare Variables */
	int i=0, nAS_all=0,
		ni=0, nB = 0;
	double psi = 0.0, pi = 0.0;
	double* exPara;
	
	exPara  = (double *) ex;
	nAS_all = ceil(exPara[0]-0.5);
	int dbl_nASp4 = 4+2*nAS_all;
	int nASp4 = nAS_all+4;
	
	/* Main Body of the Function */
	psi = exp(para[0]);
	
	// Initialize the gradient
	gr[0] = 0.0;
	
	for(i=0;i<nAS_all;i++){
		ni = ceil(exPara[4+i]-0.5);
		nB = ceil(exPara[nASp4+i]-0.5);
		pi = exPara[dbl_nASp4+i];
		
		gr[0] += (1.0/(psi*psi))*(digamma(ni+1.0/psi)+digamma(pi/psi)*pi+digamma((1.0-pi)/psi)*(1.0-pi)-
								  digamma(nB+pi/psi)*pi-digamma(ni-nB+(1.0-pi)/psi)*(1.0-pi)-digamma(1.0/psi));
	}
	
	/* Store Results */
	// Not necessary - gr is updated on its own.
	gr[0] *= -psi;
}

/****************************************************************
 * (7) compute_BBHz
 ****************************************************************/
double compute_BBHz(int nAS_Hz, double* ex, double psi){
	/* Declare Variables */
	double logL = 0.0, logLi = 0.0, pi1 = 0.0, ai_pi1=0.0;
	int nAS = 0, n0 = 0, n1 = 0, nBi, ni,
		idxnt = 0, idxnb = 0, i = 0;
	
	/* Extra Parameters */
	nAS    = ceil(ex[0]-0.5);
	n0     = ceil(ex[1]-0.5);
	n1     = ceil(ex[2]-0.5);
	
	for(i=0;i<n0;i++){
		idxnt = 4+i;
		idxnb = idxnt+nAS;

		pi1 = 0.5;
		ni  = ceil(ex[idxnt]-0.5);
		nBi = ceil(ex[idxnb]-0.5);
		
		// Compute the likelihood:
		ai_pi1 = 1.0 - pi1;
		logLi = lchoose(ni,nBi);
		
		logLi += lgammafn(nBi+pi1/psi)+lgammafn(ni-nBi+ai_pi1/psi)+lgammafn(1.0/psi)-
		  lgammafn(ni+1.0/psi)-lgammafn(pi1/psi)-lgammafn(ai_pi1/psi);
		
		logL += logLi;
	}
	
	for(i=(n0+n1);i<nAS;i++){
		idxnt = 4+i;
		idxnb = idxnt+nAS;

		pi1 = 0.5;
		ni  = ceil(ex[idxnt]-0.5);
		nBi = ceil(ex[idxnb]-0.5);
		
		// Compute the likelihood:
		ai_pi1 = 1.0 - pi1;
		logLi = lchoose(ni,nBi);
		
		logLi += lgammafn(nBi+pi1/psi)+lgammafn(ni-nBi+ai_pi1/psi)+lgammafn(1.0/psi)-
		  lgammafn(ni+1.0/psi)-lgammafn(pi1/psi)-lgammafn(ai_pi1/psi);
		
		logL += logLi;
	}
	
	return(logL);
}
