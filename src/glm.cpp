/*###################################################################
# PROG.TITLE                                                        #
#####################################################################
# PROGRAM NAME:                                                     #
#                                                                   #
# PROGRAMMER:                                                       #
#                                                                   #
# DATE CREATED:                                                     #
#                                                                   #
# LAST EDIT:                                                        #
#                                                                   #
# VERSION:                                                          #
#                                                                   #
#-------------------------------------------------------------------#
# DESCRIPTION:                                                      #
#   C program to fit a negative binomial GLM. Removes the un-       #
#   necessary checks for Poisson and allows me to extract the       #
#   regression coefficients--something that is not done by Dr. Sun's#
#   version of the code.                                            #                                                               #
###################################################################*/

#include <float.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <RcppEigen.h>
#include <time.h>
#include "TReC.h"

//#define DEBUG 1

 /*****************************************************************
  * (5) Compute nu0
  *****************************************************************/
void compute_nu0(double* X, int n, int M, double* Beta, double* ex){
	// Declare Eigen:
	Eigen::Map<const Eigen::MatrixXd> Xmat(X,n,M);
	Eigen::Map<const Eigen::VectorXd> Betav(Beta,M);
	Eigen::VectorXd Nuvec(n);
	
	int i=0;
	
	// Compute:
	Nuvec = (Xmat*Betav).array().exp();
	
	// Store:
	int dbl_np9 = 9+2*n;
	for(i=0;i<n;i++){
		ex[dbl_np9+i] = Nuvec(i);
	}
}

/********************************************************************
 * Compute Offset                                                   *
 ********************************************************************/
void comp_offset(double kappa, double eta, double gamma, double* ex1, double* offset){
	// Declare Things
	double delta=0.0, cs=0.0, zeta = 0.0;
	int n_subj=0, n0=0, n1=0,i=0,
		idxr=0;
	
	// Parameters
	n_subj = ceil(ex1[0]-0.5);
	n0 = ceil(ex1[2]-0.5);
	n1 = ceil(ex1[3]-0.5);
	
	// Function Body
	// ex1 contains our (sorted) data!
	for(i=0;i<n0;i++){
		idxr = 9+n_subj+i;
		
		delta = 1.0 - ex1[idxr]+ex1[idxr]*kappa;
		offset[i] = log(delta);
	}
	
	for(i=n0;i<(n0+n1);i++){
		idxr = 9+n_subj+i;
				
		delta = 1.0 - ex1[idxr]+ex1[idxr]*kappa;
		cs = (1.0-ex1[idxr])/delta;
		zeta = cs*eta+(1.0-cs)*gamma;
		offset[i] = log(delta)+log(((zeta+1.0)/2.0));
	}
	
	for(i=(n0+n1);i<n_subj;i++){
		idxr = 9+n_subj+i;
		
		delta = 1.0 - ex1[idxr]+ex1[idxr]*kappa;
		cs = (1.0-ex1[idxr])/delta;
		zeta = cs*eta+(1.0-cs)*gamma;
		offset[i] = log(delta)+log(zeta);
	}
}

/********************************************************************
 * Initialize Phi                                                   *
 ********************************************************************/
void phi_initialize(double* y, double* X, int n, int M, double* offset,
					  double* Beta, double* phi_curr){
	// Declare C++
	double phi_out = 0.0,
		   minTheta = 1e-5,
		   maxTheta = (1/minTheta);
	
	// Eigen Initialize
	Eigen::Map<const Eigen::MatrixXd> Xmat(X,n,M);
	Eigen::Map<const Eigen::VectorXd> Yvec(y,n);
	Eigen::Map<const Eigen::VectorXd> Offvec(offset,n);
	Eigen::Map<const Eigen::VectorXd> Betav(Beta,M);
	Eigen::VectorXd Etavec(n);
	Eigen::VectorXd Resid(n);
	Eigen::VectorXd Fitvec(n);
	
	// Compute Necessaries:
	Etavec = Xmat*Betav;
	Fitvec = (Etavec+Offvec).array().exp();
	Resid = Yvec - Fitvec;
	
	// Compute Phi:
	int i=0, n_used = 0;
	for(i=0;i<n;i++){
	  if(Fitvec(i)<1e-5){
	  } else {
	    n_used += 1;
	    phi_out += (((Resid(i)*Resid(i)/Fitvec(i))-1.0)/Fitvec(i)); 
	  }
	}
	
	if(n_used<=M){
	  Rprintf("Strange Cluster! Fewer subjects with initial Predicted Count < 1e-5 than covariates used!\n");
	  phi_out = 2.0;
	} else {
	  phi_out /= (n_used-M);
	}
	
	if(phi_out<minTheta){
		phi_out = minTheta+0.001;
	} else if(phi_out>maxTheta){
		phi_out = maxTheta-0.001;
	}
		
	// Return Results:
	*phi_curr = phi_out;
}

/********************************************************************
 * Update Phi                                                       *
 ********************************************************************/
double phi_update(double phi_init,int n_subj,const double* y, double* mu, int limit, double eps){
	// Declare Things:
	double gradL = 0.0, hessL=0.0, vphi = (1.0/phi_init),
		   vmui = 0.0,
		   vci = 0.0,
		   par_change = 1.0,
		   minTheta = 1e-5,
		   maxTheta = (1.0/minTheta),
		   ci = 0.0,
		   del = 0.0;
	int iter = 0,i=0;
	
	// Function Body:
	while((iter<limit) && (par_change>eps)){
		for(i=0;i<n_subj;i++){
			ci = y[i];
			vmui = vphi+mu[i];
			vci  = ci+vphi;
			
			gradL += R::digamma(vci)-R::digamma(vphi)+log(vphi/vmui);
			gradL += 1.0 - (vci)/(vmui);
			
			hessL += R::trigamma(vci)-R::trigamma(vphi)+(1.0/vphi);
			hessL += -(2.0/vmui)+((vci)/(vmui*vmui));
		}
		
		del = gradL/hessL;
		par_change = abs(del);
		
		vphi -= del;
		
		if(vphi>maxTheta){
			vphi = maxTheta;
			par_change=0.0;
		} else if(vphi<minTheta){
			vphi = minTheta;
			par_change=0.0;
		}
		
		iter += 1;
	}
	
	return((1.0/vphi));
}

/********************************************************************
 * VarFunc                                                          *
 ********************************************************************/
double varfunc(double mu,double phi){
	return(mu+mu*mu*phi);
}
 
 
/********************************************************************
 * dfunc                                                            *
 ********************************************************************/
double dfunc(double mu){
	return(1.0/mu);
}
 
/********************************************************************
 * IRLS algorithm                                                   *
 ********************************************************************/
/******************* NB REGRESSION **********************************/
int IRLS_NB_fit(double* y, int* yint, double* X, int n, int M, double* offset,
				 int maxIter, double eps, double convg, int init_beta,
				 double* phi_out, double* logLNB, double* Beta_out){
	// Declare C++ Values
	double par_change = 1.0, phi_curr = 0.0,
		   minTheta = 0.00001,
		   maxTheta = (1.0/minTheta);
	int iter = 1, i=0; 
	
	// Declare Eigen
	Eigen::Map<const Eigen::MatrixXd> Xmat(X,n,M);
	Eigen::Map<const Eigen::VectorXd> Yvec(y,n);
	Eigen::Map<const Eigen::VectorXd> Offvec(offset,n);
	Eigen::VectorXd Weightvec(n);
	Eigen::VectorXd Etavec(n);
	Eigen::VectorXd Zvec(n);
	Eigen::VectorXd Beta(M);
	Eigen::VectorXd Betatmp(M);
	Eigen::VectorXd Fitvec(n);
	Eigen::VectorXd Resid(n);
	Eigen::MatrixXd Xnew(n,M);
	
	// Beta and Phi should already be initialized at some value
	if(init_beta==0){
		phi_curr = *phi_out;
		for(i=0;i<M;i++){ Beta(i) = Beta_out[i];}
			
		Etavec = Xmat*Beta;
		Fitvec = (Etavec+Offvec).array().exp();
		Resid  = Yvec - Fitvec;
		for(i=0;i<n;i++){
			Weightvec(i) = (1.0/Fitvec(i))+phi_curr;
		}
	} else {
		Resid = ((Yvec.array()+0.01).log()-Offvec.array());
		Beta = Xmat.colPivHouseholderQr().solve(Resid);
		
		Etavec = Xmat*Beta;
		Fitvec = (Etavec+Offvec).array().exp();
		Resid  = Yvec - Fitvec;
		for(i=0;i<n;i++){
			Weightvec(i) = (1.0/Fitvec(i))+phi_curr;
		}
		
		phi_curr = ((((Resid.array()*Resid.array()/Fitvec.array())-1.0)/Fitvec.array()).sum()/(n-M));
		if(phi_curr<minTheta||phi_curr>maxTheta){phi_curr = 0.5;}
	}
	
	*logLNB = loglikNB_pmy(phi_curr,n,Fitvec.data(),yint);
	
	// Function Body
	while(iter < maxIter && par_change>convg){
		// Compute the Updated Z
		Zvec = Resid.array()/Fitvec.array();
		Zvec += Etavec;
		Zvec = Zvec.array()/Weightvec.array().sqrt();
		
		// Compute the Updated X
		// I wonder if this will compile, but it should.
		for(i=0;i<n;i++){Xnew.row(i) = Xmat.row(i).array()/sqrt(Weightvec(i));}

		// QR decomposition
		Betatmp = Beta;
		Beta = Xnew.colPivHouseholderQr().solve(Zvec);
		
		#ifdef DEBUG
			Rprintf("Curr Beta (%d):\n %5.5f, %5.5f, %5.5f\n",iter,Beta(0),Beta(1),Beta(2));
		#endif
			
		// Update Necessaries
		Etavec = Xmat*Beta;
		Fitvec = (Etavec+Offvec).array().exp();
		Resid = Yvec-Fitvec;
		
		// Update Phi
		phi_curr = phi_update(phi_curr,n,Yvec.data(),Fitvec.data(),maxIter,convg);
		
		// Check Weights
		for(i=0;i<n;i++){
			Weightvec(i) = (1.0/Fitvec(i))+phi_curr;
		}
		
		#ifdef DEBUG
			Rprintf("Curr Phi(%d):\n %5.5f\n",iter,phi_curr);
			Rprintf("Curr Fitvec(%d):\n %5.5f, %5.5f, %5.5f.\n",iter, Fitvec(0),Fitvec(1),Fitvec(2));
		#endif
		
		// The Weightvec is akin to variance in the linear model. If weights to small, then
		// poor model fit is likely to be the cause.
		if(Weightvec.sum()<eps){
			Rprintf("GLM NB: Model Fit Failed - Weights too Small");
			return(0);
		}
		
		// Assess the change in Parameters
		*logLNB = loglikNB_pmy(phi_curr,n,Fitvec.data(),yint);
		//par_change = abs(*logLNB-logLNB_tmp)/(abs(logLNB_tmp)+0.005);
		par_change = (Betatmp-Beta).norm()/(Betatmp.norm()+0.001);
		iter+=1;
	}

	// Return Results
	*phi_out = phi_curr;
	for(i=0;i<M;i++){Beta_out[i] = Beta(i);}
	*logLNB = loglikNB_pmy(phi_curr,n,Fitvec.data(),yint);
	return(1);
}

/********************** POISSON REGRESSION **************************/

int IRLS_Pois_fit(double* y, double* X, int n, int M, double* offset,
				 int maxIter, double eps, double convg, double* Beta_out){
	// Declare C++ Values
	double par_change = 5.0;
	int iter = 0, i=0; 
	
	// Declare Eigen 
	Eigen::Map<const Eigen::MatrixXd> Xmat(X,n,M);
	Eigen::Map<const Eigen::VectorXd> Yvec(y,n);
	Eigen::Map<const Eigen::VectorXd> Offvec(offset,n);
	Eigen::VectorXd SSQvec(n);
	Eigen::VectorXd Weightvec(n);
	Eigen::VectorXd Etavec(n);
	Eigen::VectorXd Zvec(n);
	Eigen::VectorXd Beta(M);
	Eigen::VectorXd Betatmp(M);
	Eigen::VectorXd Fitvec(n);
	Eigen::VectorXd Resid(n);
	Eigen::MatrixXd Xnew(n,M);
		
	// Beta should already be initialized at some value
	Resid = ((Yvec.array()+0.01).log()-Offvec.array());
	Beta = Xmat.colPivHouseholderQr().solve(Resid);
	//for(i=0;i<M;i++){ Beta(i) = Beta_out[i];}
	Etavec = Xmat*Beta;
	Fitvec = (Etavec+Offvec).array().exp();
	Resid  = Yvec - Fitvec;

	for(i=0;i<n;i++){
		Weightvec(i) = dfunc(Fitvec(i));
	}
	
	// Function Body
	while(iter < maxIter && par_change>convg){
		// Compute the Updated Z
		Zvec = Resid.array()/Fitvec.array();
		Zvec += Etavec;
		Zvec = Zvec.array()/Weightvec.array().sqrt();
		
		// Compute the Updated X
		// I wonder if this will compile, but it should.
		for(i=0;i<n;i++){Xnew.row(i) = Xmat.row(i).array()/sqrt(Weightvec(i));}
		
		// QR decomposition
		Betatmp = Beta;
		Beta = Xnew.colPivHouseholderQr().solve(Zvec);
		
		// Update Necessaries
		Etavec = Xmat*Beta;
		Fitvec = (Etavec+Offvec).array().exp();
		Resid = Yvec-Fitvec;
		for(i=0;i<n;i++){
			Weightvec(i) = dfunc(Fitvec(i));
		}

		// Check Weights
		if(Weightvec.sum()<eps){
			Rprintf("GLM Pois: Model Fit Failed - Weights too Small");
			return(0);
		}
		
		// Assess the change in Parameters		
		par_change = (Beta-Betatmp).norm();
		iter++;
	}

	// Return Results
	for(i=0;i<M;i++){Beta_out[i] = Beta(i);}
	return(1);
}

/**********************************************************/
// EXPECTATION FUNCTION                                   //
/***********************************************************/
/*
NAME: ASE_expFunc()
DESC: Computes the Necessary expected values from ASE for a single
subject.

INPUT DIAGRAM:

-------------------------------------------------------------
|VARIABLE          || TYPE      || DESCRIPTION              |
-------------------------------------------------------------
| ctvec            || double*   || vector 0:maxAS           |
_____________________________________________________________
| lgct             || double*   || log(ctvec!)              |
_____________________________________________________________
| rct_vec          || double*   || vector maxAS:0           |
_____________________________________________________________
| lgrct            || double*   || log(rct_vec!)            |
_____________________________________________________________
| pvec             || double*   || Placeholder vector for   |
|                  ||           || probabilities            |
_____________________________________________________________
| pi               || double    || Subject-specific pi      |
_____________________________________________________________
| Toti             || int       || Total number of AS       |
|                  ||           || reads for a subject      |
_____________________________________________________________
| vpsi             || double    || value of 1/psi           |
_____________________________________________________________
| outvec           || double*   || Address of output storage|
|                  ||           || for a single subject     |
_____________________________________________________________

OUTPUT DIAGRAM: 
(VOID) Function. Outvec serves as the return.

outvec[0] = Expected value of digamma(vpsi*pi+nBi)
outvec[1] = Expected value of digamma(vpsi*(1.0-pi)+nTi-nBi)
outvec[2] = Expected value of trigamma(vpsi*pi+nBi)
outvec[3] = Expected value of trigamma(vpsi*(1.0-pi)+nTi-nBi)

This is a single row of an output matrix which contains these
four values for all subjects in the ASE expression dataset. 

CODE:
*/

double ArrLgamma(double x){
  return R::lgammafn(x);
}

double ArrDigamma(double x){
  return R::digamma(x);
}

double ArrTrigamma(double x){
  return R::trigamma(x);
}

void ASE_ExpFunc(double* ctvec, double* lgct, double* rct_vec, double* lgrct,
                 double* pvec, double pi, int Toti, double vpsi, double* outvec,
                 double* tmpctvec_, double* tmprctvec_, int maxAS){
  // Commonly Used Values //
  int    Tot_p1 = Toti+1;
  double ai_pi, piv, ai_piv;
  
  ai_pi  = 1.0-pi;
  piv    = pi*vpsi;
  ai_piv = ai_pi*vpsi;
  
  // Eigen Declarations  //
  // Utilize this format as I can create two single vectors and minimize
  // redundant calculation across subjects.
  Eigen::Map<const Eigen::VectorXd> ctvec_e(ctvec,Tot_p1);
  Eigen::Map<const Eigen::VectorXd> lgct_e(lgct,Tot_p1);
  Eigen::Map<const Eigen::VectorXd> rctvec_e((rct_vec+(maxAS-Toti)),Tot_p1);
  Eigen::Map<const Eigen::VectorXd> lgrct_e((lgrct+(maxAS-Toti)),Tot_p1);
  Eigen::Map<Eigen::VectorXd> pvec_e(pvec,Tot_p1);
  Eigen::Map<Eigen::VectorXd> tmpctvec(tmpctvec_,Tot_p1);
  Eigen::Map<Eigen::VectorXd> tmprctvec(tmprctvec_,Tot_p1);
  
  // Probability Vector //
  tmpctvec = ctvec_e.array()+piv;
  tmprctvec = rctvec_e.array()+ai_piv;
  pvec_e = -lgct_e-lgrct_e;
  pvec_e += tmpctvec.unaryExpr(&ArrLgamma);
  pvec_e += tmprctvec.unaryExpr(&ArrLgamma);
  pvec_e = pvec_e.array()+R::lgammafn(Tot_p1)+R::lgammafn(vpsi)-R::lgammafn(piv)-R::lgammafn(ai_piv)-R::lgammafn(vpsi+Toti);
  pvec_e = pvec_e.array().exp();
  
  // Expected Values:
  outvec[0] = (pvec_e.array()*(tmpctvec.unaryExpr(&ArrDigamma)).array()).sum();
  outvec[1] = (pvec_e.array()*(tmprctvec.unaryExpr(&ArrDigamma)).array()).sum();
  outvec[2] = (pvec_e.array()*(tmpctvec.unaryExpr(&ArrTrigamma)).array()).sum();
  outvec[3] = (pvec_e.array()*(tmprctvec.unaryExpr(&ArrTrigamma)).array()).sum();
}

/**********************************************************/
// SCORE TEST FUNCTION                                    //
/**********************************************************/
/*
NAME: CisTrans_Score
DESC: Computes a score test to determine whether there is a 
differential eQTL effect found with TReC and ASE components
of the likelihood. Essentially, determines whether or not
we can trust the result of a TReCASE model or whether we 
should revert to just a TReC Model. 

INPUT DIAGRAM:

-------------------------------------------------------------
|VARIABLE          || TYPE      || DESCRIPTION              |
-------------------------------------------------------------
| score_val        || double*   || Address of cisTrans score|
||           || test statistic.          |
_____________________________________________________________
| pval             || double*   || Address of cisTrans score|
|                  ||           || test p-value.            |
_____________________________________________________________
| ScoreVec_        || double*   || Address of cisTrans score|
|		   ||           || vector for the alpha     |
|                  ||           || components.              |
_____________________________________________________________
| mu_              || double*   || vector 0:n containing    |
|                  ||           || computed means for each  |
|                  ||           || subject.                 |
_____________________________________________________________
| musq_            || double*   || vector 0:n which will    |
|                  ||           || contain the squared mean |
|                  ||           || for each subject.        |
_____________________________________________________________
| Delta1_          || double*   || Matrix Placeholder       |
_____________________________________________________________
| Delta2_          || double*   || Matrix Placeholder       |
_____________________________________________________________
| Delta3_          || double*   || Matrix Placeholder       |
_____________________________________________________________
| Ibb_             || double*   || Matrix Placeholder       |
_____________________________________________________________
| Ibe_             || double*   || Matrix Placeholder       |
_____________________________________________________________
| Iee_             || double*   || Matrix Placeholder       |
_____________________________________________________________
| Iep_             || double*   || Matrix Placeholder       |
_____________________________________________________________
| Iea_             || double*   || Matrix Placeholder       |
_____________________________________________________________
| Ipp_             || double*   || Matrix Placeholder       |
_____________________________________________________________
| Ipa_             || double*   || Matrix Placeholder       |
_____________________________________________________________
| Iaa_             || double*   || Matrix Placeholder       |
_____________________________________________________________
| M1_              || double*   || Matrix Placeholder       |
_____________________________________________________________
| M2_              || double*   || Matrix Placeholder       |
_____________________________________________________________
| ex1              || double*   || Extra Parameters (TReC)  |
_____________________________________________________________
| ex2              || double*   || Extra Parameters (ASE)   |
_____________________________________________________________
| rhosAS           || double*   || Tumor Purities for the AS|
|                  ||           || subjects.                |
_____________________________________________________________
| M                || int       || Number of covariates.    |
_____________________________________________________________
| n                || int       || Number of subjects.      |
_____________________________________________________________
| ctvec            || double*   || vector 0:maxAS           |
_____________________________________________________________
| lgct             || double*   || log(ctvec!)              |
_____________________________________________________________
| rct_vec          || double*   || vector maxAS:0           |
_____________________________________________________________
| lgrct            || double*   || log(rct_vec!)            |
_____________________________________________________________
| pvec             || double*   || Placeholder vector for   |
|                  ||           || probabilities            |
_____________________________________________________________
| outvec           || double*   || Address of output storage|
|                  ||           || for a single subject     |
_____________________________________________________________
| Kappa            || double    || Estimated Kappa Value    |
_____________________________________________________________
| Eta              || double    || Estimated Eta Value      |
_____________________________________________________________
| Gamma            || double    || Estimated Gamma value    |
_____________________________________________________________
| Phi              || double    || Estimated Phi Value      |
_____________________________________________________________
| Psi              || double    || Estimated Psi Value      |
_____________________________________________________________


OUTPUT DIAGRAM: 
(VOID Function) score_val and pval serve as the returns.

outScore = The score test statistic
outSPval = The score test p-value using a chi-squared distribution
with 2 degrees of freedom

This is a single row of an output matrix which contains these
four values for all subjects in the ASE expression dataset. 

CODE:
*/

void CisTrans_Score(double* score_val, double* pval, double* ScoreVec_,
                    double* mu_, double* musq_,double* Delta1_,double* Delta2_,double* Delta3_,
                    double* Ibb_, double* Ibe_, double* Iee_, double* Iep_, double* Iea_,
                    double* Ipp_, double* Ipa_, double* Iaa_, double* M1_, double* M2_,  
                    double* ex1, double* ex2, double* rhosAS, int M, int n, double* Xmat_,
                    double* ctvec, double* lgct, double* rct_vec, double* lgrct, double* pvec, double* tmpctvec, double* tmprctvec,
                    double* Dmu_, double* deps_, double* offset_, double* Beta_,
                    double* outvec, double psi, double phi, double kappa, double eta, double gamma, int maxAS, double* ZeroMat_F_, double* OIMat_){
  // Eigen: Mean and Variance Pieces //
  Eigen::Map<Eigen::VectorXd> mu(mu_,n);
  Eigen::Map<Eigen::VectorXd> musq(musq_,n);
  Eigen::VectorXd vari(n);
  
  // Eigen: X Matrix Input Pieces //
  Eigen::Map<Eigen::MatrixXd> Xmat(Xmat_,n,M);
  Eigen::Map<Eigen::VectorXd> Beta(Beta_,M);
  Eigen::Map<Eigen::MatrixXd> Delta1(Delta1_,n,n);
  Eigen::Map<Eigen::MatrixXd> Delta2(Delta2_,n,n);
  Eigen::Map<Eigen::MatrixXd> Delta3(Delta3_,n,n);
  
  // Eigen: Information Matrix Pieces //
  Eigen::Map<Eigen::MatrixXd> Ibb(Ibb_,M,M);
  Eigen::Map<Eigen::MatrixXd> Ibe(Ibe_,M,3);
  Eigen::Map<Eigen::MatrixXd> Iee(Iee_,3,3);
  Eigen::Map<Eigen::MatrixXd> Iep(Iep_,3,1);
  Eigen::Map<Eigen::MatrixXd> Iea(Iea_,3,2);
  Eigen::Map<Eigen::MatrixXd> Ipp(Ipp_,1,1);
  Eigen::Map<Eigen::MatrixXd> Ipa(Ipa_,1,2);
  Eigen::Map<Eigen::MatrixXd> Iaa(Iaa_,2,2);
  Eigen::Map<Eigen::MatrixXd> M1(M1_,(M+3+1),(M+3+1));
  Eigen::Map<Eigen::MatrixXd> M2(M2_,(M+3+1),2);
  Eigen::Map<Eigen::Matrix<double,2,1> > ScoreVec(ScoreVec_);
  
  Eigen::Map<Eigen::MatrixXd> Dmu(Dmu_,n,3);
  Eigen::Map<Eigen::MatrixXd> deps(deps_,n,3);
  
  Eigen::Map<Eigen::MatrixXd> ZeroMat_F(ZeroMat_F_,n,n);
  Eigen::Map<Eigen::MatrixXd> OIMat(OIMat_,(M+6),(M+6));
    
    Ibb.block(0,0,M,M) = ZeroMat_F.block(0,0,M,M);
    Iee.block(0,0,3,3) = ZeroMat_F.block(0,0,3,3);
    Ipp.block(0,0,1,1) = ZeroMat_F.block(0,0,1,1);
    Ipa.block(0,0,1,2) = ZeroMat_F.block(0,0,1,2);
    Iaa.block(0,0,2,2) = ZeroMat_F.block(0,0,2,2);
    Iea.block(0,0,3,2) = ZeroMat_F.block(0,0,3,2);
    Iep.block(0,0,3,1) = ZeroMat_F.block(0,0,3,1);
    M1.block(0,0,(M+5),(M+5)) = ZeroMat_F.block(0,0,(M+5),(M+5));
    M2.block(0,0,(M+5),2) = ZeroMat_F.block(0,0,(M+5),2);
    ScoreVec.block(0,0,2,1) = ZeroMat_F.block(0,0,2,1);
    Delta1.block(0,0,n,n) = ZeroMat_F.block(0,0,n,n);
    Delta2.block(0,0,n,n) = ZeroMat_F.block(0,0,n,n);
    Delta3.block(0,0,n,n) = ZeroMat_F.block(0,0,n,n);
    OIMat.block(0,0,(M+7),(M+7)) = ZeroMat_F.block(0,0,(M+6),(M+6));
    
  
  //TReC: Incidentals //
  double pi, nu0, delta, ci, dci_dkappa, epsi;
  int n0   = ceil(ex1[2]-0.5),
      n1   = ceil(ex1[3]-0.5),
      i;
  int np9  = ceil(ex1[0]-0.5)+9;
  int dnp9 = np9+ceil(ex1[0]-0.5); 
  
  // Compute the means and variances //
  compute_nu0(Xmat_,n,M,Beta_,ex1);
  comp_offset(kappa,eta,gamma,ex1,offset_);
  
  Eigen::Map<Eigen::VectorXd> OffVec(offset_,n);
  
  mu = (Xmat*Beta+OffVec).array().exp();
  
  musq = mu.array()*mu.array();
  vari = mu+phi*musq;
  
  // Check 1: What do the means, squared means, and variances look like?
  /*Rprintf("Cis-Trans Score Test: Checking");
  Rprintf("Mean: %.10f,%.10f,%.10f,%.10f\n",mu(0),mu(1),mu(2),mu(3));
  Rprintf("MnSq: %.10f,%.10f,%.10f,%.10f\n",musq(0),musq(1),musq(2),musq(3));
  Rprintf("Vari: %.10f,%.10f,%.10f,%.10f\n",vari(0),vari(1),vari(2),vari(3));
  
  std::cout << "Mean Vec: " << mu.head(20) << std::endl;
  std::cout << "Mean Sqd: " << musq.head(20) << std::endl;
  std::cout << "Variance: " << vari.head(20) << std::endl;*/
  
  for(i=0;i<n0;i++){
    pi  = ex1[np9+i];
    nu0 = ex1[dnp9+i];
    delta = 1.0 - pi + pi*kappa;
    
    ci = (pi*kappa)/delta;
    dci_dkappa = ci*(1.0-ci)*(1.0/kappa);
    
    // Setting up Deps
    epsi      = (1.0-ci)*eta+ci*gamma;
    deps(i,0) = (gamma-eta)*dci_dkappa;
    deps(i,1) = 1.0-ci;
    deps(i,2) = ci;
    
    // Setting up the mean derivatives
    Dmu(i,0) = nu0*pi;
    Dmu(i,1) = 0;
    Dmu(i,2) = 0;
    
  }
  
  for(i=n0;i<(n0+n1);i++){
    pi  = ex1[np9+i];
    nu0 = ex1[dnp9+i];
    delta = 1.0 - pi + pi*kappa;
    
    ci = (pi*kappa)/delta;
    dci_dkappa = ci*(1.0-ci)*(1.0/kappa);
    
    // Setting up Deps
    epsi      = (1.0-ci)*eta+ci*gamma;
    deps(i,0) = (gamma-eta)*dci_dkappa;
    deps(i,1) = 1.0-ci;
    deps(i,2) = ci;
    
    // Setting up the mean derivatives
    Dmu(i,0) = nu0*(pi*0.5*(1.0+epsi)+delta*0.5*deps(i,0));
    Dmu(i,1) = nu0*delta*0.5*deps(i,1);
    Dmu(i,2) = nu0*delta*0.5*deps(i,2);
  }
  
  for(i=(n0+n1);i<n;i++){
    pi  = ex1[np9+i];
    nu0 = ex1[dnp9+i];
    delta = 1.0 - pi + pi*kappa;
    
    ci = (pi*kappa)/delta;
    dci_dkappa = ci*(1.0-ci)*(1.0/kappa);
    
    // Setting up Deps
    epsi      = (1.0-ci)*eta+ci*gamma;
    deps(i,0) = (gamma-eta)*dci_dkappa;
    deps(i,1) = 1.0-ci;
    deps(i,2) = ci;
    
    // Setting up the mean derivatives
    Dmu(i,0) = nu0*(pi*epsi+delta*deps(i,0));
    Dmu(i,1) = nu0*delta*deps(i,1);
    Dmu(i,2) = nu0*delta*deps(i,2);
  }
  
  // TREC: Components //
  Delta1.diagonal() = (1.0/vari.array());
  Delta2.diagonal() = mu.array()*Delta1.diagonal().array();
  Delta3.diagonal() = musq.array()*Delta1.diagonal().array();
  
  Ibb = Xmat.transpose()*Delta3*Xmat;
  Ibe = Xmat.transpose()*Delta2*Dmu;
  Iee = Dmu.transpose()*Delta1*Dmu;
  
  // Checking Values for ease of use //
  /*std::cout << "Ibb matrix: \n" << Ibb << std::endl;
  std::cout << "Ibe matrix: \n" << Ibe << std::endl;
  std::cout << "Iee matrix (Stage 1): " << Iee << std::endl;*/
  
  // ASE: Expectations and Incidentals //
  double vpsi = 1.0/psi,
    piv,
    ai_pi,
    ai_piv,
    dgpv,
    dgap,
    dgvp = R::digamma(vpsi),
    trpv,
    trap,
    trvp=R::trigamma(vpsi),
    vepsi;
  double d2ci_dkappa2, dpi_dkappa, dpi_deta, dpi_dgamma, dpi_dae, dpi_dag,
  dpi_dkappa2, dpi_dkappa_deta, dpi_dkappa_dgamma, dpi_dkappa_dae, dpi_dkappa_dag,
  dpi_deta2, dpi_deta_dgamma, dpi_deta_dae, dpi_deta_dag,
  dpi_dgamma2, dpi_dgamma_dae, dpi_dgamma_dag,
  dpi_dae2, dpi_dae_dag,
  dpi_dag2,
  depsi_de, depsi_dg, depsi_dkappa_dg, depsi_dkappa_de;
  double B, dB_dp, dB_dpi, B_score;
  int nAS = ceil(ex2[0]-0.5),
      nAS_0 = ceil(ex2[1]-0.5),
      nAS_1 = ceil(ex2[2]-0.5),
      nASp4 = nAS+4,
      nTi,
      nBi;
  
  for(i=0;i<nAS_0;i++){
    // Finding Expectations
    nTi = ceil(ex2[4+i]-0.5);
    
    ASE_ExpFunc(ctvec,lgct,rct_vec,lgrct,pvec,0.5,nTi,vpsi,outvec,tmpctvec,tmprctvec,maxAS);
    
    //Plugging In:
    Ipp(0) -= (2)*vpsi*vpsi*vpsi*(dgvp-R::digamma((0.5*vpsi))-
      R::digamma(vpsi+nTi)+0.5*(outvec[0]+outvec[1]));
    Ipp(0) -= vpsi*vpsi*vpsi*vpsi*(trvp-0.5*R::trigamma(0.5*vpsi)-
      R::trigamma(vpsi+nTi)+0.25*(outvec[2]+outvec[3]));
    
  }
  
  for(i=nAS_0;i<(nAS_0+nAS_1);i++){
    // Common Definitions and Derivatives
    nTi = ceil(ex2[4+i]-0.5);
    nBi = ceil(ex2[nASp4+i]-0.5);
    
    delta        = 1-rhosAS[i]+rhosAS[i]*kappa;
    ci           = (rhosAS[i]*kappa)/delta;
    dci_dkappa   = ci*(1-ci)*(1.0/kappa);
    d2ci_dkappa2 = (1.0-2*ci)*(1.0/kappa)*dci_dkappa-ci*(1.0-ci)*(1.0/(kappa*kappa));
    epsi         = (1.0-ci)*eta+ci*gamma; 
    vepsi        = 1.0/(1.0+epsi);
    
    pi = epsi*vepsi;
    ai_pi = 1.0 - pi;
    piv = pi*vpsi;
    dgpv = R::digamma(piv);
    trpv = R::trigamma(piv);
    
    ai_piv = ai_pi*vpsi;
    dgap = R::digamma(ai_piv);
    trap = R::trigamma(ai_piv);
    
    //dgvp = R::digamma(vpsi);
    //trvp = R::trigamma(vpsi);
    //dgpn = R::digamma(vpsi+nTi);
    //trpn = R::trigamma(vpsi+nTi);
    
    depsi_dkappa_de = -dci_dkappa;
    depsi_dkappa_dg = dci_dkappa;
    depsi_de = (1.0-ci);
    depsi_dg = ci;
    
    dpi_dkappa        = vepsi*vepsi*(gamma-eta)*dci_dkappa;
    dpi_dkappa2       = vepsi*vepsi*(gamma-eta)*(d2ci_dkappa2-2*vepsi*(gamma-eta)*dci_dkappa*dci_dkappa);
    dpi_dkappa_deta   = vepsi*vepsi*(depsi_dkappa_de-2*vepsi*(gamma-eta)*dci_dkappa*depsi_de);
    dpi_dkappa_dgamma = vepsi*vepsi*(depsi_dkappa_dg-2*vepsi*(gamma-eta)*dci_dkappa*depsi_dg);
    dpi_dkappa_dae    = vepsi*vepsi*(depsi_dkappa_de-2*vepsi*(gamma-eta)*dci_dkappa*depsi_de);
    dpi_dkappa_dag    = vepsi*vepsi*(depsi_dkappa_dg-2*vepsi*(gamma-eta)*dci_dkappa*depsi_dg);
    
    dpi_deta        = vepsi*vepsi*depsi_de;
    dpi_deta2       = -2*vepsi*vepsi*vepsi*depsi_de*depsi_de;
    dpi_deta_dgamma = -2*vepsi*vepsi*vepsi*depsi_de*depsi_dg;
    dpi_deta_dae    = -2*vepsi*vepsi*vepsi*depsi_de*depsi_de;
    dpi_deta_dag    = -2*vepsi*vepsi*vepsi*depsi_de*depsi_dg;
    
    dpi_dgamma     = vepsi*vepsi*depsi_dg;
    dpi_dgamma2    = -2*vepsi*vepsi*vepsi*depsi_dg*depsi_dg;
    dpi_dgamma_dae = -2*vepsi*vepsi*vepsi*depsi_dg*depsi_de;
    dpi_dgamma_dag = -2*vepsi*vepsi*vepsi*depsi_dg*depsi_dg;
    
    dpi_dae     = dpi_deta;
    dpi_dae2    = -2*vepsi*vepsi*vepsi*depsi_de*depsi_de;
    dpi_dae_dag = -2*vepsi*vepsi*vepsi*depsi_de*depsi_dg;
    
    dpi_dag  = dpi_dgamma;
    dpi_dag2 = -2*vepsi*vepsi*vepsi*depsi_dg*depsi_dg;
    
    // Finding Expectations
    
    ASE_ExpFunc(ctvec,lgct,rct_vec,lgrct,pvec,pi,nTi,vpsi,outvec,tmpctvec,tmprctvec,maxAS);
    
    //Plugging In:
    B = outvec[0]-outvec[1]-dgpv+dgap;
    B_score = R::digamma(vpsi*pi+nBi)-R::digamma(vpsi*ai_pi+nTi-nBi)-dgpv+dgap;
    dB_dp = pi*outvec[2]-ai_pi*outvec[3]-trpv*pi+ai_pi*trap;
    dB_dpi = outvec[2]+outvec[3]-trpv-trap;
    
    Iep(0,0) -= -vpsi*vpsi*(B+vpsi*dB_dp)*dpi_dkappa;                      // kappa
    Iep(1,0) -= -vpsi*vpsi*(B+vpsi*dB_dp)*dpi_deta;	                       // eta
    Iep(2,0) -= -vpsi*vpsi*(B+vpsi*dB_dp)*dpi_dgamma;                      // gamma
    Ipa(0,0) -= -vpsi*vpsi*(B+vpsi*dB_dp)*dpi_dae;                         // alpha_eta
    Ipa(0,1) -= -vpsi*vpsi*(B+vpsi*dB_dp)*dpi_dag;                         // alpha_gamma
    Iea(0,0) -= vpsi*B*dpi_dkappa_dae+vpsi*vpsi*dB_dpi*dpi_dkappa*dpi_dae; // kappa alpha_eta
    Iea(0,1) -= vpsi*B*dpi_dkappa_dag+vpsi*vpsi*dB_dpi*dpi_dkappa*dpi_dag; // kappa alpha_gamma
    Iea(1,0) -= vpsi*B*dpi_deta_dae+vpsi*vpsi*dB_dpi*dpi_deta*dpi_dae;     // eta alpha_eta
    Iea(1,1) -= vpsi*B*dpi_deta_dag+vpsi*vpsi*dB_dpi*dpi_deta*dpi_dag;     // eta alpha_gamma
    Iea(2,0) -= vpsi*B*dpi_dgamma_dae+vpsi*vpsi*dB_dpi*dpi_dgamma*dpi_dae; // gamma alpha_eta
    Iea(2,1) -= vpsi*B*dpi_dgamma_dag+vpsi*vpsi*dB_dpi*dpi_dgamma*dpi_dag; // gamma alpha_gamma
    
    Iee(0,0) -= vpsi*B*dpi_dkappa2 + vpsi*vpsi*dB_dpi*dpi_dkappa*dpi_dkappa;       // kappa kappa
    Iee(0,1) -= vpsi*B*dpi_dkappa_deta + vpsi*vpsi*dB_dpi*dpi_dkappa*dpi_deta;     // kappa eta
    Iee(0,2) -= vpsi*B*dpi_dkappa_dgamma + vpsi*vpsi*dB_dpi*dpi_dkappa*dpi_dgamma; // kappa gamma
    Iee(1,1) -= vpsi*B*dpi_deta2 + vpsi*vpsi*dB_dpi*dpi_deta*dpi_deta;             // eta eta
    Iee(1,2) -= vpsi*B*dpi_deta_dgamma + vpsi*vpsi*dB_dpi*dpi_deta*dpi_dgamma;     // eta gamma
    Iee(2,2) -= vpsi*B*dpi_dgamma2 + vpsi*vpsi*dB_dpi*dpi_dgamma*dpi_dgamma;       // gamma gamma		
    
    
    Iaa(0,0) -= vpsi*B*dpi_dae2 + vpsi*vpsi*dB_dpi*dpi_dae*dpi_dae;    // alpha_eta alpha_eta
    Iaa(0,1) -= vpsi*B*dpi_dae_dag + vpsi*vpsi*dB_dpi*dpi_dae*dpi_dag; // alpha_eta alpha_gamma
    Iaa(1,1) -= vpsi*B*dpi_dag2 + vpsi*vpsi*dB_dpi*dpi_dag*dpi_dag;    // alpha_gamma alpha_gamma
    
    Ipp(0) -= (2)*vpsi*vpsi*vpsi*(R::digamma(vpsi)+pi*outvec[0]+ai_pi*outvec[1]-
      pi*R::digamma(piv)-ai_pi*R::digamma(ai_piv)-R::digamma(vpsi+nTi));
    Ipp(0) -= vpsi*vpsi*vpsi*vpsi*(R::trigamma(vpsi)+pi*pi*outvec[2]+ai_pi*ai_pi*outvec[3]-
      pi*pi*R::trigamma(piv)-ai_pi*ai_pi*R::trigamma(ai_piv)-
      R::trigamma(vpsi+nTi));
    
    // Score Function Contribution:
    ScoreVec(0,0) += vpsi*B_score*dpi_dae;
    ScoreVec(1,0) += vpsi*B_score*dpi_dag;
    
    // Checking where the NaN is introducted ( only want to call attention to the subject and their values):
    /*std::cout << "Subject " << i << "\n Pertinent Values: nTi = " << nTi << "// Rho = " << rhosAS[i] << std::endl;
      std::cout << "Kappa = " << kappa << " // Eta = " << eta << " // gamma = " << gamma << "\n B = " << B << " // dpi_dae = " << dpi_dae << std::endl;
      std::cout << "OutVec[0//1//2//3] = " << outvec[0] << "//" << outvec[1] << "//" << outvec[2] << "//" << outvec[3] << std::endl;*/
  }
  
  for(i=(nAS_0+nAS_1);i<nAS;i++){
    // Finding Expectations
    nTi = ceil(ex2[4+i]-0.5);
    
    ASE_ExpFunc(ctvec,lgct,rct_vec,lgrct,pvec,0.5,nTi,vpsi,outvec,tmpctvec,tmprctvec,maxAS);
    
    //Plugging In:
    Ipp(0) -= (2)*vpsi*vpsi*vpsi*(R::digamma(vpsi)-R::digamma((0.5*vpsi))-
      R::digamma(vpsi+nTi)+0.5*(outvec[0]+outvec[1]));
    Ipp(0) -= vpsi*vpsi*vpsi*vpsi*(R::trigamma(vpsi)-0.5*R::trigamma(0.5*vpsi)-
      R::trigamma(vpsi+nTi)+0.25*(outvec[2]+outvec[3]));
  }
  
  Iaa(1,0) = Iaa(0,1);
  Iee(1,0) = Iee(0,1);
  Iee(2,0) = Iee(0,2);
  Iee(2,1) = Iee(1,2);
  
  // Reporting the Values of matrices so that I can utilize them for checking
  /*std::cout << "Iee matrix: \n" << Iee << std::endl;
  std::cout << "Iep matrix: \n" << Iep << std::endl;
  std::cout << "Iea matrix: \n" << Iea << std::endl;
  std::cout << "Ipp matrix: \n" << Ipp << std::endl;
  std::cout << "Ipa matrix: \n" << Ipa << std::endl;
  std::cout << "Iaa matrix: \n" << Iaa << std::endl;
  std::cout << "Score Vect: \n" << ScoreVec << std::endl;*/
  
  // Reporting Components to see if M1 is properly specified:
  //Debugging purposes only!
//  std::ofstream out_data;
//  out_data.open("Test_out.txt");
//  out_data << "Ibb matrix: \n" << Ibb << "\n";
//  out_data << "Ibe matrix: \n" << Ibe << "\n";
//  out_data << "Iee matrix: \n" << Iee << "\n";
//  out_data << "Iep matrix: \n" << Iep << "\n";
//  out_data << "Ipp value : \n" << Ipp << "\n";
//  out_data << "Iea Matrix: \n" << Iea << "\n";
//  out_data << "Ipa matrix: \n" << Ipa << "\n";
//  out_data << "Iaa matrix: \n" << Iaa << "\n";
  
  // Final Information Matrix Components//
  M1.block(0,0,M,M)     = Ibb;
  M1.block(0,M,M,3)     = Ibe;
  M1.block(M,0,3,M)     = Ibe.transpose();
  M1.block<3,3>(M,M)     = Iee;
  M1.block<3,1>(M,(M+3)) = Iep;
  M1.block<1,3>((M+3),M) = Iep.transpose();
  M1((M+3),(M+3))        = Ipp(0);
  
  M2.block<3,2>(M,0)     = Iea;
  M2.block<1,2>((M+3),0) = Ipa;
  
  // Checking the formation of M1 and M2:
  // Only for Debugging Purposes!
//  out_data << "M1 matrix: \n" << M1 << "\n";
//  out_data << "M2 matrix: \n" << M2 << "\n";
//  out_data << "Iaa-Blah : \n" << Iaa-M2.transpose()*M1.inverse()*M2 << "\n";
//  out_data.close();
  
  OIMat.block(0,0,(M+4),(M+4)) = M1;
  OIMat.block(0,(M+4),(M+4),2) = M2;
  OIMat.block((M+4),0,2,(M+4)) = M2.transpose();
  OIMat.block((M+4),(M+4),2,2) = Iaa;
  
  // Checking Invertibility
  Eigen::JacobiSVD<Eigen::MatrixXd> svd1(M1);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd2(Iaa-M2.transpose()*M1.inverse()*M2);
  double cond1 = svd1.singularValues()(0)/svd1.singularValues()(svd1.singularValues().size()-1);
  double cond2 = svd2.singularValues()(0)/svd2.singularValues()(svd2.singularValues().size()-1);
  
  //Rprintf("Checking Condition Numbers: Cond1- %.10f // Cond2- %.10f\n",cond1,cond2);
  
  // Checking Positive Definiteness:
  Eigen::LLT<Eigen::MatrixXd> lltOfA(OIMat);
  
  if(lltOfA.info()==Eigen::NumericalIssue){
    *score_val = -6.0;
    *pval = 1.0;
  } else if((cond1>(1/DBL_EPSILON))||(cond2>(1/DBL_EPSILON))){
    *score_val = -5.0;
    *pval = 1.0;
  } else {
    *score_val = ScoreVec.transpose()*((Iaa-M2.transpose()*M1.inverse()*M2).inverse())*ScoreVec;
    *pval = R::pchisq(*score_val,2,0,0);
  }
  
  if(std::isnan(*score_val)){
    *score_val = -6.0;
    *pval = 1.0;
  }
  // This way if the score is -5.0, we know that we have invertibility issues.
}


void CisTrans_ObsScore(double* score_val, double* pval, double* ScoreVec_,
                       double* mu_, double* musq_,double* Delta1_,double* Delta2_,double* Delta3_,
                       double* Ibb_, double* Ibe_, double* Iee_, double* Iep_, double* Iea_,
                       double* Ipp_, double* Ipa_, double* Iaa_, double* M1_, double* M2_,  
                       double* ex1, double* ex2, double* rhosAS, int M, int n, double* Xmat_,
                       double* ctvec, double* lgct, double* rct_vec, double* lgrct, double* pvec, double* tmpctvec, double* tmprctvec,
                       double* Dmu_, double* deps_, double* offset_, double* Beta_,
                       double* outvec, double psi, double phi, double kappa, double eta, double gamma, int maxAS,
                       double* Ibt_, double* Iet_, double* Itt_,double* ZeroMat_F_, double* OIMat_){
  // Eigen: Mean and Variance Pieces //
  Eigen::Map<Eigen::VectorXd> Yv((ex1+9),n);
  Eigen::Map<Eigen::VectorXd> mu(mu_,n);
  Eigen::Map<Eigen::VectorXd> musq(musq_,n);
  Eigen::VectorXd vari(n);
  //Eigen::MatrixXd inverseMat(2,2);
  
  // Eigen: X Matrix Input Pieces //
  Eigen::Map<Eigen::MatrixXd> Xmat(Xmat_,n,M);
  Eigen::Map<Eigen::VectorXd> Beta(Beta_,M);
  Eigen::Map<Eigen::MatrixXd> Delta1(Delta1_,n,n);
  Eigen::Map<Eigen::MatrixXd> Delta2(Delta2_,n,n);
  Eigen::Map<Eigen::MatrixXd> Delta3(Delta3_,n,n);
  
  // Eigen: Information Matrix Pieces //
  Eigen::Map<Eigen::MatrixXd> Ibb(Ibb_,M,M);
  Eigen::Map<Eigen::MatrixXd> Ibe(Ibe_,M,3);
  Eigen::Map<Eigen::MatrixXd> Iee(Iee_,3,3);
  Eigen::Map<Eigen::MatrixXd> Ibt(Ibt_,M,1);
  Eigen::Map<Eigen::MatrixXd> Iet(Iet_,3,1);
  Eigen::Map<Eigen::MatrixXd> Itt(Itt_,1,1);
  Eigen::Map<Eigen::MatrixXd> Iep(Iep_,3,1);
  Eigen::Map<Eigen::MatrixXd> Iea(Iea_,3,2);
  Eigen::Map<Eigen::MatrixXd> Ipp(Ipp_,1,1);
  Eigen::Map<Eigen::MatrixXd> Ipa(Ipa_,1,2);
  Eigen::Map<Eigen::MatrixXd> Iaa(Iaa_,2,2);
  Eigen::Map<Eigen::MatrixXd> M1(M1_,(M+5),(M+5));
  Eigen::Map<Eigen::MatrixXd> M2(M2_,(M+5),2);
  Eigen::Map<Eigen::Matrix<double,2,1> > ScoreVec(ScoreVec_);
  Eigen::Map<Eigen::MatrixXd> OIMat(OIMat_,(M+7),(M+7));
  
  Eigen::Map<Eigen::MatrixXd> Dmu(Dmu_,n,3);
  Eigen::Map<Eigen::MatrixXd> deps(deps_,n,3);
  
  Eigen::Map<Eigen::MatrixXd> ZeroMat_F(ZeroMat_F_,n,n);
  
  // This initialization works for any time the test is used after the first time
  // Need to ensure that the original matrices from which these arise are set to 
  // 0 for the first implementation.
  Ibb.block(0,0,M,M) = ZeroMat_F.block(0,0,M,M);
  Iee.block(0,0,3,3) = ZeroMat_F.block(0,0,3,3);
  Itt.block(0,0,1,1) = ZeroMat_F.block(0,0,1,1);
  Ipp.block(0,0,1,1) = ZeroMat_F.block(0,0,1,1);
  Ipa.block(0,0,1,2) = ZeroMat_F.block(0,0,1,2);
  Iaa.block(0,0,2,2) = ZeroMat_F.block(0,0,2,2);
  Iea.block(0,0,3,2) = ZeroMat_F.block(0,0,3,2);
  Iep.block(0,0,3,1) = ZeroMat_F.block(0,0,3,1);
  M1.block(0,0,(M+5),(M+5)) = ZeroMat_F.block(0,0,(M+5),(M+5));
  M2.block(0,0,(M+5),2) = ZeroMat_F.block(0,0,(M+5),2);
  ScoreVec.block(0,0,2,1) = ZeroMat_F.block(0,0,2,1);
  Delta1.block(0,0,n,n) = ZeroMat_F.block(0,0,n,n);
  Delta2.block(0,0,n,n) = ZeroMat_F.block(0,0,n,n);
  Delta3.block(0,0,n,n) = ZeroMat_F.block(0,0,n,n);
  OIMat.block(0,0,(M+7),(M+7)) = ZeroMat_F.block(0,0,(M+7),(M+7));
  
  //TReC: Incidentals //
  double pi, nu0, delta, ci, dci_dkappa, epsi, C_comp,
         vphi = 1.0/phi,
         dci_dkappa2,
         depsi_dkappa2,
         mfact;
  int n0   = ceil(ex1[2]-0.5),
    n1   = ceil(ex1[3]-0.5),
    i;
  int np9  = ceil(ex1[0]-0.5)+9;
  int dnp9 = np9+ceil(ex1[0]-0.5); 
  
  // Compute the means and variances //
  compute_nu0(Xmat_,n,M,Beta_,ex1);
  comp_offset(kappa,eta,gamma,ex1,offset_);
  
  Eigen::Map<Eigen::VectorXd> OffVec(offset_,n);
  
  mu = (Xmat*Beta+OffVec).array().exp();
  
  musq = mu.array()*mu.array();
  vari = mu+phi*musq;
  
  // Check 1: What do the means, squared means, and variances look like?
  /*Rprintf("Cis-Trans Score Test: Checking");
  Rprintf("Mean: %.10f,%.10f,%.10f,%.10f\n",mu(0),mu(1),mu(2),mu(3));
  Rprintf("MnSq: %.10f,%.10f,%.10f,%.10f\n",musq(0),musq(1),musq(2),musq(3));
  Rprintf("Vari: %.10f,%.10f,%.10f,%.10f\n",vari(0),vari(1),vari(2),vari(3));
  
  std::cout << "Mean Vec: " << mu.head(20) << std::endl;
  std::cout << "Mean Sqd: " << musq.head(20) << std::endl;
  std::cout << "Variance: " << vari.head(20) << std::endl;*/
  
  for(i=0;i<n0;i++){
    pi  = ex1[np9+i];
    nu0 = ex1[dnp9+i];
    delta = 1.0 - pi + pi*kappa;
    
    ci = (pi*kappa)/delta;
    dci_dkappa = ci*(1.0-ci);
    
    // Setting up Deps
    epsi      = (1.0-ci)*eta+ci*gamma;
    deps(i,0) = (gamma-eta)*dci_dkappa;
    deps(i,1) = (1.0-ci)*eta;
    deps(i,2) = (ci)*gamma;
    
    // Setting up the mean derivatives
    Dmu(i,0) = nu0*pi*kappa;
    Dmu(i,1) = 0;
    Dmu(i,2) = 0;
    
    // Adding to Iee
    mfact = ((Yv(i)-mu(i))/vari(i));
    Iee(0,0) -= mfact*nu0*(pi*kappa);
    
    // Adding to Itt
    C_comp = R::digamma(Yv(i)+vphi)-R::digamma(vphi)-log(1.0+phi*mu(i));
    Itt(0) -= (2*vphi*vphi*vphi*C_comp + vphi*vphi*vphi*vphi*(R::trigamma(Yv(i)+vphi)-R::trigamma(vphi))+
            2*vphi*vphi*(mu(i)/(1.0+phi*mu(i)))+(vphi+Yv(i))*(musq(i)*musq(i)/(vari(i)*vari(i)))-vphi*vphi*Yv(i));
  }
  
  for(i=n0;i<(n0+n1);i++){
    pi  = ex1[np9+i];
    nu0 = ex1[dnp9+i];
    delta = 1.0 - pi + pi*kappa;
    
    ci = (pi*kappa)/delta;
    dci_dkappa = ci*(1.0-ci);
    dci_dkappa2 = (1.0-2*ci)*dci_dkappa;
    depsi_dkappa2 = (gamma-eta)*dci_dkappa2;
    
    // Setting up Deps
    epsi      = (1.0-ci)*eta+ci*gamma;
    deps(i,0) = (gamma-eta)*dci_dkappa;
    deps(i,1) = (1.0-ci)*eta;
    deps(i,2) = ci*gamma;
    
    // Setting up the mean derivatives
    Dmu(i,0) = nu0*(pi*kappa*0.5*(1.0+epsi)+delta*0.5*deps(i,0));
    Dmu(i,1) = nu0*delta*0.5*deps(i,1);
    Dmu(i,2) = nu0*delta*0.5*deps(i,2);
    
    // Adding to Iee
    mfact = ((Yv(i)-mu(i))/vari(i));
    
    Iee(0,0) -= mfact*nu0*(pi*kappa*(1.0+epsi)*0.5+pi*kappa*deps(i,0)+delta*0.5*depsi_dkappa2);
    Iee(0,1) -= mfact*nu0*(pi*kappa*0.5*deps(i,1)-delta*0.5*dci_dkappa*eta);
    Iee(0,2) -= mfact*nu0*(pi*kappa*0.5*deps(i,2)+delta*0.5*dci_dkappa*gamma);
    Iee(1,1) -= mfact*nu0*delta*0.5*(1-ci)*eta;
    Iee(2,2) -= mfact*nu0*delta*0.5*ci*gamma;

    // Adding to Itt
    C_comp = R::digamma(Yv(i)+vphi)-R::digamma(vphi)-log(1.0+phi*mu(i));
    Itt(0) -= (2*vphi*vphi*vphi*C_comp + vphi*vphi*vphi*vphi*(R::trigamma(Yv(i)+vphi)-R::trigamma(vphi))+
      2*vphi*vphi*(mu(i)/(1.0+phi*mu(i)))+(vphi+Yv(i))*(musq(i)*musq(i)/(vari(i)*vari(i)))-vphi*vphi*Yv(i));
  }
  
  for(i=(n0+n1);i<n;i++){
    pi  = ex1[np9+i];
    nu0 = ex1[dnp9+i];
    delta = 1.0 - pi + pi*kappa;
    
    ci = (pi*kappa)/delta;
    dci_dkappa = ci*(1.0-ci);
    dci_dkappa2 = (1.0-2*ci)*dci_dkappa;
    depsi_dkappa2 = (gamma-eta)*dci_dkappa2;
    
    // Setting up Deps
    epsi      = (1.0-ci)*eta+ci*gamma;
    deps(i,0) = (gamma-eta)*dci_dkappa;
    deps(i,1) = (1.0-ci)*eta;
    deps(i,2) = ci*gamma;
    
    // Setting up the mean derivatives
    Dmu(i,0) = nu0*(pi*kappa*epsi+delta*deps(i,0));
    Dmu(i,1) = nu0*delta*deps(i,1);
    Dmu(i,2) = nu0*delta*deps(i,2);
    
    // Adding to Iee
    mfact = ((Yv(i)-mu(i))/vari(i));
    
    Iee(0,0) -= mfact*nu0*(pi*kappa*(epsi)+2*pi*kappa*deps(i,0)+delta*depsi_dkappa2);
    Iee(0,1) -= mfact*nu0*(pi*kappa*deps(i,1)-delta*dci_dkappa*eta);
    Iee(0,2) -= mfact*nu0*(pi*kappa*deps(i,2)+delta*dci_dkappa*gamma);
    Iee(1,1) -= mfact*nu0*delta*(1-ci)*eta;
    Iee(2,2) -= mfact*nu0*delta*ci*gamma;
    
    // Adding to Itt
    C_comp = R::digamma(Yv(i)+vphi)-R::digamma(vphi)-log(1.0+phi*mu(i));
    Itt(0) -= (2*vphi*vphi*vphi*C_comp + vphi*vphi*vphi*vphi*(R::trigamma(Yv(i)+vphi)-R::trigamma(vphi))+
      2*vphi*vphi*(mu(i)/(1.0+phi*mu(i)))+(vphi+Yv(i))*(musq(i)*musq(i)/(vari(i)*vari(i)))-vphi*vphi*Yv(i));
  }
  
  // Update the Iee components to match
  Iee(1,0) = Iee(0,1);
  Iee(2,0) = Iee(0,2);
  Iee(2,1) = Iee(1,2);
  
  // TREC: Components //
  Delta1.diagonal() = (1.0/vari.array());
  Delta2.diagonal() = mu.array()*Delta1.diagonal().array();
  Delta3.diagonal() = musq.array()*Delta1.diagonal().array();
  
  Ibb = Xmat.transpose()*Delta3*Xmat;
  Ibe = Xmat.transpose()*Delta2*Dmu;
  Iee += Dmu.transpose()*Delta1*Dmu;
  
  // Add the "Observed Pieces" for the other components
  Delta3.diagonal() = (Yv-mu).array()*mu.array()*musq.array()*Delta1.diagonal().array()*Delta1.diagonal().array();
  Delta2.diagonal() = (Yv-mu).array()*musq.array()*Delta1.diagonal().array()*Delta1.diagonal().array();
  Delta1.diagonal() = (Yv-mu).array()*(1.0+2.0*phi*mu.array()).array()/(vari.array()*vari.array());
  
  Ibb += phi*Xmat.transpose()*Delta3*Xmat;
  Ibe += phi*Xmat.transpose()*Delta2*Dmu;
  Iee += Dmu.transpose()*Delta1*Dmu;
  
  Ibt = Xmat.transpose()*Delta3.diagonal();
  Iet = Dmu.transpose()*Delta2.diagonal();
  
  //Testing again
  /*std::cout << "Iee matrix: \n" << Iee << std::endl;*/
  
  // Checking Values for ease of use //
  /*std::cout << "Ibb matrix: \n" << Ibb << std::endl;
  std::cout << "Ibe matrix: \n" << Ibe << std::endl;
  std::cout << "Iee matrix (Stage 1): " << Iee << std::endl;*/
  
  // ASE: Expectations and Incidentals //
  double vpsi = 1.0/psi,
    piv,
    ai_pi,
    ai_piv,
    dgpv,
    dgap,
    dgvp = R::digamma(vpsi),
    trpv,
    trap,
    trvp=R::trigamma(vpsi),
    vepsi;
  double d2ci_dkappa2, dpi_dkappa, dpi_deta, dpi_dgamma, dpi_dae, dpi_dag,
  dpi_dkappa2, dpi_dkappa_deta, dpi_dkappa_dgamma, dpi_dkappa_dae, dpi_dkappa_dag,
  dpi_deta2, dpi_deta_dgamma, dpi_deta_dae, dpi_deta_dag,
  dpi_dgamma2, dpi_dgamma_dae, dpi_dgamma_dag,
  dpi_dae2, dpi_dae_dag,
  dpi_dag2,
  depsi_de, depsi_dg, depsi_dkappa_dg, depsi_dkappa_de, depsi_dkappa_dag, depsi_dkappa_dae,
  depsi_dae, depsi_dag;
  double B, dB_dp, dB_dpi;
  int nAS = ceil(ex2[0]-0.5),
    nAS_0 = ceil(ex2[1]-0.5),
    nAS_1 = ceil(ex2[2]-0.5),
    nASp4 = nAS+4,
    nTi,
    nBi;
  
  for(i=0;i<nAS_0;i++){
    // Finding Expectations
    nTi = ceil(ex2[4+i]-0.5);
    nBi = ceil(ex2[nASp4+i]-0.5);
    
    //ASE_ExpFunc(ctvec,lgct,rct_vec,lgrct,pvec,0.5,nTi,vpsi,outvec,tmpctvec,tmprctvec,maxAS);
    outvec[0] = R::digamma(vpsi*0.5+nBi);
    outvec[1] = R::digamma(vpsi*0.5+nTi-nBi);
    outvec[2] = R::trigamma(vpsi*0.5+nBi);
    outvec[3] = R::trigamma(vpsi*0.5+nTi-nBi);
    
    //Plugging In:
    Ipp(0) -= vpsi*(dgvp-R::digamma((0.5*vpsi))-
      R::digamma(vpsi+nTi)+0.5*(outvec[0]+outvec[1]));
    Ipp(0) -= vpsi*vpsi*(trvp-0.5*R::trigamma(0.5*vpsi)-
      R::trigamma(vpsi+nTi)+0.25*(outvec[2]+outvec[3]));
    
  }
  
  for(i=nAS_0;i<(nAS_0+nAS_1);i++){
    // Common Definitions and Derivatives
    nTi = ceil(ex2[4+i]-0.5);
    nBi = ceil(ex2[nASp4+i]-0.5);
  
    delta        = 1-rhosAS[i]+rhosAS[i]*kappa;
    ci           = (rhosAS[i]*kappa)/delta;
    dci_dkappa   = ci*(1-ci);
    d2ci_dkappa2 = (1.0-2*ci)*dci_dkappa;
    epsi         = (1.0-ci)*eta+ci*gamma; 
    vepsi        = 1.0/(1.0+epsi);
    
    pi = epsi*vepsi;
    ai_pi = 1.0 - pi;
    piv = pi*vpsi;
    dgpv = R::digamma(piv);
    trpv = R::trigamma(piv);
    
    ai_piv = ai_pi*vpsi;
    dgap = R::digamma(ai_piv);
    trap = R::trigamma(ai_piv);
    
    //dgvp = R::digamma(vpsi);
    //trvp = R::trigamma(vpsi);
    //dgpn = R::digamma(vpsi+nTi);
    //trpn = R::trigamma(vpsi+nTi);
    
    depsi_dkappa_de  = -dci_dkappa*eta;
    depsi_dkappa_dae = -dci_dkappa;
    depsi_dkappa_dg  = dci_dkappa*gamma;
    depsi_dkappa_dag = dci_dkappa;
    depsi_de = (1.0-ci)*eta;
    depsi_dae = (1.0-ci);
    depsi_dg = ci*gamma;
    depsi_dag = ci;
    
    dpi_dkappa        = vepsi*vepsi*(gamma-eta)*dci_dkappa;
    dpi_dkappa2       = vepsi*vepsi*((gamma-eta)*d2ci_dkappa2-2*vepsi*(gamma-eta)*(gamma-eta)*dci_dkappa*dci_dkappa);
    dpi_dkappa_deta   = vepsi*vepsi*(depsi_dkappa_de-2*vepsi*(gamma-eta)*dci_dkappa*depsi_de);
    dpi_dkappa_dgamma = vepsi*vepsi*(depsi_dkappa_dg-2*vepsi*(gamma-eta)*dci_dkappa*depsi_dg);
    dpi_dkappa_dae    = vepsi*vepsi*(depsi_dkappa_dae-2*vepsi*(gamma-eta)*dci_dkappa*depsi_dae);
    dpi_dkappa_dag    = vepsi*vepsi*(depsi_dkappa_dag-2*vepsi*(gamma-eta)*dci_dkappa*depsi_dag);
    
    dpi_deta        = vepsi*vepsi*depsi_de;
    dpi_deta2       = -2*vepsi*vepsi*vepsi*depsi_de*depsi_de+vepsi*vepsi*depsi_de;
    dpi_deta_dgamma = -2*vepsi*vepsi*vepsi*depsi_de*depsi_dg;
    dpi_deta_dae    = -2*vepsi*vepsi*vepsi*depsi_de*depsi_dae;
    dpi_deta_dag    = -2*vepsi*vepsi*vepsi*depsi_de*depsi_dag;
    
    dpi_dgamma     = vepsi*vepsi*depsi_dg;
    dpi_dgamma2    = -2*vepsi*vepsi*vepsi*depsi_dg*depsi_dg+vepsi*vepsi*depsi_dg;
    dpi_dgamma_dae = -2*vepsi*vepsi*vepsi*depsi_dg*depsi_dae;
    dpi_dgamma_dag = -2*vepsi*vepsi*vepsi*depsi_dg*depsi_dag;
    
    dpi_dae     = vepsi*vepsi*depsi_dae;
    dpi_dae2    = -2*vepsi*vepsi*vepsi*depsi_dae*depsi_dae;
    dpi_dae_dag = -2*vepsi*vepsi*vepsi*depsi_dae*depsi_dag;
    
    dpi_dag  = vepsi*vepsi*depsi_dag;
    dpi_dag2 = -2*vepsi*vepsi*vepsi*depsi_dag*depsi_dag;
    
    // Finding Expectations
    
    //ASE_ExpFunc(ctvec,lgct,rct_vec,lgrct,pvec,0.5,nTi,vpsi,outvec,tmpctvec,tmprctvec,maxAS);
    outvec[0] = R::digamma(vpsi*pi+nBi);
    outvec[1] = R::digamma(vpsi*(1.0-pi)+nTi-nBi);
    outvec[2] = R::trigamma(vpsi*pi+nBi);
    outvec[3] = R::trigamma(vpsi*(1.0-pi)+nTi-nBi);
    
    //Plugging In:
    B = outvec[0]-outvec[1]-dgpv+dgap;
    dB_dp = pi*outvec[2]-ai_pi*outvec[3]-trpv*pi+ai_pi*trap;
    dB_dpi = outvec[2]+outvec[3]-trpv-trap;
    
    Iep(0,0) -= -vpsi*(B+vpsi*dB_dp)*dpi_dkappa;                      // kappa
    Iep(1,0) -= -vpsi*(B+vpsi*dB_dp)*dpi_deta;	                       // eta
    Iep(2,0) -= -vpsi*(B+vpsi*dB_dp)*dpi_dgamma;                      // gamma
    Ipa(0,0) -= -vpsi*(B+vpsi*dB_dp)*dpi_dae;                         // alpha_eta
    Ipa(0,1) -= -vpsi*(B+vpsi*dB_dp)*dpi_dag;                         // alpha_gamma
    Iea(0,0) -= vpsi*B*dpi_dkappa_dae+vpsi*vpsi*dB_dpi*dpi_dkappa*dpi_dae; // kappa alpha_eta
    Iea(0,1) -= vpsi*B*dpi_dkappa_dag+vpsi*vpsi*dB_dpi*dpi_dkappa*dpi_dag; // kappa alpha_gamma
    Iea(1,0) -= vpsi*B*dpi_deta_dae+vpsi*vpsi*dB_dpi*dpi_deta*dpi_dae;     // eta alpha_eta
    Iea(1,1) -= vpsi*B*dpi_deta_dag+vpsi*vpsi*dB_dpi*dpi_deta*dpi_dag;     // eta alpha_gamma
    Iea(2,0) -= vpsi*B*dpi_dgamma_dae+vpsi*vpsi*dB_dpi*dpi_dgamma*dpi_dae; // gamma alpha_eta
    Iea(2,1) -= vpsi*B*dpi_dgamma_dag+vpsi*vpsi*dB_dpi*dpi_dgamma*dpi_dag; // gamma alpha_gamma
    
    Iee(0,0) -= vpsi*B*dpi_dkappa2 + vpsi*vpsi*dB_dpi*dpi_dkappa*dpi_dkappa;       // kappa kappa
    Iee(0,1) -= vpsi*B*dpi_dkappa_deta + vpsi*vpsi*dB_dpi*dpi_dkappa*dpi_deta;     // kappa eta
    Iee(0,2) -= vpsi*B*dpi_dkappa_dgamma + vpsi*vpsi*dB_dpi*dpi_dkappa*dpi_dgamma; // kappa gamma
    Iee(1,1) -= vpsi*B*dpi_deta2 + vpsi*vpsi*dB_dpi*dpi_deta*dpi_deta;             // eta eta
    Iee(1,2) -= vpsi*B*dpi_deta_dgamma + vpsi*vpsi*dB_dpi*dpi_deta*dpi_dgamma;     // eta gamma
    Iee(2,2) -= vpsi*B*dpi_dgamma2 + vpsi*vpsi*dB_dpi*dpi_dgamma*dpi_dgamma;       // gamma gamma		
    
    Iaa(0,0) -= vpsi*B*dpi_dae2 + vpsi*vpsi*dB_dpi*dpi_dae*dpi_dae;    // alpha_eta alpha_eta
    Iaa(0,1) -= vpsi*B*dpi_dae_dag + vpsi*vpsi*dB_dpi*dpi_dae*dpi_dag; // alpha_eta alpha_gamma
    Iaa(1,1) -= vpsi*B*dpi_dag2 + vpsi*vpsi*dB_dpi*dpi_dag*dpi_dag;    // alpha_gamma alpha_gamma
    
    Ipp(0) -= vpsi*(R::digamma(vpsi)+pi*outvec[0]+ai_pi*outvec[1]-
      pi*R::digamma(piv)-ai_pi*R::digamma(ai_piv)-R::digamma(vpsi+nTi));
    Ipp(0) -= vpsi*vpsi*(R::trigamma(vpsi)+pi*pi*outvec[2]+ai_pi*ai_pi*outvec[3]-
      pi*pi*R::trigamma(piv)-ai_pi*ai_pi*R::trigamma(ai_piv)-
      R::trigamma(vpsi+nTi));
    
    // Score Function Contribution:
    ScoreVec(0,0) += vpsi*B*dpi_dae;
    ScoreVec(1,0) += vpsi*B*dpi_dag;
    
    // Checking where the NaN is introducted ( only want to call attention to the subject and their values):
    /*std::cout << "dpi_dpar: Kappa (kappa eta gamma a_eta a_gamma)\n" << dpi_dkappa2  << " // " << dpi_dkappa_deta << " // " << dpi_dkappa_dgamma << " // " << dpi_dkappa_dae << " // " << dpi_dkappa_dag << std::endl;
    std::cout << "dpi_dpar: Eta (kappa eta gamma a_eta a_gamma)\n" << dpi_dkappa_deta << " // " << dpi_deta2 << " // " << dpi_deta_dgamma << " // " << dpi_deta_dae << " // " << dpi_deta_dag << std::endl;
    std::cout << "dpi_dpar: Gamma (kappa eta gamma a_eta a_gamma)\n" << dpi_dkappa_dgamma << " // " << dpi_deta_dgamma << " // " << dpi_dgamma2 << " // " << dpi_dgamma_dae <<" // " << dpi_dgamma_dag << std::endl;
    std::cout << "dpi_dpar: A Eta (kappa eta gamma a_eta a_gamma)\n" << dpi_dkappa_dae << " // " << dpi_deta_dae << " // " << dpi_dgamma_dae << " // " << dpi_dae2 << " // " << dpi_dae_dag << std::endl;
    std::cout << "dpi_dpar: A Gamma (kappa eta gamma a_eta a_gamma)\n" << dpi_dkappa_dag << " // " << dpi_deta_dag << " // " << dpi_dgamma_dag << " // " << dpi_dae_dag << " // " << dpi_dag2 << std::endl;*/
    //std::cout << "Subject " << i << "\n Pertinent Values: nTi = " << nTi << "// Rho = " << rhosAS[i] << std::endl;
    //std::cout << "Kappa = " << kappa << " // Eta = " << eta << " // gamma = " << gamma << "\n B = " << B << " // dpi_dae = " << dpi_dae << std::endl;
    //std::cout << "OutVec[0//1//2//3] = " << outvec[0] << "//" << outvec[1] << "//" << outvec[2] << "//" << outvec[3] << std::endl;
  }
  
  for(i=(nAS_0+nAS_1);i<nAS;i++){
    // Finding Expectations
    nTi = ceil(ex2[4+i]-0.5);
    nBi = ceil(ex2[nASp4+i]-0.5);
    
    //ASE_ExpFunc(ctvec,lgct,rct_vec,lgrct,pvec,0.5,nTi,vpsi,outvec,tmpctvec,tmprctvec,maxAS);
    outvec[0] = R::digamma(vpsi*0.5+nBi);
    outvec[1] = R::digamma(vpsi*0.5+nTi-nBi);
    outvec[2] = R::trigamma(vpsi*0.5+nBi);
    outvec[3] = R::trigamma(vpsi*0.5+nTi-nBi);
    
    //Plugging In:
    Ipp(0) -= vpsi*(R::digamma(vpsi)-R::digamma((0.5*vpsi))-
      R::digamma(vpsi+nTi)+0.5*(outvec[0]+outvec[1]));
    Ipp(0) -= vpsi*vpsi*(R::trigamma(vpsi)-0.5*R::trigamma(0.5*vpsi)-
      R::trigamma(vpsi+nTi)+0.25*(outvec[2]+outvec[3]));
  }
  
  Iaa(1,0) = Iaa(0,1);
  Iee(1,0) = Iee(0,1);
  Iee(2,0) = Iee(0,2);
  Iee(2,1) = Iee(1,2);
  
  // Reporting the Values of matrices so that I can utilize them for checking
  /*std::cout << "Iee matrix: \n" << Iee << std::endl;
  std::cout << "Iep matrix: \n" << Iep << std::endl;
  std::cout << "Iea matrix: \n" << Iea << std::endl;
  std::cout << "Ipp matrix: \n" << Ipp << std::endl;
  std::cout << "Ipa matrix: \n" << Ipa << std::endl;
  std::cout << "Iaa matrix: \n" << Iaa << std::endl;
  std::cout << "Score Vect: \n" << ScoreVec << std::endl;*/
  
  // Final Information Matrix Components//
  M1.block(0,0,M,M)      = Ibb;
  M1.block(0,M,M,3)      = Ibe;
  M1.block(M,0,3,M)      = Ibe.transpose();
  M1.block(0,(M+3),M,1)  = Ibt; 
  M1.block((M+3),0,1,M)  = Ibt.transpose();
  M1.block<3,3>(M,M)     = Iee;
  M1.block<3,1>(M,(M+3)) = Iet;
  M1.block<1,3>((M+3),M) = Iet.transpose();
  M1.block<3,1>(M,(M+4)) = Iep;
  M1.block<1,3>((M+4),M) = Iep.transpose();
  M1((M+3),(M+3))        = Itt(0);
  M1((M+4),(M+4))        = Ipp(0);
  
  M2.block<3,2>(M,0)     = Iea;
  M2.block<1,2>((M+4),0) = Ipa;
  
  OIMat.block(0,0,(M+5),(M+5)) = M1;
  OIMat.block(0,(M+5),(M+5),2) = M2;
  OIMat.block((M+5),0,2,(M+5)) = M2.transpose();
  OIMat.block((M+5),(M+5),2,2) = Iaa;
  
  // Checking the formation of M1 and M2:
  //std::cout << "M1 Matrix: \n" << M1 << std::endl;
  //std::cout << "M2 Matrix: \n" << M2 << std::endl;
  
  // Checking Invertibility
  Eigen::JacobiSVD<Eigen::MatrixXd> svd1(M1);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd2(Iaa-M2.transpose()*M1.inverse()*M2);
  double cond1 = svd1.singularValues()(0) / svd1.singularValues()(svd1.singularValues().size()-1);
  double cond2 = svd2.singularValues()(0) / svd2.singularValues()(svd2.singularValues().size()-1);
  
  // Checking Positive Definiteness:
  Eigen::LLT<Eigen::MatrixXd> lltOfA(OIMat);
  
  if(lltOfA.info()==Eigen::NumericalIssue){
    *score_val = -6.0;
    *pval = 1.0;
  } else if((cond1>(1/DBL_EPSILON))||(cond2>(1/DBL_EPSILON))){
    *score_val = -5.0;
    *pval = 1.0;
  } else {
    *score_val = ScoreVec.transpose()*((Iaa-M2.transpose()*M1.inverse()*M2).inverse())*ScoreVec;
    *pval = R::pchisq(*score_val,2,0,0);
  }
  
  if(std::isnan(*score_val)){
    *score_val = -6.0;
    *pval = 1.0;
  }
  // This way if the score is -5.0, we know that we have invertibility issues.
}
