#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <RcppEigen.h>

/******************************************************************
 *
 ******************************************************************/
 
void DJ_Rcpp_sort(double* Xs, double* rhos, int* y, int* y1, int* y2, int* z,
				  double* Xt, double* rhot, int* yt, int* y1t, int* y2t,
				  int* n0, int* n1, int* n3, int* n4, int* nAS_all, int* nAS_h,
				  const int n_subj, const int M, const int min_ASE_Total,
				  int* idx0, int* idx1, int* idx3, int* idx4){
	// Declare RcppEigen classes:
	Eigen::Map<const Eigen::MatrixXd> Ximat(Xs,n_subj,M);
	Eigen::Map<Eigen::MatrixXd> Xomat(Xt,n_subj,M);
	
	// Declare C++ objects:
	int i=0, i_prime=0, AS_ind = 0, nTot_i = 0;
	*n0 = 0;
	*n1 = 0;
	*n3 = 0;
	*n4 = 0;
	
	*nAS_all = 0;
	*nAS_h   = 0;
	
	// Index Finding:
	for(i=0;i<n_subj;i++){
		AS_ind = 0;
		nTot_i = y1[i]+y2[i];
		if(nTot_i>=min_ASE_Total){
			*nAS_all += 1;
			AS_ind = 1;
		}
		
		if(z[i]==0){
			idx0[*n0] = i;
			*n0 += 1;
		} else if(z[i]==1){
			idx1[*n1] = i;
			*n1 += 1;
			*nAS_h += AS_ind;
		} else if(z[i]==3){
			idx3[*n3] = i;
			*n3 += 1;
			*nAS_h += AS_ind;
		} else if(z[i]==4){
			idx4[*n4] = i;
			*n4 += 1;
		}
	}
	
	// Placing Output:
	for(i=0;i<*n0;i++){
		Xomat.row(i) = Ximat.row(idx0[i]);
		rhot[i] = rhos[idx0[i]];
		
		yt[i]  = y[idx0[i]];
		y1t[i] = y1[idx0[i]];
		y2t[i] = y2[idx0[i]];	
		
		idx0[i] = -1;
	}
	
	i_prime = 0;
	for(i=*n0;i<(*n0+*n1);i++){
		Xomat.row(i) = Ximat.row(idx1[i_prime]);
		rhot[i] = rhos[idx1[i_prime]];
		
		yt[i]  = y[idx1[i_prime]];
		y1t[i] = y1[idx1[i_prime]];
		y2t[i] = y2[idx1[i_prime]];
		
		idx1[i_prime] = -1;
		i_prime += 1;
	}
	
	i_prime = 0;
	for(i=(*n0+*n1);i<(*n0+*n1+*n3);i++){
		Xomat.row(i) = Ximat.row(idx3[i_prime]);
		rhot[i] = rhos[idx3[i_prime]];
		
		yt[i]  = y[idx3[i_prime]];
		y1t[i] = y2[idx3[i_prime]];
		y2t[i] = y1[idx3[i_prime]];
		
		idx3[i_prime] = -1;
		i_prime += 1;
	}
	
	i_prime = 0;
	for(i=(*n0+*n1+*n3);i<n_subj;i++){
		Xomat.row(i) = Ximat.row(idx4[i_prime]);
		rhot[i] = rhos[idx4[i_prime]];
		
		yt[i]  = y[idx4[i_prime]];
		y1t[i] = y1[idx4[i_prime]];
		y2t[i] = y2[idx4[i_prime]];
		
		idx4[i_prime] = -1;
		i_prime += 1;
	}
}

/******************************************************************
 * Setting up the extra Parameters
 ******************************************************************/
void set_extras(int* y, int* y1, int* y2, double* rhos, int n0, int n1, 
				int n_subj, int nAS_all, int nAS_h, 
				int min_ASE_Total, int min_nASE,
				double* ex1, double* ex2, double* rhos_ASall){
	// Declaring values:
	int i = 0, idxy = 0, idxr1 = 0, idxn = 0, idxnt = 0,
		idxnb = 0, idxr2 = 0;
	int nTot_i = 0, i_prime1 = 0, i_prime2 = 0, trpn_p9 = 3*n_subj+9;
	
	// ex1 Preliminaries:
	ex1[0] = n_subj;       // n_subj
	ex1[2] = n0;           // n0
	ex1[3] = n1;           // n1
	ex1[4] = n_subj-n0-n1; // n2
	ex1[5] = 0;            // H0
	ex1[6] = 0;			   // Hess
	ex1[7] = 0.5;          // phi
	ex1[8] = 0.5;          // psi
	
	// Filling in the extras:
	for(i=0;i<n0;i++){
		// Fill in ex1:
		idxy  = 9+i;
		idxr1 = idxy+n_subj;
		idxn  = idxr1+n_subj;
		
		ex1[idxy]  = y[i];
		ex1[idxr1] = rhos[i];
		ex1[idxn]  = 0;
		
		// Fill in ex2:
		nTot_i = y1[i]+y2[i];
		if(nTot_i>=min_ASE_Total){
			idxnt = 4+i_prime2;
			idxnb = idxnt+nAS_all;
			idxr2 = idxnb+nAS_all;
			
			ex2[idxnt] = nTot_i;
			ex2[idxnb] = y2[i];
			ex2[idxr2] = rhos[i];
			rhos_ASall[i_prime2] = rhos[i];
			
			i_prime2 += 1;
		}
	}
	
	ex2[1] = i_prime2;
	
	for(i=n0;i<(n0+n1);i++){
		idxy  = 9+i;
		idxr1 = idxy+n_subj;
		idxn  = idxr1+n_subj;
		
		ex1[idxy]  = y[i];
		ex1[idxr1] = rhos[i];
		ex1[idxn]  = 0;
		
		// Fill in ex2 and ex1 AS:
		nTot_i = y1[i]+y2[i];
		if(nTot_i>=min_ASE_Total){
			idxnt = 4+i_prime2;
			idxnb = idxnt+nAS_all;
			idxr2 = idxnb+nAS_all;
			
			ex2[idxnt] = nTot_i;
			ex2[idxnb] = y2[i];
			ex2[idxr2] = rhos[i];
			rhos_ASall[i_prime2] = rhos[i];
			
			idxnt = trpn_p9 + i_prime1;
			idxnb = idxnt+nAS_h;
			idxr2 = idxnb+nAS_h;
			
			ex1[idxnt] = nTot_i;
			ex1[idxnb] = y2[i];
			ex1[idxr2] = rhos[i];
			
			i_prime2 += 1;
			i_prime1 += 1;
		}
	} 
	
	ex2[2] = i_prime2 - ex2[1];
	
	for(i=(n0+n1);i<n_subj;i++){
		idxy  = 9+i;
		idxr1 = idxy+n_subj;
		idxn  = idxr1+n_subj;
		
		ex1[idxy]  = y[i];
		ex1[idxr1] = rhos[i];
		ex1[idxn]  = 0;
		
		// Fill in ex2:
		nTot_i = y1[i]+y2[i];
		if(nTot_i>=min_ASE_Total){
			idxnt = 4+i_prime2;
			idxnb = idxnt+nAS_all;
			idxr2 = idxnb+nAS_all;
			
			ex2[idxnt] = nTot_i;
			ex2[idxnb] = y2[i];
			ex2[idxr2] = rhos[i];
			rhos_ASall[i_prime2] = rhos[i];
			
			i_prime2 += 1;
		}
	}
	
	ex1[1] = i_prime1;
	ex2[3] = i_prime2 - ex2[2] - ex2[1];
	
	// ex2 beginning values:
	ex2[0] = i_prime2;
}

/**** Parameter Descriptions ****/
 /*
 // denotes a line break
 ex1 : n,nAS,n0,n1,n2,H0,Hess, phi, psi// [9 elts]
       Y   // [n elts]
       rho // [n elts]
       nu  // [n elts]
 
 ex2 : nAS,nAS_0,nAS_1,nAS_2// [4 elts]
       nTotal // [nAS elts]
       nB     // [nAS elts]
       pi     // [nAS elts]
*/
