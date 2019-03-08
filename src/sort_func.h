void DJ_Rcpp_sort(double* Xs, double* rhos, int* y, int* y1, int* y2, int* z,
				  double* Xt, double* rhot, int* yt, int* y1t, int* y2t,
				  int* n0, int* n1, int* n3, int* n4, int* nAS_all, int* nAS_h,
				  const int n_subj, const int M, const int min_ASE_Total,
				  int* idx0, int* idx1, int* idx3, int* idx4);
void set_extras(int* y, int* y1, int* y2, double* rhos, int n0, int n1, 
				int n_subj, int nAS_all, int nAS_h, 
				int min_ASE_Total, int min_nASE,
				double* ex1, double* ex2, double* rhos_ASall);
