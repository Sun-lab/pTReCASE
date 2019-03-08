void compute_nu0(double* X, int n, int M, double* Beta, double* ex);
void comp_offset(double kappa, double eta, double gamma, double* ex1,
				 double* offset);
double phi_update(double phi_init,int n_subj, double* y, double* mu, int limit, double eps);
double varfunc(double mu,double phi);
double dfunc(double mu);
void phi_initialize(double* y, double* X, int n, int M, double* offset,
					  double* Beta, double* phi_curr);
int IRLS_NB_fit(double* y, int* yint, double* X, int n, int M, double* offset,
				 int iter, double eps, double convg, int init_beta_,
				 double* phi_out, double* logLNB, double* Beta_out);
int IRLS_Pois_fit(double* y, double* X, int n, int M, double* offset,
				 int iter, double eps, double convg, double* Beta_out);
