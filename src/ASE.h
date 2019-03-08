#ifdef __cplusplus
extern "C"{
#endif

double loglikBB(int n, int* nB, int* nTotal, double* pis, double psi);
double neg_logLASE_TS(int n, double* para, void* ex, SEXP x1);
void neg_gradASE_TS(int n, double* para, double* gr, double* Hessian,void* ex, SEXP x1);
void compute_pis(double kappa, double eta, double gamma, double* rhos, double* ex);
double logLASE_psi(int n, double* para, void* ex, SEXP x1);
double compute_BBHz(int nAS_Hz, double* ex, double psi_curr);
void gradASE_psi(int n, double* para, double* gr, void* ex, SEXP x1);

#ifdef __cplusplus
}
#endif
