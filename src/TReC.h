#ifdef __cplusplus
extern "C" {
#endif

double loglikNB_pmy (double phi, int n, double* mu, int* y);
double neg_logLTReC_TS(int n, double* para, void* ex, SEXP x1);
double neg_gradTReC_TS(int n, double* para, double* gr, double* Hessian, void* ex, SEXP x1);

#ifdef __cplusplus
}
#endif
