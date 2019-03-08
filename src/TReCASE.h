#ifdef __cplusplus
extern "C"{
#endif

	double neg_logLik_TS(int n, double* para, void* ex, SEXP x1);
	void neg_gradLik_TSnh(int n, double* para, double* gr, void* ex, SEXP x1);
	void KEG_optroutine(int npara,int lmm,double* initPara,double* lower,
				double* upper, int* nbd, double* Fmin, int* fail, double* ex,
				double factr, double pgtol, int* fncount, int* grcount, int maxit1,
				char* msg, int nREPORT, double* wa, int* iwa, double* g1, SEXP x1);
	void KEG_boptroutine(int npara,double* initPara,double* Fmin, int maxit1,
                      int* mask, double abstol, double reltol, int nReport,
                      double* ex, int* fncount, int*grcount, int* fail,
                      double* gv, double* tv, double* Xv, double* cv, double** Bv,
                      SEXP x1);
	void KEG_optroutine2(int npara,int lmm,double* initPara,double* lower,
                     double* upper, int* nbd, double* Fmin, int* fail, double* ex,
                     double factr, double pgtol, int* fncount, int* grcount, int maxit1,
                     char* msg, int nREPORT, double* wa, int* iwa, double* g1, SEXP x1);
	void KEG_boptroutine2(int npara,double* initPara,double* Fmin, int maxit1,
                      int* mask, double abstol, double reltol, int nReport,
                      double* ex, int* fncount, int*grcount, int* fail,
                      double* gv, double* tv, double* Xv, double* cv, double** Bv, SEXP x1);
	void psi_optroutine(int npara,int lmm,double* initPara,double* lower,
				double* upper, int* nbd, double* Fmin, int* fail, double* ex,
				double factr, double pgtol, int* fncount, int* grcount, int maxit1,
				char* msg, int nREPORT, double* wa, int* iwa, double* g1, SEXP x1);
	void psi_boptroutine(int npara,double* initPara,double* Fmin, int maxit1,
                      int* mask, double abstol, double reltol, int nReport,
                      double* ex, int* fncount, int*grcount, int* fail,
                      double* gv, double* tv, double* Xv, double* cv, double** Bv, SEXP x1);
	
#ifdef __cplusplus
}
#endif
