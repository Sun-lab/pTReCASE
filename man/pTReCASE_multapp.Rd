\name{pTReCASE_multapp}
\alias{pTReCASE_multapp}
\title{
pTReCASE Multiple Application
}
\description{
This function takes in total read count (TReC) and allele specific expression (ASE) data for several genes, genotype data for several potential SNPs, and tumor purity estimates. It returns pTReCASE model results for each local gene-SNP pair. Function is designed to be applied to only a single chromosome at a time.
}
\usage{
pTReCASE_multapp(Y, Y1, Y2, Z, X, rhos, F_out, 
                geno_pos, gene_start, gene_end, 
                min_ASE_Total = 8, min_nASE = 10, 
                eps = 1e-05, convg = 1e-05, 
                maxIter = 200, maxIter2 = 200, 
                Perform_CT_co = 0.001, CT_pco = 0.001, 
                Obs_Score = 1, useASE = 1, 
                cis_window = 1e+05, gene_add = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  A matrix of size N-by-G of TReC values for G genes and N subjects. Data is stored as one column per gene and one row per subject.
}
  \item{Y1}{
  A matrix of size N-by-G of read counts mapping to arbitrarily defined allele 1. Data is stored as one column per gene and one row per subject.
}
  \item{Y2}{
  A matrix of size N-by-G of read counts mapping to arbitrarily defined allele 1. Data is stored as one column per gene and one row per subject.
}
  \item{Z}{
  A matrix of size N-by-S of genotype data for N subjects at S potential SNPs. Data is stored as one row per subject, one column per SNP. A single item of Z takes values 1,2,3, or 4 corresponding to genotypes AA, AB, BA, and BB.
}
  \item{X}{
  A matrix of size N-by-M of covariates for N subjects and M measured covariates. Data is stored as one column per covariate, one row per subject.
}
  \item{rhos}{
  A vector of length N containing the estimated tumor purities for each subject. The order of data in \code{rhos} must match the order of the rows in Y, Y1, Y2, and Z.
}
  \item{F_out}{
  Vector of length 2. First element contains a character string that records the file name at which the output is stored. The other element is a place holder. 
}
  \item{geno_pos}{
  A vector of length S containing the location of the each SNP in the matrix Z. Order in geno_pos must match the order of the columns in Z.
}
  \item{gene_start}{
  A vector of length G containing the starting location for each gene in matrices Y, Y1, Y2. Order must match the order of columns in Y, Y1, Y2.
}
  \item{gene_end}{
  A vector of length G containing the ending location for each gene in matrices Y, Y1, Y2. Order must match the order of columns in Y, Y1, Y2.
}
  \item{min_ASE_Total}{
  The minimum number of Allele Specific reads a subject must have to be included in the ASE analysis. Recommended value of 8.
}
  \item{min_nASE}{
  The minimum number of subjects that must exceed the min_ASE_Total threshold in order for ASE data to be included.
}
  \item{eps}{
  A goodness of fit value for NB regression. If computed model weights are less than eps, then there are fit issues. Default of 1e-5.
}
  \item{convg}{
  Specifies whether pTReCASE model has converged. Maximum change in any of kappa, eta, gamma, phi or psi. Default value of 1e-5.
}
  \item{maxIter}{
  Maximum number of iterations for pTReC and pTReCASE model fits. Default value of 200.
}
  \item{maxIter2}{
  Maximum number of iterations for Negative Binomial regression model fits. Default value of 200.
}
  \item{Perform_CT_co}{
  The p-value for a tumor-specific eQTL below which triggers the computation of the Cis-Trans score test for model identification. Default value is 0.001.
}
  \item{CT_pco}{
  The p-value cutoff for the Cis-Trans Score test which determines whether to use the pTReC or pTReCASE model. If p<CT_pco, use the pTReC model. If p>CT_pco, then use the pTReCASE model. Results for both models are output to the final dataset.  
}
  \item{Obs_Score}{
  Indicator of whether or not the observed information matrix should be used to compute the Cis-Trans Score test (1=Observed Information, 0=Approximate Fisher's Information). The default is 1.
}
  \item{useASE}{
  Indicator of whether or not to use ASE data for any gene (1=Use ASE,0=TReC only). Default is set to 1.
}
  \item{cis_window}{
  Value determining local eQTL. A SNP is considered local to a gene if and only if it is located within the gene body or within cis_window bp from either the start or end of the gene.
}
  \item{gene_add}{
  Value to assist in the storage of results when tests are conducted in parallel across different submitted jobs. Say we split Y,Y1, and Y2 into 3 sub-matrices each consisting of columns [1,..,10], [11,...,20], and [21,30].  
}
}
\value{
All output is stored in the file specified in F_out[1] in tab-delimited format. For each Gene-SNP pair there are at most two records: 

(1) Results from the pTReCASE model if sufficient sample size for ASE analysis.

(2) Results from pTReC model if pTReCASE fails the Cis-Trans score test.

$SNP: Potential SNP modeled and stored using 0-index. Thus, if record states SNP = 4, this corresponds to SNP in the (4+1 = 5)th column of Z.

$Gene: The modeled gene, stored using 0-index. Thus, if the record states Gene=4, this corresponds to the Gene in the (4+1 = 5)th column of Y, Y1, and Y2.

$Tot_Dur: The total duration, in seconds, of the model fit.

$KEG_Dur: Total time (in seconds) spent updating the Kappa, Eta, and Gamma parameters.

$NB_Dur: Total time (in seconds) spent updating the covariate parameters Beta and Phi

$Psi_Dur: Total time (in seconds) spent updating the Psi parameter which models overdispersion in the ASE data.

$LRT_Eta: Value of the test statistic for the likelihood ratio test of Eta=1 (no normal eQTL).

$P_Eta: P-value from the likelihood ratio test of Eta=1 (no normal eQTL).

$LRT_Gamma: Value of the test statistic for the likelihood ratio test of Gamma=1 (no tumor eQTL).

$P_Gamma: P-value from the likelihood ratio test of Gamma=1 (no tumor eQTL).

$Lik_Full: Likelihood of the full model fit

$Lik_Eta: Likelihood of the Eta=1 model fit

$Lik_Gamma: Likelihood of the Gamma=1 model fit.

$Kappa: Estimate of the Kappa parameter in the full model fit. (Tumor Overexpression effect)

$Eta: Estimate of the Eta parameter in the full model fit. (Normal eQTL effect)

$Gamma: Estimate of the Gamma parameter in the full model fit. (Tumor eQTL effect)

$Phi: Estimate of overdispersion in the Negative Binomial model for TReC.

$Psi: Estimate of overdispersion in Beta-Binomial model for ASE in the full model. Is set equal to -1 if the model fit is pTReC. >0 if pTReCASE

$Fail_Full: Failure Indicator for Full Model fit.

$Fail_Full_Eta: Failure Indicator for full model fit initialized at the Eta=1 model parameters.

$Fail_Full_Gamma: Failure indicator for the full model fit initialized at the Gamma=1 model parameters.

$Fail_Eta: Failure indicator for the model fit under hypothesis Eta=1

$Fail_Gamma: Failure indicator for the model fit under hypothesis Gamma=1

[Model Fit Failure Codes]:
 - 0: Converged
 - 1: Optimization Routine Failure (likelihood decrease)
 - 2: Iteration Limit Reached
 - 3: KEG optimization failure causing a cessation in parameter updates
 - 4: Results from a call of 100 in the Full, Eta=1, or Gamma=1 model Fits.
 - 100: Optimization routine in LBFGS (infinite likelihood)
 

$Coeff_J: Coefficient estimate for the J-th covariate in the negative binomial TReC model, stored in the J-th column of matrix X.

$CT_Score: The value of the test statistic for the Cis-Trans Score Test

$CT_Pval: The p-value of the Cis-Trans Score test

$CT_Fail: Cis-Trans score test failure indicator.

[Cis-Trans Score Test Failure Codes]
- 0: Used Observed Information in score test
- 1: Used Expected Fisher's Information in score test
- 2: Invertibility issues with the Hessian in both the Observed and Expected Information matrices.
- 5: CisTrans Score test not conducted.
}
\references{
Douglas Wilson, Joseph Ibrahim, Wei Sun "Mapping Tumor-Specific Expression QTLs in Impure Tumor Samples", JASA, in revision. 
}
\author{
Douglas Roy Wilson, Jr.
}
\examples{
# Change Working Directory on your personal machine
# setwd()
data("calibrationData")
currFile = "./Calibration_Out.txt"

pTReCASE_multapp(Y = as.matrix(calibrationData$y),
                 Y1 = as.matrix(calibrationData$y1),
                 Y2 = as.matrix(calibrationData$y2),
                 Z = as.matrix(calibrationData$z),
                 X = as.matrix(calibrationData$Xs),
                 rhos = calibrationData$rhos,
                 maxIter = 400, convg=1e-5,
                 geno_pos = c(1), gene_start = c(1),
                 gene_end = c(1), useASE=1, 
                 cis_window = 100000, Obs_Score = 1,
                 F_out = currFile, Perform_CT_co=0.5)

pTReCASE_cal = read.table(file = currFile, header = TRUE, 
                          sep = "\t", stringsAsFactors = FALSE, 
                          fill = TRUE)
}
