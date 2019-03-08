#' Calibration data for the pTReCASE package
#'
#' A simulated dataset containing total read counts (TReC),
#' Allele-specific read counts (AS) at a single gene across 500 subjects. 
#' This dataset is provided to demonstrate the use of pTReCASE_multapp.
#'
#' These data were simulated from a scenario where Kappa=1.2, Eta = 1.0,
#' Gamma = 1.8, psi=0.2, and phi=0.2. Running pTReCASE on this data should
#' provide: (Kappa = 1.17, Eta = 0.99, Gamma = 1.64, psi = 0.18, phi = 0.22).
#'
#' @format A data frame of 6 elements: \describe{ \item{y,y1,y2}{Vectors of the
#' simulated read counts from a single gene (y = TReC, y1 = allele-specific reads 
#' mapping to haplotype 1, y2 = allele-specific reads mapping to haplotype 2)}
#' \item{z}{Vector of simulated genotypes (1,2,3,4)}
#' \item{rhos}{Simulated tumor purities}
#' \item{Xs}{Intercept-only covariate matrix}}
"calibrationData"
