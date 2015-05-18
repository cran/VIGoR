\name{hyperpara}
\alias{hyperpara}
\title{
Calculation of hyperparameter values
}
\description{
This function determines the hyperparameter values of regression methods, based on several assumptions. 
}
\usage{
hyperpara(Geno, Mvar, Method = c("BL", "EBL", "wBSR", "BayesB", "BayesC", "SSVS", "MIX"),
          Kappa, A = 0.9, Xtype="Geno", f = 0, BL.Phi = 1, EBL.Phi = 0.1, 
          EBL.Omega = 0.1, Psi = 1, Nu = 5, Printinfo = FALSE)
}
\arguments{
  \item{Geno}{
  An (N x P) matrix, where N and P denote the number of individuals and markers, respectively. Marker genotypes should be coded as 0 (AA), 1 (AB), and 2 (BB). Because missing values are not allowed, missing genotypes should be imputed before analysis. Doubles between 0 and 2 are allowed. Marker genotypes are used to calculate sum (2*Q*(1-Q)*(1+f)) where Q is a vector of allele frequencies and f is the inbreeding coefficient.
}
  \item{Mvar}{
  A scalar denoting the assumed proportion of phenotypic variance that can be explained by the markers. Mvar is < 1.0 in BL and EBL, and <= 1.0 in the other methods.
}
  \item{Method}{
  One of the seven regression methods
}
  \item{Kappa}{
  A scalar or vector containing the assumed proportion of markers with non-zero effects. In MIX and SSVS, Kappa is < 1; in the other methods, 0 < Kappa <= 1.
}
  \item{A}{
  In the SSVS and MIX models. A is a scalar or vector denoting the assumed proportion of Mvar explained by the markers assigned to the normal prior distribution with the larger variance. For example, given Mvar = 0.5, A = 0.9, and Kappa = 0.01, then 0.45 (0.5 x 0.9) of the phenotypic variance is assumed to be explained by (P x 0.01) markers, and 0.05 (0.5 x (1-0.9)) of variance is assumed to be explained by (P x (1 - 0.01)) markers.
}
  \item{Xtype}{
  Allowed Xtypes are "Geno" and "Var". Enter "Geno" when Geno contains the marker genotypes and "Var" when Geno contains variables other than marker genotypes.
}
  \item{f}{
  A scalar representing the inbreeding coefficient. Enter 1 for inbred species. 
}
  \item{BL.Phi}{
  A scalar or vector containing Phi values of BL.
}
  \item{EBL.Phi}{
  A scalar or vector containing Phi values of EBL.
}
  \item{EBL.Omega}{
  A scalar or vector containing Omega values of EBL.
}
  \item{Psi}{
  A scalar or vector containing Psi values of EBL.
}
  \item{Nu}{
  A scalar or vector containing Nu values of wBSR, BayesB, BayesC, SSVS, and MIX.
}
  \item{Printinfo}{
  Specifies whether to print information (TRUE) or not (FALSE).
}
}
\details{
  To run vigor, users must specify the following hyperparameter values.\cr
  \itemize{
    \item BL: Phi, Omega
    \item EBL: Phi, Omega, Psi, Theta
    \item wBSR: Nu, S2, Kappa
    \item BayesB: Nu, S2, Kappa      
    \item BayesC: Nu, S2, Kappa  
    \item SSVS: c, Nu, S2, Kappa 
    \item MIX: c, Nu, S2, Kappa
    }
  This function calculates the Omega of BL; Theta of EBL; S2 of wBSR, BayesB, and BayesC; and c and S2 of SSVS and MIX. Mvar, Kappa, and the other hyperparameters required by each method are specified by the user. The SSVS and MIX models also require A. The definitions of Mvar, Kappa, and A are intuitively understandable and relatively easy to specify (see Arguments). We recommend the default values for the other hyperparameters. When the arguments of hyperpara are vectors, all value-combinations are returned as a matrix. The hyperparameters are explained in the details of vigor and in the pdf manual of VIGoR (Onogi 2015).
}
\value{
This function returns a vector when yielding a single hyperparameter set, and a matrix when yielding multiple hyperparameter sets. The rows and columns of the matrix correspond to the sets (value combinations) and the hyperparameters, respectively. See examples below.\cr
}
\seealso{
  vigor
}
\references{
  Onogi & Iwata, VIGoR: variational Bayesian inference for genome-wide regression, in prep.\cr
  Onogi A, 2015, Documents for VIGoR (May 2015).https://github.com/Onogi/VIGoR
}
\examples{
#data
data(sampledata)
dim(Geno)#100 samples and 1000 markers
unique(Geno[1:(100*1000)])#coded as 0, 1, 2

#A single Kappa value is assumed for BL. A vector is returned.
hyperpara(Geno,0.5,"BL",0.01,Printinfo=TRUE)

#Phi is set to 1 as default. To change Phi, use BL.Phi.
hyperpara(Geno,0.5,"BL",0.01,BL.Phi=5)

#Calculate multiple hyperparameter value sets of BayesC assuming that Kappa is 0.01, 0.1, and 1.
#A matrix is returned.
hyperpara(Geno,0.5,"BayesC",c(0.01,0.1,1))

#The output vector can be used as the argument of vigor
Result<-vigor(Pheno$Height,Geno,"wBSR",hyperpara(Geno,0.5,"wBSR",0.01))

#Calculate multiple hyperparameter sets of SSVS 
#assuming that Mvar is 0.5, Kappa is 0.01 and 0.1, and A is 0.9 and 0.99.
hyperpara(Geno,0.5,"SSVS",c(0.01,0.1),c(0.9,0.99))
#2 x 2 sets are created.

#Calculate hyperparameter values of BayesB 
#assuming that Mvar is 0.5, and Kappa is 0.01. Inbred lines are analyzed.
hyperpara(Geno,0.5,"BayesB",0.01,f=1)

#Calculate hyperparameter values of EBL
#assuming that Mvar is 0.5, and Kappa is 0.01. Consider marker genotypes as general variables.
hyperpara(Geno,0.5,"EBL",0.01,Xtype="Var")
}


