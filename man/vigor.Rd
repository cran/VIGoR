\name{vigor}
\alias{vigor}
\title{Variational Bayesian inference for genome-wide regression}
\description{
This function performs Bayesian genome-wide regression using variational Bayesian algorithms.
The available regression methods are Bayesian lasso (BL), extended Bayesian lasso (EBL), 
weighted Bayesian shrinkage regression (wBSR), BayesB, BayesC, 
stochastic search variable selection (SSVS), and Bayesian mixture model (MIX).
}
\usage{
vigor(Pheno, Geno, Method = c("BL", "EBL", "wBSR", "BayesB", "BayesC", "SSVS", "MIX"), 
      Hyperparameters, Function = "fitting", Nfold = 10, CVFoldTuning = 5, 
      Partition=NULL, Covariates = "Intercept", Threshold = 2+log10(ncol(Geno)), 
      Maxiterations=1000, RandomIni=FALSE, Printinfo=TRUE) 
}
\arguments{
  \item{Pheno}{
  An N-length Vector of response variables (e.g., phenotypic values), where N is the number of individuals. Missing data (coded as NA) are allowed.
}
  \item{Geno}{
  An N x P matrix of marker genotypes, where P is the number of markers. Both integers and doubles can be used. Because missing values are not allowed, missing genotypes should be imputed before analysis.
}
  \item{Method}{
  String representing the selected regression method. See details below.
}
  \item{Hyperparameters}{
  A vector or matrix of hyperparameter values. When multiple combinations of hyperparameter values (hyperparameter sets) are used, the columns of the matrix correspond to the hyperparameters, and the rows correspond to the combinations (sets). See details below.}
  \item{Function}{
  One of the strings "fitting", "tuning", and "cv". See details below.
  } 
  \item{Nfold}{
  An integer value. When n > 1, n-fold cross-validation (CV) is performed on randomly partitioned individuals. When the integer is -1, leave-one-out CV is conducted. When the integer is -9, the partitioning for the CV is defined by the argument Partition. Used when Function = "cv".
  } 
  \item{CVFoldTuning}{
  An integer specifying the fold number of the CV in hyperparameter tuning. Used when Function = "cv" or "tuning" and the number of hyperparameter sets (number of rows of Hyperparameters) > 1.
  } 
  \item{Partition}{
  A matrix defining the partitions of CV. See details and examples below. Used when Function = "cv" or "tuning" and Nfold = -9.
  }
  \item{Covariates}{
  If Covariates = "Intercept", the intercept is automatically added to the regression model. An N x F matrix where F is the number of covariates also can be input as the covariates. Both integers and doubles are permitted. Missing values are not allowed. See details below.
}
  \item{Threshold}{
  Specifies the convergence threshold. Calculation terminates if the convergence metric is < 1e-Threshold. See the pdf manual of VIGoR (Onogi 2015) for the metric.
}
  \item{Maxiterations}{
  Maximum number of iterations.
}
  \item{RandomIni}{
  If TRUE, the initial values of the marker effects are randomly determined. Otherwise, they are set to 0.
}
  \item{Printinfo}{
  If TRUE, print the run information to the console.
}
}
\details{
  For details of vigor, the user is referred to the pdf manual (Onogi 2015).\cr\cr
  \bold{Regression methods}\cr  
  Vigor assumes the following linear model for individual i;\cr\cr  
  Pheno[i] = sum(Covariates[i,]*Alpha) + sum(Gamma*Geno[i,]*Beta) + Ei\cr\cr  
  where Pheno[i] is the response variable, Alpha contains the covariate coefficients, and Gamma contains the variables that indicate whether the corresponding markers are included in the model (1) or not (0). Beta contains the marker effects, and Ei is the residual. Gamma is set to 1 except in wBSR, which infers Gamma. Ei is assumed to follow a Normal (0, 1/Tau02) distribution, where Tau02 is assumed to follow 1/Tau02.\cr
  The methods assume different prior distributions of the marker effect p (Beta[p]).\cr
  \itemize{
    \item BL\cr
    Beta[p] ~ normal (0, sqrt(1/Tau2[p]/Tau02))\cr
    Tau2[p] ~ inverse gamma (1, Lambda2/2)\cr
    Lambda2 ~ gamma (Phi, Omega)
 
    \item EBL\cr
    Beta[p] ~ normal (0, sqrt(1/Tau2[p]/Tau02))\cr   
    Tau2[p] ~ inverse gamma (1, Delta2*Eta2[p]/2)\cr    
    Delta2 ~ gamma (Phi, Omega)\cr    
    Eta2[p] ~ gamma (Psi, Theta)

    \item wBSR\cr
    Beta[p] ~ normal (0, sqrt(Sigma2[p]))\cr    
    Sigma2[p] ~ scaled-inverse-chi square (Nu, S2)\cr    
    Gamma[p] ~ Bernoulli (Kappa)

    \item BayesB\cr
    Beta[p] ~ normal (0, sqrt(Sigma2[p])) if Rho[p]=1, and 0 if Rho[p]=0\cr    
    Sigma2[p] ~ scaled-inverse-chi square (Nu, S2)\cr    
    Rho[p] ~ Bernoulli (Kappa)
    
    \item BayesC\cr
    Beta[p] ~ normal (0, sqrt(Sigma2)) if Rho[p]=1, and 0 if Rho[p]=0\cr    
    Sigma2 ~ scaled-inverse-chi square (Nu, S2)\cr    
    Rho[p] ~ Bernoulli (Kappa)
   
    \item SSVS\cr 
    Beta[p] ~ normal (0, sqrt(Sigma2)) if Rho[p]=1, and normal(0, sqrt(c*Sigma2)) if Rho[p]=0\cr   
    Sigma2 ~ scaled-inverse-chi square (Nu, S2)\cr  
    Rho[p] ~ Bernoulli (Kappa)
  
    \item MIX\cr
    Beta[p] ~ normal (0, sqrt(Sigma2[1])) if Rho[p]=1, and normal(0, sqrt(Sigma2[2])) if Rho[p]=0\cr    
    Sigma2[1] ~ scaled-inverse-chi square (Nu, S2)\cr    
    Sigma2[2] ~ scaled-inverse-chi square (Nu, c*S2)\cr    
    Rho[p] ~ Bernoulli (Kappa)
  }
  \bold{Hyperparameters}\cr
  To run vigor, the following hyperparameter values must be declared as arguments for Hyperparameters.\cr
  \itemize{
    \item BL: Phi, Omega
    \item EBL: Phi, Omega, Psi, Theta
    \item wBSR: Nu, S2, Kappa
    \item BayesB: Nu, S2, Kappa     
    \item BayesC: Nu, S2, Kappa  
    \item SSVS: c, Nu, S2, Kappa 
    \item MIX: c, Nu, S2, Kappa
    }
  The argument Hyperparameters is a Nh-length vector, or a (Nh x Nset) matrix, where Nh and Nset are the number of hyperparameters and the number of hyperparameter sets, respectively. For example, when the regression method is BL, the matrix\cr\cr
    1 0.001\cr
    1 0.01\cr
    1 0.1\cr\cr
  indicates that Phi = 1 and Omega = 0.001 in the first set, Phi = 1 and Omega = 0.01 in the second set, and Phi = 1 and Omega = 0.1 in the third set. As another example, when the regression model is BayesB, the vector\cr\cr
    4 0.5 0.001\cr\cr
  indicates that Nu = 4, S2 = 0.5, and Kappa = 0.001. Vectors or matrices for Hyperparameters can be created using hyperpara.\cr\cr
 \bold{Functions}\cr
  The functions of vigor are "fitting", "tuning", and "cv". "Fitting" fits the selected regression model to the data. When Hyperparameters includes multiple hyperparameter sets (i.e., is input as a matrix), only the first set is used.\cr
  "Tuning" selects the hyperparameter set with the lowest MSE in CV. This set is then used in the model fitting. When tuning the hyperparameters, CV is performed on randomly partitioned data. The number of folds is determined by CVFoldTuning.\cr
  The "cv" function conducts CV, and returned the predicted values. When Hyperparameters includes multiple hyperparameter sets, tuning is performed at each fold of the CV.\cr\cr 
\bold{Partition matrix}\cr
  The following is a possible Partition of 20 individuals evaluated in a five-fold CV:\cr\cr
    14 11  3  2  7\cr
     5  4 20 10  9\cr
     6  8 16 15 12\cr
    18 13 17  1 19\cr\cr
  Individuals (row numbers in Pheno/Geno/Covariates) 14, 5, 6, and 18 are removed from the training set at the first fold of the five-fold CV. Samples 11, 4, 8, and 13 are removed at the next fold. This process is repeated up to the fold number of the CV. If the number of individuals N is 19, the gap is filled with -9. For example,\cr\cr
     8  6  3 14 18\cr
    12  4  1 15  5\cr
    17  9 13 11 10\cr
    19 16  7  2 -9\cr\cr
  An example of random sampling validation in which individuals can be sampled more than once is shown below.\cr\cr
    18  3 11 16 13\cr
    17  8 13 13 18\cr
     7 15 14 19  7\cr
     1 13 12  7  2\cr\cr
  Individuals 18, 13, and 7 are repeatedly used as testing samples.\cr\cr
  Random partitioning (i.e., Nfold = n) outputs a Partition matrix, which can be input as the Partition matrix in subsequent analysis.\cr\cr
  \bold{Intercept}\cr
  If Covariates = "Intercept", vigor automatically adds the intercept to the regression model. If Covariates is a user-specified matrix, vigor uses this matrix as the covariates, regarding the first column as the intercept. Thus, the first column of Covariates should be filled with 1s. However, when Covariates is a matrix of population-assignment probabilities when correcting for population stratification, the intercept is not necessary.\cr\cr
  \bold{Standardization}\cr 
  Vigor standardizes Pheno (response variables). Although the returned values are generally scaled back to the original scale, some estimates are reported in the standardized scale. See value.\cr
}
\value{
  When Function = "fitting" or "tuning", a list containing the following elements is returned.
  \item{$LB}{Lower bound of the marginal log likelihood of Pheno.} 
  \item{$ResidualVar}{Residual variances (1/Tau02) of the iterations. Reported in the standardized scale.}   
  \item{$Beta}{Posterior means of the marker effects (E[Beta|Pheno]).\cr
  	The wBSR model returns E[Beta|Pheno]*E[Gamma|Pheno].} 
  \item{$Sd.beta}{Posterior uncertainty (standard deviation) of the marker effects (sqrt(Var[Beta|Pheno])).\cr 
  	For wBSR,\cr 
  	sqrt(E[Beta^2|Pheno]*Var[Gamma|Pheno]+Var[Beta|Pheno]*E[Gamma|Pheno]^2)\cr is returned.}
  \item{$Tau2}{Posterior mean of Tau2.}
  \item{$Sigma2}{Posterior mean of Sigma2.}
  \item{$Alpha}{Posterior means of Alpha.}
  \item{$Sd.alpha}{Posterior uncertainty of Alpha.}
  \item{$Lambda2}{Posterior mean of Lambda2. Reported in the standardized scale. Returned only by the BL model.}  
  \item{$Delta2}{Posterior mean of Delta2. Reported in the standardized scale. Returned only by the EBL model.}   
  \item{$Eta2}{Posterior means of Eta2. Reported in the standardized scale. Returned only by the EBL model.}  
  \item{$Gamma}{Posterior means of Gamma. Returned only by the wBSR model.} 
  \item{$Rho}{Posterior means of Rho. Returned only by the BayesB, BayesC, SSVS, and MIX models.}
  \item{$MSE}{A data frame with (2 + Nh) columns, where Nh is the number of hyperparameters. Returned when Function = "tuning". \cr
    \itemize{
  \item Set: Hyperparameter set number. This number corresponds to the row number of Hyperparameters.
  \item (Hyperparameters): Hyperparameter values of the set.
  \item MSE: The MSE of the hyperparameter set.
  }}
  When Function = "cv", a list containing the following elements is returned.
  \item{$Prediction}{A data frame with 4 columns, Test, Y, Yhat, and BV:
    \itemize{
  \item Test: Tested individuals (row numbers in Pheno/Geno/Covariates).
  \item Y: Phenotypic values of the tested individuals (true values).
  \item Yhat: Predicted values. Summation of the marker and covariate effects.
  \item BV: Breeding values. Summation of marker effects.
  }}
  \item{$MSE}{A data frame with (3 + Nh) columns, returned when analyzing multiple hyperparameter sets.
  	\itemize{  	
  \item Fold: Fold number of CV
  \item ChosenSet: Chosen hyperparameter set at the fold. Set numbers correspond to the row numbers of Hyperparameters.
  \item (Hyperparameters): Hyperparameter values of the chosen set.
  \item MSE: The MSE of the chosen set.
  }}
  \item{$Partition}{A matrix representing the partition used in random partitioning. This matrix can be used as the argument Partition in subsequent analyses.}  
}
\references{
  Onogi & Iwata, VIGoR: variational Bayesian inference for genome-wide regression, in prep.\cr
  Onogi A, 2015, Documents for VIGoR (May 2015).https://github.com/Onogi/VIGoR
}
\author{
  Akio Onogi\cr
  Hiroyoshi Iwata
}
\examples{
#data
data(sampledata)
dim(Geno)#100 samples and 1000 markers
dim(Pheno)#100 samples and 3 traits
dim(Covariates)#100 samples and 2 covariates (including the intercept)

#Use BL. Draw a simple Manhattan plot.
Result<-vigor(Pheno$Height,Geno,"BL",c(1,1),Covariates=Covariates)
plot(abs(Result$Beta),pch=20)
which((abs(Result$Beta)-1.96*Result$Sd.beta)>0) #Significant markers (P<0.05)

#Use BayesC without covariates.
Result<-vigor(Pheno$Height,Geno,"BayesC",matrix(c(4,1,0.01,4,1,0.001),nr=2,byrow=TRUE))
plot(abs(Result$Beta),pch=20)
which((abs(Result$Beta)-1.96*Result$Sd.beta)>0)
print(Result$Alpha)#intercept is automatically added.

#Tuning hyperparameters. Two hyperparameter sets are given as a matrix. Use BayesB.
H<-matrix(c(5,1,0.001,5,1,0.01),nc=3,byrow=TRUE)
print(H)
Result<-vigor(Pheno$Height,Geno,"BayesB",H,Function="tuning",Covariates=Covariates)
plot(abs(Result$Beta),pch=20)
print(Result$MSE)#the set with the lowest MSE was used
#When Function of vigor is "fitting", only the first set is used for regression.
#to repeat analyses under the different sets, for example, 
Result<-as.list(numeric(2))
for(set in 1:2){Result[[set]]<-vigor(Pheno$Height,Geno,"wBSR",H[set,],Covariates=Covariates)}

#Perform cross-validation. Use BL. Number of hyperparameter sets is 2.
#the first set is c(1,0.01) and the second is c(1,0.1)
#6-fold CV
Result<-vigor(Pheno$Height,Geno,"BL",
matrix(c(1,0.01,1,0.1),ncol=2,byrow=TRUE),Function="cv",Nfold=6,Covariates=Covariates)
plot(Result$Prediction$Y,Result$Prediction$Yhat) #plot true and predicted values
cor(Result$Prediction$Y,Result$Prediction$Yhat) #accuracy
print(Result$MSE) #see which the set used at each fold.
print(Result$Partition) #see the partition of CV
#Perform CV using the same partition. Use BayesC.
H<-matrix(c(5,1,0.01,5,1,0.1),nc=3,byrow=TRUE)
Result2<-vigor(Pheno$Height,Geno,"BayesC",H,Function="cv",Nfold=-9,
Partition=Result$Partition,Covariates=Covariates) 
cor(Result2$Prediction$Y,Result2$Prediction$Yhat) #accuracy
}
