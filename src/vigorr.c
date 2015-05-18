#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "MT.h"
#include "vigorr.h"

/*
	VIGoR function for R package.

	Copyright (C) 2015 Akio Onogi and Hiroyoshi Iwata

	Released under the MIT license
	http://opensource.org/licenses/mit-license.php

*/

void vigorr (int *Priortype, int *Methodcode, int *CondResidual, int *P, int *F, int *N, int *Maxiteration, int *RandomIni, double *Thresholdvalue,
	double *Observations, double *Genotypes, double *Covariates, double *Hyperparameters, double *Tau0, double *LBmonitor, double *Rmonitor,
	double *Expectations, double *Uncertainty, double *Variance, double *ExpectationsQ, double *UncertaintyQ, double *Gamma, double *expDelta2, double *Eta2)
{
	/* Observations, genotypes, covariates, and hyperparameters */
	Ystruct *Y;
	Xstruct *X, *Q;
	Hstruct *H;

	/* repeat statement */
	int	locus;

	/* Calculation time */
	double	FittingTime=0.0;

	/* Printinfo */
	int		Printinfo=0;

	/* Temporary use */
	double temp;

	/* used seed */
	double usedseed;

	/* Use variational Bayesian inference */
	int	Algorithm=1;

	/* Copy */
	Y = (Ystruct*) calloc (1, sizeof(Ystruct));
	Y[0].observations = (double*) calloc (N[0], sizeof(double));
	memcpy(Y[0].observations, Observations, sizeof(double)*N[0]);
	X = (Xstruct*) calloc ( P[0], sizeof(Xstruct));
	for(locus=0; locus<P[0]; locus++)
	{
		X[locus].covariates = (double*)calloc(N[0], sizeof(double));
		memcpy(X[locus].covariates, Genotypes+locus*N[0], sizeof(double)*N[0]);
	}
	Q = (Xstruct*) calloc ( F[0], sizeof(Xstruct));
	for(locus=0; locus<F[0]; locus++)
	{
		Q[locus].covariates = (double*)calloc(N[0], sizeof(double));
		memcpy(Q[locus].covariates, Covariates+locus*N[0], sizeof(double)*N[0]);
	}
	H = (Hstruct*) calloc (1, sizeof(Hstruct));
	switch (Methodcode[0]){
		case 1: H[0].deltaShape = Hyperparameters[0]; H[0].deltaRate = Hyperparameters[1]; break;
		case 2: H[0].deltaShape = Hyperparameters[0]; H[0].deltaRate = Hyperparameters[1]; H[0].etaShape = Hyperparameters[2]; H[0].etaRate = Hyperparameters[3]; break;
		case 3: H[0].v = Hyperparameters[0]; H[0].S2 = Hyperparameters[1]; H[0].Pi = Hyperparameters[2]; break;
		case 4: H[0].v = Hyperparameters[0]; H[0].S2 = Hyperparameters[1]; H[0].Pi = Hyperparameters[2]; break;
		case 5: H[0].c = Hyperparameters[0]; H[0].v = Hyperparameters[1]; H[0].S2 = Hyperparameters[2]; H[0].Pi = Hyperparameters[3]; break;
		case 6: H[0].c = Hyperparameters[0]; H[0].v = Hyperparameters[1]; H[0].S2 = Hyperparameters[2]; H[0].Pi = Hyperparameters[3]; break;
		case 7: H[0].v = Hyperparameters[0]; H[0].S2 = Hyperparameters[1]; H[0].Pi = Hyperparameters[2]; break;
	}

	/* Initialize the seed */
	usedseed = (int)time(NULL);
	init_genrand(usedseed);

	/* Calculation */
	GenomeWideRegression (Algorithm, Priortype[0], Methodcode[0], CondResidual[0], P[0], F[0], N[0], Thresholdvalue[0], Maxiteration[0], 
	Y, X, Q, H, expDelta2, Tau0, LBmonitor, Rmonitor, &FittingTime, Printinfo, RandomIni[0]);

	/* Copy */
	for(locus=0; locus<F[0]; locus++){ ExpectationsQ[locus]=Q[locus].expEffect;}
	switch(Methodcode[0]){
	case 1:
		for(locus=0; locus<P[0]; locus++)
		{
			Expectations[locus] = X[locus].expEffect;
			Variance[locus] = X[locus].expTau2;
		}
		break;
	case 2:
		for(locus=0; locus<P[0]; locus++)
		{
			Expectations[locus] = X[locus].expEffect;
			Variance[locus] = X[locus].expTau2;
			Eta2[locus] = X[locus].expEta2;
		}
		break;
	case 3:
		for(locus=0; locus<P[0]; locus++)
		{
			Expectations[locus] = X[locus].expEffect*X[locus].expGamma;
			Variance[locus] = X[locus].expSigma2;
			Gamma[locus] = X[locus].expGamma;
		}
		break;
	case 4:
		for(locus=0; locus<P[0]; locus++)
		{
			Expectations[locus] = X[locus].expEffect*X[locus].expGamma;
			Gamma[locus] = X[locus].expGamma;
		}
		Variance[0] = X[0].expSigma2;
		break;
	case 5:
		for(locus=0; locus<P[0]; locus++)
		{
			Expectations[locus] = X[locus].expEffect;
			Gamma[locus] = X[locus].expGamma;
		}
		Variance[0] = X[0].expSigma2;
		break;
	case 6:
		for(locus=0; locus<P[0]; locus++)
		{
			Expectations[locus] = X[locus].expEffect;
			Gamma[locus] = X[locus].expGamma;
		}
		Variance[0] = X[0].expSigma2;
		Variance[1] = X[1].expSigma2;
	case 7:
		for(locus=0; locus<P[0]; locus++)
		{
			Expectations[locus] = X[locus].expEffect*X[locus].expGamma;
			Variance[locus] = X[locus].expSigma2;
			Gamma[locus] = X[locus].expGamma;
		}
		break;
		break;
	}

	if(Algorithm==1)
	{	
		for(locus=0; locus<F[0]; locus++) {UncertaintyQ[locus]=sqrt(Q[locus].varEffect);}
		switch(Methodcode[0]){
		case 1: for(locus=0; locus<P[0]; locus++) Uncertainty[locus] = sqrt(X[locus].varEffect); break;
		case 2: for(locus=0; locus<P[0]; locus++) Uncertainty[locus] = sqrt(X[locus].varEffect); break;
		case 3: for(locus=0; locus<P[0]; locus++)
				{
					temp = X[locus].expGamma * (1.0-X[locus].expGamma) * X[locus].varEffect 
							+ X[locus].expGamma * (1.0-X[locus].expGamma) * pow(X[locus].expEffect, 2.0) 
							+ pow(X[locus].expGamma,2.0) * X[locus].varEffect;
					Uncertainty[locus] = sqrt(temp);
				}
				break;
		case 4: for(locus=0; locus<P[0]; locus++)
				{
					temp = pow(X[locus].expEffect,2.0)*X[locus].expGamma*(1.0-X[locus].expGamma)+X[locus].varEffect*X[locus].expGamma;
					Uncertainty[locus] = sqrt(temp);
				}
				break;
		case 5: for(locus=0; locus<P[0]; locus++) Uncertainty[locus] = sqrt(X[locus].varEffect); break;
		case 6: for(locus=0; locus<P[0]; locus++) Uncertainty[locus] = sqrt(X[locus].varEffect); break;
		case 7: for(locus=0; locus<P[0]; locus++)
				{
					temp = pow(X[locus].expEffect,2.0)*X[locus].expGamma*(1.0-X[locus].expGamma)+X[locus].varEffect*X[locus].expGamma;
					Uncertainty[locus] = sqrt(temp);
				}
				break;
		}	
	}

	/* free */
	free (Y[0].observations); free(Y);
	for(locus=0; locus<P[0]; locus++) free(X[locus].covariates); free(X);
	for(locus=0; locus<F[0]; locus++) free(Q[locus].covariates); free(Q);
	free(H);

}
