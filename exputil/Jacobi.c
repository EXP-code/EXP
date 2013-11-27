/* Copyright 1985 Pierre Asselin.
/*
/* This package may be used and distributed freely, but may not be sold.
/* The present notice must be retained in all copies.
/* If incorporated in a commercial product,
/*  1) the package's origin and availability must be acknowledged
/*     prominently in the commercial product's documentation;
/*  2) the source code and documentation of the package must be
/*     made available to purchasers on request and at no extra cost.
/**/

/* Translation of module Jacobi (import, implement).
/**/

#include "GaussCore.h"
#include "Jacobi.h"

#define sqr(x) ((x)*(x))

Jacobi(n, alpha, beta, abscis, weight)
int n;
double alpha, beta;
double abscis[], weight[];
{
    int k;

    GaussCheck(alpha);
    GaussCheck(beta);
    GaussMaster(n, alpha, beta, NOCONFLUENT, abscis, weight);
    for (k=0; k<n; k++) abscis[k] /= (1 + abscis[k]);
}


Radau_Jacobi(n, alpha, beta, abscis, weight, leftw)
int n;
double alpha, beta;
double abscis[], weight[], *leftw;
{
    int k;
    double temp;

    GaussCheck(alpha);
    GaussCheck(beta);
    temp= 1 + alpha;
    Jacobi(n, temp, beta, abscis, weight);
    temp /= (1 + temp + beta);

    for (k=0; k<n; k++) weight[k] *= (temp/abscis[k]);

    if (n == 0) *leftw=1;
    else {
	temp= 1 + n + alpha;
	temp *= (temp + beta);
	temp = (1 + beta)/temp;
	for (k=2; k<=n; k++)
	    temp *= ((k * (k+beta)) / ((k+alpha) * (k+alpha+beta)));
	*leftw= temp;
    };
};


Lobatto_Jacobi(n, alpha, beta, abscis, weight, leftw, rightw)
int n;
double alpha, beta;
double abscis[], weight[], *leftw, *rightw;
{
    int k;
    double temp1, temp2;

    GaussCheck(alpha);
    GaussCheck(beta);
    alpha += 1;  beta += 1;
    GaussMaster(n, alpha, beta, NOCONFLUENT, abscis, weight);

    temp1= (alpha * beta)/((alpha+beta) * (1+alpha+beta));
    for (k=0; k<n; k++) {
	temp2= 1 + abscis[k];
	weight[k] *= ( temp1 * sqr(temp2) ) / abscis[k];
	abscis[k] /= temp2;
    };

    temp1= beta/(alpha + beta);
    for (k=1; k <= n; k++)
        temp1 *= ((k * (k+beta))/((k+alpha) * (k+alpha+beta)));
    *leftw= temp1;

    temp1= alpha/(alpha + beta);
    for (k=1; k <= n; k++)
        temp1 *= ((k * (k+alpha))/((k+beta) * (k+beta+alpha)));
    *rightw= temp1;
}
