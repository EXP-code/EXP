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

/* Translation of module Hermite (import, implement).
/**/

#include <math.h>		/* square root */
#include "GaussCore.h"
#include "Laguerre.h"
#include "Hermite.h"

/*
/* Macro for Pascal's odd().
/**/
#define odd(n) ((unsigned)(n) & 01)

Odd_Hermite(n, alpha, abscis, weight, w0)
int n;
double alpha, abscis[], weight[], *w0;
{
    int k;

    GaussCheck(alpha);
    n /= 2;
    alpha = (alpha-1)/2;
    Radau_Laguerre(n, alpha, abscis, weight, w0);

    for (k=0; k<n; k++) {
	weight[k] /= 2;
	abscis[k] = sqrt(abscis[k]);
    };
};


Even_Hermite(n, alpha, abscis, weight)
int n;
double alpha, abscis[], weight[];
{
    int k;

    GaussCheck(alpha);
    n /= 2;
    alpha= (alpha-1)/2;
    Laguerre(n, alpha, abscis, weight);

    for (k=0; k<n; k++) {
	weight[k] /= 2;
	abscis[k] = sqrt(abscis[k]);
    };
};


Hermite(n, alpha, abscis, weight)
int n;
double alpha, abscis[], weight[];
{
    int n1, k;

    n1= n / 2;

    if ( ! odd(n) )
        Even_Hermite(n, alpha, &abscis[n1], &weight[n1]);
    else {
	Odd_Hermite(n, alpha, &abscis[1+n1], &weight[1+n1], &weight[n1]);
	abscis[n1]= 0;
    };

    for (k=0; k<n1; k++) {
	abscis[k] = -abscis[n-k-1];
	weight[k] =  weight[n-k-1];
    };
}
