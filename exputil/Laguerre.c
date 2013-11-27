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

/* Translation of module Laguerre (import, implement).
/**/

#include "GaussCore.h"
#include "Laguerre.h"

Laguerre(n, alpha, abscis, weight)
int n;
double alpha, abscis[], weight[];
{
    GaussCheck(alpha);
    GaussMaster(n, alpha, 0.0, CONFLUENT, abscis, weight);
};


Radau_Laguerre(n, alpha, abscis, weight, leftw)
int n;
double alpha, abscis[], weight[], *leftw;
{
    int k;
    double temp;

    GaussCheck(alpha);
    temp= 1 + alpha;
    Laguerre(n, temp, abscis, weight);
    for (k=0; k<n; k++) weight[k] *= (temp/abscis[k]);

    temp= 1;
    for (k=1; k<=n; k++) temp *= (k / (1+k+alpha));
    *leftw= temp;
};
