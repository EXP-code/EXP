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

#define   CONFLUENT 1
#define NOCONFLUENT 0

extern GaussMaster();
extern GaussCheck();

/* ARGUMENT LISTS:
/* 
/* GaussCheck(value)
/* double value;
/* 
/* GaussMaster(n, alpha, beta, conflag, abscis, weight)
/* int n;
/* double alpha, beta;
/* int conflag;
/* double abscis[], weight[];
/* 
/**/
