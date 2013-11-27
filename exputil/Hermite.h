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

/* Translation of module Hermite (export).
/**/

extern Odd_Hermite();
extern Even_Hermite();
extern Hermite();

/* ARGUMENT LISTS:
/* 
/* Odd_Hermite(n, alpha, abscis, weight, w0)
/* int n;
/* double alpha, abscis[], weight[], *w0;
/* 
/* Even_Hermite(n, alpha, abscis, weight)
/* int n;
/* double alpha, abscis[], weight[];
/* 
/* Hermite(n, alpha, abscis, weight)
/* int n;
/* double alpha, abscis[], weight[];
/* 
/**/
