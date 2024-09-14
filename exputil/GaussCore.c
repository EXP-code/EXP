/* Copyright 1985 Pierre Asselin.

   This package may be used and distributed freely, but may not be sold.
   The present notice must be retained in all copies.
   If incorporated in a commercial product,

    1) the package's origin and availability must be acknowledged
       prominently in the commercial product's documentation;

    2) the source code and documentation of the package must be
       made available to purchasers on request and at no extra cost.
*/

/* 
   Translation of module GaussKernel (implement).
*/

#include <stdio.h>		/* for error messages */
#include <stdlib.h>
#include "GaussCore.h"

#define GaussEPS 1.0e-12
#define abs(x) ( ((x) < 0) ? -(x) : (x) )
#define sqr(x) ((x)*(x))

/*
   Forward declaration for function QQp:
*/
static int QQp(double, double*, double*);

/*
   Function to test a real value for exceeding -1,
   and quit on an error if condition not met.
*/
void GaussCheck(double value)
{
    if (value <= (-1.0)) {
	fprintf(stderr, "Gauss package: parameter out of range: %g\n", value);
	exit(1);
    }
}



/*
   Static variables needed to imitate Pascal's block-structured
   scope rules.
*/
static double alf1, bet1;
static int confl1;
static int n1;

/*
   The workhorse.
*/
void GaussMaster(int n, double alpha, double beta, int conflag, double abscis[], double weight[])
{
#define FALSE 0
#define  TRUE 1
    typedef int bool;
    int k, m;
    int below;
    double t, min, max, Glob;
    double delta, Qp;
    double temp;
    bool ok;
#define junk1 &delta
#define junk2 &Qp
    
/* Construction to share storage.
   Check if porting to another machine.
*/
    typedef union {
	double d;
	int i;
    } workcell;

    workcell *rank = (workcell *) weight;     /* rank[k].i is on weight[k] */
    double *upbnd = abscis;		      /* upbnd[k]  is on abscis[k] */


  /*
     Give a copy of formal parameters to QQp.
  */
    alf1= alpha;  bet1= beta;  confl1= conflag;  n1= n;

  /*
     Begin first stage: find upper bound to all roots.
     Do not waste the intermediate information generated along the way.
     upbnd[k]= best known upper bound to root #(k+1),
     rank[k].i= # roots below upbnd[k] ( k+1 or less ).
  */
    for (k=0; k<n; k++) rank[k].i= 0;

    t= 1.0;			/* turns out to be a good choice */
    m= 0;
    do {
	t *= 2;				/* double t */
	below= QQp(t, junk1, junk2);	/* rank it */

	for (k=m; k<below; k++) {	/* mark known bounds */
	    rank[k].i= below;		/* of roots less than t */
	    upbnd[k]= t;
	};
	m= below;
    } while (below < n);	/* stops when t larger than all roots */


  /*
     Begin second stage: Isolate the roots by bisection.
     Crunch until (rank[k].i == (k+1)) for all k, at which point
     root #(k+1) lies between upbnd[k-1] and upbdn[k] for all k.
     (define upbnd[-1] as 0.0).
  */
    min= 0.0;			/* lower bound of lowest root */

    for (k=0; k<n; k++) {
	while ((rank[k].i) > (k+1) ) {

	  /* bisect, test, update. */
 	    t= (min + upbnd[k])/2;
	    below= QQp(t, junk1, junk2);
	    if (below <= k) min= t;
	    for (m=k; m<below; m++) {
	        upbnd[m]= t;
	        rank[m].i = below;
	    };
	};
	min= upbnd[k];  	 /* prepare l. b. of next root */
    };


  /*
     Before Newton refinement:
     compute the global factor Glob.
  */
    Glob= (alpha + 1.0);
    if (conflag)
        for (k=2; k<= n; k++) Glob *= ((k+alpha)/k);
    else {
	Glob *= (beta + 1.0);
	for (k=2; k <= n; k++)
	    Glob *= (((k+alpha)*(k+beta))/(k*(k+alpha+beta)));
    };

  /*
     Begin third stage: Newton refinement.
     Cautiously refine each root to accuracy.
       Newton's method needs a good starting approximation;
       revert to bisection if things look bad.
       The construct "do {... break...break...} while false"
       is used to attempt a single Newton step and abort at the
       first sign of trouble.
  */
    min= 0;				/* lower bound of lowest root */

    for (k=0; k<n; k++) {
	max= upbnd[k];			/* upper bound of current root */
	t= (min + max)/2;		/* starting approximation */
	do {
	    below= QQp(t, &delta, &Qp);	/* Compute values, rank */

	  /* Squeeze the bounds now.  Bug discovered by
	     todd@cougarxp.Princeton.edu (Todd Jay Mitty).
	  */ 
	    if (below <= k) min= t; else max= t;

	  /* One cautious Newton step.
	     -Desirable improvement:
	      include machine roundoff in equality tests.
	  */
	    ok= FALSE;
	    do {
		if (Qp==0)  break;
		delta/= -Qp;
		if (conflag) {
		    temp= t + delta;
		} else {
		    if (delta == 1)  break;
		    temp= (t + delta)/(1 - delta);
		    delta*= (1+t);
		}
		ok = (min <= temp) && (temp <= max);
	    } while (FALSE);

	    if (ok) t=temp;		/* use Newton iterate... */
	    else t= (min + max)/2;	/* ...or bisect to recover */

	} while (!ok || (abs(delta) >= abs(t)*GaussEPS) );

	min= upbnd[k];     /*  prepare lower bound for next root */

	(void) QQp(t, junk1, &Qp);	/* derivative needed for weight */
	abscis[k]= t;			/* clobbers upbnd[k]  */
	weight[k]= Glob/(t * sqr(Qp));	/* clobbers rank[k].i */
    }
#undef junk1
#undef junk2
#undef  TRUE
#undef FALSE
}


/*   The polynomial computer...
     alf1, bet1, n1, confl1 are passed by GaussMaster.
*/
static int QQp(double t, double* poly, double* deriv)
{
    int k;
    double tmp, tmpa, tmpb, tmpab;
    int negative;
    int rank;

    *poly= 1.0;  *deriv= 0.0;  rank=0;

    for (k=1; k<=n1; k++) {
	tmpa= k + alf1;
        negative= (*poly < 0);

	if (confl1) {
	    *deriv -= *poly;
	    *poly = ( (*poly * tmpa) +  (*deriv * t) )/k;
	}
	else {
	    tmpb= k + bet1;
	    tmpab= tmpa + tmpb;

	    tmp= *deriv/(tmpa + bet1);
	    *deriv= (tmpb - t*tmpa)*tmp -  *poly * tmpab;
	    *poly=  *poly * (tmpa - t*tmpb) + t*tmpab*tmp;

	    tmp= t+1;
	    *deriv /= tmp;
	    *poly /= (k * tmp);
	};

	if (negative != (*poly < 0)) ++rank;
    };
    return rank;
}
