
/*
        QADAPT - a recursive adaptive quadrature program,
                callable from C.

        Calling syntax:

                qadapt(a, b, func, tolerance)

        Arguments: 

                double a  --- lower limit of integration
                double b  --- upper limit of integration
                func      --- the function you want to 
                              integrate
                double tolerance --- error tolerance


        Returns: 
                double 






        QADAPT is an adaptive quadrature routine based upon the
        midpoint rule. The heart of the algorithm is the 
        function quadsplit(), which takes an interval and divides
        it into three parts, and then uses the midpoint rule to estimate 
        the integral over the whole interval and over each subinterval.
        If these two estimates differ significantly, quadsplit() is called
        for each subinterval in turn.

        If the number of subdivisions becomes excessively large
        (currently the maximum is about 25 -- a dynamic range of 10^12),
        qadapt() will flash a warning signal and then return its
        most recent guess at the integral.
        If this occurs, first make sure the interval is supposed
        to converge at all! If your integral is convergent
        then you may be able to fix the trouble by either (1) relaxing
        your error tolerance or (2) subdividing the interval ``by hand''
        and calling qadapt() on each piece.

        Remember that if your function has a particularly narrow spike
        in an arbitrary place, qadapt() might not find it on the
        first iterations and hence might think it has converged.
        If you are faced with a function having spikes at places
        unknown to you, be sure to call qadapt() on intervals small
        enough to give it a chance to detect the spikes.



                Author: Kevin Long 
                Written: 2/25/91
                Revised: 3/11/91
                ANSI Version: 9/30/91

*/




#include <math.h>
#include <stdio.h>
#include <numerical.H>








/* 
	some numerical constants that will be needed
*/
#define ONE_SIXTH 	0.16666666666667
#define FIVE_SIXTHS 	0.83333333333333

/* this choice of TINY allows about 25 subdivisions */


#ifndef TINY
#define TINY 1.e-16
#endif


/* 
	for compilers with small default stack, increase stack size
	to prevent recursion from overflowing stack.
*/
#ifdef __TURBOC__
	extern int _stklen = 32767;
#endif


/* 
	a global variable holding the initial interval width
*/

static double qadapt_initial_width;







/*
	top level function --- takes limits a and b, a function func,
	and an error tolerance tol. It evaluates the function
	at the midpoint and then subdivides the interval.
	The subdivision function then calls itself as many times 
	as necessary to achieve error<tol.
*/

double quadsplit(double, double, std::function<double(double)>, double, double);


double qadapt(double a, double b, std::function<double(double)> func, double tol)
{
	double fmid;

	/* evaluate function at midpoint as first guess for integral */
	fmid = func((a+b)/2.0);

	/* store initial width as global variable */
	qadapt_initial_width = b-a;

	/* call subdivision routine */
	return (b-a)*quadsplit(a, b, func, fmid, tol);
}




double quadsplit(double a, double b, std::function<double(double)> func,
		 double fmid, double tol)
{
	double flow, fhigh, fmean, delta;
	double qlow, qmid, qhigh;


	/* evaluate the function at the low and high sample points.
	No need to do the midpoint again, since it has been passed along
	from the calling program */
	
	flow = func(a + ONE_SIXTH*(b-a));
	fhigh = func(a + FIVE_SIXTHS*(b-a));
	

	/* compute average of low, mid, and high function evaluations */
	
	fmean = (flow + fmid + fhigh)/3.0;
	

	/* 
	see if splitting the interval made any difference; if it did
	not, call it a day and return to previous level.
	*/
	
	if (fabs(b-a)*fabs(fmean - fmid) < tol) 
	{
		return fmean;
	}
		

	/* 
	the average and midpoint values did not compare well, so we
	have to look at the interval more closely.
	*/
	
	delta = (b-a)/3.0;

	/* 
	if the interval size gets too small, give up.
	This will happen if the integral does not converge,
	and in some cases if it is converging too slowly.
	*/

	if (fabs(delta/qadapt_initial_width) <= TINY) 
	{
		puts("I'm giving up: too many subdivisions");
		return fmean;
	}

	/* recursively split all three intervals */

	qlow = quadsplit(a, a + delta, func, flow, tol);
	qmid = quadsplit(a+delta, a+2.0*delta, func, fmid, tol);
	qhigh = quadsplit(a+2.0*delta, b, func, fhigh, tol);

	/* return the average of the three integrals */

	return (qlow + qmid + qhigh)/3.0;
}
	
	

#undef TINY	


	


