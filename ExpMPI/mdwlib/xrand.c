#include <math.h>

/* Random number generators:
 *
 *  rnd_init (unsigned seed) 
 *			: initializes the generator
 *
 *  rnd_i ()		: returns positive integers [0,0x7fffffff]
 *  rnd_u ()		: returns unsigned's        [0,0xffffffff]
 *  rnd_ri (long n)	: returns positive integers [0,n-1]
 *  rnd_01d ()		: returns doubles	    [0.0,1.0)
 *			  Note: ")" is no typo - rnd_01d will not return a 1.0,
 *                              but can return the next smaller FP number.
 *  rnd_ned (double lam): returns neg. exponential distributed doubles [0.0,+inf)
 *
 *  Algorithm M as describes in Knuth's "Art of Computer Programming", Vol 2. 1969
 *  is used with a linear congruential generator (to get a good uniform
 *  distribution) that is permuted with a Fibonacci additive congruential
 *  generator to get good independence.
 *
 *  Bit, byte, and word distributions were extensively tested and pass
 *  Chi-squared test near perfect scores (>7E8 numbers tested, Uniformity
 *  assumption holds with probability > 0.999)
 *
 *  Passed Kolmogorov-Smirnov tests of 1000 trials as suggested by Knuth
 *
 *  Run-up tests for on 7E8 numbers confirm independence with
 *  probability > 0.97.
 *
 *  Plotting random points in 2d reveals no apparent structure.
 *
 *  Autocorrelation on sequences of 5E5 numbers (A(i) = SUM X(n)*X(n-i), i=1..512)
 *  results in no obvious structure (A(i) ~ const).
 *
 *  On a SUN 3/60, rnd_u() takes about 19.4 usec per call, which is about 44%
 *  slower than Berkeley's random() (13.5 usec/call).
 *
 *  Except for speed and memory requirements, this generator outperforms
 *  random() for all tests. (random() scored rather low on uniformity tests,
 *  while independence test differences were less dramatic).
 *
 *  Adopted from a routine written by A. Nowatzyk at CMU
 *  by MDW
 *
 */

/* LC-parameter selection follows recommendations in 
 * "Handbook of Mathematical Functions" by Abramowitz & Stegun 10th, edi.
 */
#define LC_A 66049		    /* = 251^2, ~= sqrt(2^32)			*/
#define LC_C 3907864577		    /* result of a long trial & error series    */

#define Xrnd(x) (x * LC_A + LC_C)   /* the LC polynomial			*/
			
static unsigned long Fib[55];	    /* will use X(n) = X(n-55) - X(n-24)	*/
static int Fib_ind;		    /* current index in circular buffer		*/
static unsigned long Xrnd_var;	    /* LCA - recurrence variable		*/
static unsigned long auxtab[256];   /* temporal permutation table		*/
static unsigned long prmtab[64] = { /* spatial permutation table		*/
    0xffffffff, 0x00000000,  0x00000000,  0x00000000,  /* 3210 */
    0x0000ffff, 0x00ff0000,  0x00000000,  0xff000000,  /* 2310 */
    0xff0000ff, 0x0000ff00,  0x00000000,  0x00ff0000,  /* 3120 */
    0x00ff00ff, 0x00000000,  0xff00ff00,  0x00000000,  /* 1230 */

    0xffff0000, 0x000000ff,  0x00000000,  0x0000ff00,  /* 3201 */
    0x00000000, 0x00ff00ff,  0x00000000,  0xff00ff00,  /* 2301 */
    0xff000000, 0x00000000,  0x000000ff,  0x00ffff00,  /* 3102 */
    0x00000000, 0x00000000,  0x00000000,  0xffffffff,  /* 2103 */

    0xff00ff00, 0x00000000,  0x00ff00ff,  0x00000000,  /* 3012 */
    0x0000ff00, 0x00000000,  0x00ff0000,  0xff0000ff,  /* 2013 */
    0x00000000, 0x00000000,  0xffffffff,  0x00000000,  /* 1032 */
    0x00000000, 0x0000ff00,  0xffff0000,  0x000000ff,  /* 1023 */

    0x00000000, 0xffffffff,  0x00000000,  0x00000000,  /* 0321 */
    0x00ffff00, 0xff000000,  0x00000000,  0x000000ff,  /* 0213 */
    0x00000000, 0xff000000,  0x0000ffff,  0x00ff0000,  /* 0132 */
    0x00000000, 0xff00ff00,  0x00000000,  0x00ff00ff   /* 0123 */
};

union hack {			    /* used to access doubles as unsigneds	*/
    double d;
    unsigned long u[2];
};

static union hack man;		    /* mantissa bit vector			*/

rnd_init (unsigned int seed)			    /* modified: seed 0-31 use precomputed stuff */
                  
{
    register unsigned long u;
    register int i;
    double x, y;
    union hack t;

    static unsigned seed_tab[32] = {
		0xbdcc47e5, 0x54aea45d, 0xec0df859, 0xda84637b,
		0xc8c6cb4f, 0x35574b01, 0x28260b7d, 0x0d07fdbf,
		0x9faaeeb0, 0x613dd169, 0x5ce2d818, 0x85b9e706,
		0xab2469db, 0xda02b0dc, 0x45c60d6e, 0xffe49d10,
		0x7224fea3, 0xf9684fc9, 0xfc7ee074, 0x326ce92a,
		0x366d13b5, 0x17aaa731, 0xeb83a675, 0x7781cb32,
		0x4ec7c92d, 0x7f187521, 0x2cf346b4, 0xad13310f,
		0xb89cff2b, 0x12164de1, 0xa865168d, 0x32b56cdf  };

    if (seed < 32)
	u = seed_tab[seed];
    else
	u = seed ^ seed_tab[seed & 31];

    for (i = 55; i--;)		    /* set up Fibonacci additive congruential	*/
	Fib[i] = u = Xrnd(u);

    for (i = 256; i--;)
	auxtab[i] = u = Xrnd(u);

    Fib_ind = u % 55;		    /* select a starting point			*/

    Xrnd_var = u;

    if (sizeof(x) != 2 * sizeof(unsigned long)) {
	x = 0.0;
	y = 1.0;
	y /= x;			    /*** intentional divide by 0: rnd_01d will
					 not work because a double doesn't fit
					 in 2 unsigned longs on your machine! ***/
	exit(-101);
    };

    x = 1.0;
    y = 0.5;
    do {			    /* find largest fp-number < 2.0 */
	t.d = x;
	x += y;
	y *= 0.5;
    } while (x != t.d && x < 2.0);  /* one or the other should be true . . .
				       if not, the RNG will fail miserably
				       for doubles */

    man.d = 1.0;
    man.u[0] ^= t.u[0];
    man.u[1] ^= t.u[1];		    /* man is now 1 for each mantissa bit	*/
}

long rnd_i (void)
/*
 * returns a positive, uniformly distributed random number in [0,0x7fffffff]
 */
{ 
    register unsigned long i, j, *t = Fib;

    i = Fib_ind;
    j = t[i];				    /* = X(n-55) */
    j -= (i >= 24) ? t[i - 24] : t[i + 21]; /* = X(n-24) */
    t[i] = j;
    if (++i >= 55) i = 0;
    Fib_ind = i;

    t = &auxtab[(j >> 24) & 0xff];
    i = *t;
    Xrnd_var = *t = Xrnd(Xrnd_var);
    t = &prmtab[j & 0x3c];

    j =  *t++ & i;
    j |= *t++ & ((i << 24) | ((i >>  8) & 0x00ffffff));
    j |= *t++ & ((i << 16) | ((i >> 16) & 0x0000ffff));
    j |= *t   & ((i <<  8) | ((i >> 24) & 0x000000ff));
    
    return j & 0x7fffff;
}

unsigned long rnd_u (void)
/*
 * same as rnd_i, but gives full 32 bit range
 */
{ 
    register unsigned long i, j, *t = Fib;

    i = Fib_ind;
    j = t[i];				    /* = X(n-55) */
    j -= (i >= 24) ? t[i - 24] : t[i + 21]; /* = X(n-24) */
    t[i] = j;
    if (++i >= 55) i = 0;
    Fib_ind = i;

    t = &auxtab[(j >> 24) & 0xff];
    i = *t;
    Xrnd_var = *t = Xrnd(Xrnd_var);
    t = &prmtab[j & 0x3c];

    j =  *t++ & i;
    j |= *t++ & ((i << 24) | ((i >>  8) & 0x00ffffff));
    j |= *t++ & ((i << 16) | ((i >> 16) & 0x0000ffff));
    j |= *t   & ((i <<  8) | ((i >> 24) & 0x000000ff));
    
    return j;
}

long rnd_ri (long int rng)
             
/*
 * randint: Return a random integer in a given Range [0..rng-1]
 *          Note:  0 < rng
 */
{
    register unsigned long  r, a;

    do {
	r = rnd_i();
	a = (r / rng) + 1;
	a *= rng;
    } while (a >= 0x7ffffff);
    
    a--;
    return a - r;
}

double rnd_01d (void)
/*
 * returns a uniformly distributed double in the range of [0..1)
 *         or  0.0 <= rnd_01d() < 1.0 to be precise
 *
 * Note: this code assumes that 2 'unsigned long's can hold a 'double'
 *       (works on SUN-3's, SUN-4's, MIPS, VAXen, IBM RT's)
 */
{
    union hack t;

    t.d = 1.0;

    t.u[0] |= rnd_u() & man.u[0];	      /* munch in 1st part   */
    t.u[1] |= rnd_u() & man.u[1];	      /* munch in 2nd part   */

    return t.d - 1.0;
}

double rnd_ned (double lam)
               
/*
 * returns a neg. exponential distributed double in the range of [0..+infinity)
 *         or  0.0 <= rnd_neg() < +infinity to be precise
 *
 * Note: this code assumes that 2 'unsigned long's can hold a 'double'
 *       it also assumes that 'log()' behaves as advertised.
 *
 */
{
    union hack t;

    t.d = 1.0;

    t.u[0] |= rnd_u() & man.u[0];	      /* munch in 1st part   */
    t.u[1] |= rnd_u() & man.u[1];	      /* munch in 2nd part   */

    return -log(2.0 - t.d) / lam;
}
