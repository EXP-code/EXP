#include <Random.h>
#include <Gamma.h>


double Gamma::getNormal()
{
  if (haveCachedNormal == 1) {
    haveCachedNormal = 0;
    return(cachedNormal);
  } else {
        
    for(;;) {
      double u1 = pGenerator -> asDouble();
      double u2 = pGenerator -> asDouble();
      double v1 = 2 * u1 - 1;
      double v2 = 2 * u2 - 1;
      double w = (v1 * v1) + (v2 * v2);
      
//
//      We actually generate two IID normal distribution variables.
//      We cache the one & return the other.
// 
      if (w <= 1) {
	double y = sqrt( (-2 * log(w)) / w);
	double x1 = v1 * y;
	double x2 = v2 * y;
               
	haveCachedNormal = 1;
	cachedNormal = x2;
	return(x1);
      }
    }
  }
}

static const float q1 =  0.04166669;
static const float a3 =  0.2000062;
static const float a4 =  -0.1662921;
static const float a5 =  0.1423657;
static const float a6 =  -0.1367177;
static const float a7 =  0.1233795;
static const float e1 =  1.0;
static const float e2 =  0.4999897;
static const float e3 =  0.166829;
static const float e4 =  0.0407753;
static const float e5 =  0.010293;
static const float q2 =  0.02083148;
static const float sqrt32 =  5.656854;
static const float q3 =  0.00801191;
static const float q4 =  0.00144121;
static const float q5 =  -7.388e-5;
static const float q6 =  2.4511e-4;
static const float q7 =  2.424e-4;
static const float a1 =  0.3333333;
static const float a2 =  -0.250003;

/*
  From:
     ALGORITHM 599, COLLECTED ALGORITHMS FROM ACM.
     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.9, NO. 2,
     JUN., 1983, P. 255-257.
*/
double Gamma::operator()()
{
  /* Initialized data */


  static float aa = (float)0.0;
  static float aaa = (float)0.0;

  /* Local variables */
  float ret_val, rr1;
  static float b, cc, dd, e, p, q, rr, s, t, u, v, w, x;
  static float q0, s2;
  static float si;
  

/*     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR */
/*             A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION */
/*     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION */

/*     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K)) */
/*     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K) */
/*     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K) */


/*     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A" */
/*     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380 */


    if (pAlpha == aa) {
	goto L1;
    }
    if (pAlpha < (float)1.) {
	goto L12;
    }

/*     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED */

    aa = pAlpha;
    s2 = pAlpha - (float).5;
    s = sqrt(s2);
    dd = sqrt32 - s * (float)12.;

/*     STEP  2:  T=STANDARD NORMAL DEVIATE, */
/*               X=(S,1/2)-NORMAL DEVIATE. */
/*               IMMEDIATE ACCEPTANCE (I) */

L1:
    t = getNormal();
    x = s + t * (float).5;
    ret_val = x * x;
    if (t >= (float)0.) {
	return ret_val;
    }

/*     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S) */

    u = pGenerator -> asDouble();
    if (dd * u <= t * t * t) {
	return ret_val;
    }

/*     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY */

    if (pAlpha == aaa) {
	goto L4;
    }
    aaa = pAlpha;
    rr = (float)1. / pAlpha;
    q0 = ((((((q7 * rr + q6) * rr + q5) * rr + q4) * rr + q3) * rr + q2) 
	    * rr + q1) * rr;

/*               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A */
/*               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND */
/*               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS */

    if (pAlpha <= (float)3.686) {
	goto L3;
    }
    if (pAlpha <= (float)13.022) {
	goto L2;
    }

/*               CASE 3:  A .GT. 13.022 */

    b = (float)1.77;
    si = (float).75;
    cc = (float).1515 / s;
    goto L4;

/*               CASE 2:  3.686 .LT. A .LE. 13.022 */

L2:
    b = s2 * (float).0076 + (float)1.654;
    si = (float)1.68 / s + (float).275;
    cc = (float).062 / s + (float).024;
    goto L4;

/*               CASE 1:  A .LE. 3.686 */

L3:
    b = s + (float).463 - s2 * (float).178;
    si = (float)1.235;
    cc = (float).195 / s - (float).079 + s * (float).016;

/*     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE */

L4:
    if (x <= (float)0.) {
	goto L7;
    }

/*     STEP  6:  CALCULATION OF V AND QUOTIENT Q */

    v = t / (s + s);
    if (fabs(v) <= (float).25) {
	goto L5;
    }
    q = q0 - s * t + t * (float).25 * t + (s2 + s2) * log(v + (float)1.);
    goto L6;
L5:
    q = q0 + t * (float).5 * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + 
	    a3) * v + a2) * v + a1) * v;

/*     STEP  7:  QUOTIENT ACCEPTANCE (Q) */

L6:
    if (log((float)1. - u) <= q) {
	return ret_val;
    }

/*     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE */
/*               U= 0,1 -UNIFORM DEVIATE */
/*               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE */

L7:
    e = -log(1 - pGenerator -> asDouble());
    u = pGenerator -> asDouble();
    u = u + u - (float)1.;
    rr1 = si * e;
    t = b + copysign(rr1, u);

/*     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719 */

    if (t < (float)-.7187449) {
	goto L7;
    }

/*     STEP 10:  CALCULATION OF V AND QUOTIENT Q */

    v = t / (s + s);
    if (fabs(v) <= (float).25) {
	goto L8;
    }
    q = q0 - s * t + t * (float).25 * t + (s2 + s2) * log(v + (float)1.);
    goto L9;
L8:
    q = q0 + t * (float).5 * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + 
	    a3) * v + a2) * v + a1) * v;

/*     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8) */

L9:
    if (q <= (float)0.) {
	goto L7;
    }
    if (q <= (float).5) {
	goto L10;
    }
    w = exp(q) - (float)1.;
    goto L11;
L10:
    w = ((((e5 * q + e4) * q + e3) * q + e2) * q + e1) * q;

/*               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8 */

L11:
    if (cc * fabs(u) > w * exp(e - t * (float).5 * t)) {
	goto L7;
    }
    x = s + t * (float).5;
    ret_val = x * x;
    return ret_val;

/*     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.)) */

L12:
    aa = (float)0.;
    b = pAlpha * (float).3678794 + (float)1.;
L13:
    p = b * pGenerator -> asDouble();
    if (p >= (float)1.) {
	goto L14;
    }
    ret_val = exp(log(p) / pAlpha);
    if (-log(1 - pGenerator -> asDouble()) < ret_val) {
	goto L13;
    }
    return ret_val;
L14:
    ret_val = -log((b - p) / pAlpha);
    if (-log(1 - pGenerator -> asDouble()) < ((float)1. - pAlpha) * log(ret_val)) {
	goto L13;
    }
    return ret_val;
}

