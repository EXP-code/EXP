***********************************************************************
*                                                                     *
*                     Routines for SLEDGE                             *
*     (Sturm-Liouville Estimates Determined by Global Errors)         *
*                                                                     *
***********************************************************************
C
C     Release 2.2     04/12/93
C
C     Steven Pruess,  Colorado School of Mines
C                     spruess@mines.colorado.edu
C     Charles Fulton, Florida Institute of Technology
C                     fulton@zach.fit.edu
C
***********************************************************************
*      These routines estimate eigenvalues, eigenfunctions and/or     *
*      spectral density functions for Sturm-Liouville problems.       *
*      The differential equation has the form:                        *
*                                                                     *
*           -(p(x)u')' + q(x)u  =  EV*r(x)u       for x in [A,B]      *
*                                                                     *
*      with boundary conditions (at regular points)                   *
*                                                                     *
*           A1*u - A2*(pu')  =  EV*(A1'*u - A2'*(pu'))    at A        *
*           B1*u + B2*(pu')  =  0                         at B .      *
*                                                                     *
*      The functions p(x) and r(x) are assumed to be positive in      *
*      the open interval (A,B).                                       *
***********************************************************************
C
C      Possible outputs are:
C         a set of eigenvalues;
C         a set of eigenvalues and tables of values for their eigen-
C            functions;
C         a table of the spectral density function (for cases with
C            continuous spectrum).
C         a classification of the problem (regular or singular; if
C            singular then limit point or limit circle, oscillatory
C            or nonoscillatory).
C
C      The code can find eigenvalues and eigenfunctions for problems 
C      in spectral category 1 (both endpoints NONOSC), spectral
C      category 2 (one endpoint NONOSC and the other O-NO), and
C      those discrete eigenvalues below the essential spectrum in
C      spectral category 10 (both endpoints O-NO).  Here OSC at an
C      endpoint means the Sturm-Liouville equation is oscillatory for
C      all real values of EV at that endpoint, NONOSC at an endpoint
C      means the equation is nonoscillatory for all real values of EV
C      at that endpoint, and O-NO means there is a `cutoff' value EV'
C      such that the equation is nonoscillatory for real values of 
C      EV < EV' and oscillatory for real values of EV > EV'.  For 
C      problems in other spectral categories an error return will 
C      be generated.  The manner in which SLEDGE classifies singular
C      endpoints of Sturm-Liouville problems as LP/LC (Limit Point/ 
C      Limit Circle), OSC/NONOSC/O-NO, and uses this information to 
C      determine the spectral category is explained in detail in 
C      reference [2].
C
C      There is one subroutine called SLEDGE of direct interest to the
C      user; additionally, a secondary routine INTERV is available which
C      determines the indices of eigenvalues located in a specified 
C      subinterval of the real line.
C
C      The names of other routines in this package are AITKEN, ASYMEV,
C      ASYMR, BRCKET, CLASS, CLSEND, DENSEF, DSCRIP, EXTRAP, GETEF,
C      GETRN, MESH, POWER, PQRINT, REGULR, SHOOT, START, STEP, and
C      ZZERO.
C
C      There are 4 blocks of labeled COMMON with the names SLREAL,
C      SLINT, SLLOG, and SLCLSS.
C
C      This is the double precision version of the code; all floating
C      point variables should be declared DOUBLE PRECISION in the
C      calling program.  In these subprograms all such local
C      variables and constants have been explicitly declared; also,
C      FORTRAN77 generic intrinsic functions have been used, so
C      conversion to single precision should be straightforward, if
C      desired.
C
C      ACKNOWLEDGMENT:  This work was partially supported by the
C      National Science Foundation under grants DMS-8813113 and DMS-
C      8905202 to Florida Institute of Technology and DMS-8800839 and
C      DMS-8905232 to the Colorado School of Mines.
C
C      References
C
C      The following papers are available from the authors on request:
C
C      [1]. Pruess & Fulton, Mathematical software for Sturm-Liouville
C           problems, ACM Trans. on Math. Software, to appear, 1993.
C
C      [2]. Fulton, Pruess & Xie, The automatic classification of Sturm-
C           Liouville problems, submitted, 1992.
C
C      [3]. Pruess, Fulton & Xie, An asymptotic numerical method for a
C           class of singular Sturm-Liouville problems, submitted, 1992.
C
C      [4]  Fulton and Pruess, Eigenvalue and eigenfunction asymptotics
C           for regular Sturm-Liouville problems, Jour. Math. Anal. and
C           Appls., to appear. 
C
C      [5]  Fulton and Pruess, Numerical Approximation of singular
C           spectral functions arising from the Fourier-Jacobi problem
C           on a half line with continuous spectra, Sixth International
C           Workshop in Analysis and its Applications, June, 1992.

C      [6]. Pruess, Fulton & Xie, Performance of the Sturm-Liouville 
C           software package SLEDGE, Colo. School of Mines, Dept. of
C           Math. and Comp. Sci., MCS-91-19, 1991. Revision 12/92.
C-----------------------------------------------------------------------
C      Brief overview of algorithms:
C
C         The code constructs (or takes from input) an initial mesh,
C      called the level 0 mesh.  Subsequent meshes (for level 1,2,...)
C      are unions of the previous level's mesh with its midpoints.  A 
C      sequence of estimates for desired eigenvalues and eigenfunctions 
C      is constructed, one set for each level.  These estimates (the 
C      eigenvalue is called EvHat) are exact solutions (up to the 
C      requested tolerance) of a Sturm-Liouville problem which is an 
C      approximation to the original one; this approximation results 
C      from replacing the given coefficient functions with step function
C      approximations relative to the current level's mesh.  The eigen-
C      functions of the resulting ODE's are piecewise trigonometric
C      (circular or hyperbolic) functions.
C         If estimates for the spectral density function are reqested,
C      these are computed as limits of a sequence of spectral density
C      functions of approximating regular problems.  For these regular
C      problems the spectral density function is a step function, and
C      is computed directly from the definition making use of computed
C      eigenvalues and the norm reciprocals of the corresponding eigen-
C      functions.  If verbose output is rquested by the user, there
C      will be displayed iterations (corresponding to the sequence of
C      approximating regular intervals which the code automatically 
C      selects) and within each iteration there will be levels 
C      (corresponding to increasingly finer meshes as described above).
C      A step spectral density function will be printed at each level
C      of each iteration.  The spectral density function displayed at
C      the end of each iteration is the result of an h-squared extra-
C      polation over the regular step functions generated at each level 
C      of this iteration.  The condition for stopping at a given iter-
C      ation is a straightforward comparison of the spectral function 
C      data for the current iteration with the previous iteration. There
C      is no extrapolation over the sequence of regular approximating
C      intervals as no extrapolation theory for the approximation of
C      the singular spectral function by regular step spectral functions
C      is known.  (To achieve closer approximation of the regular step
C      spectral functions to the singular spectral function, it is
C      actually the piecewise linear function obtained by joining the
C      midpoints of successive steps by a straight line which is used as
C      the `regular' spectral function for the purpose of generating the
C      actual data used for the h-squared extrapolation.)
C         The classification is determined by applying standard theory
C      to an approximating problem, each of whose coefficient functions,
C      in a small neighborhood of each endpoint, consists of the leading
C      term in a power-like asymptotic development.  For this reason
C      there are many problems, particularly those with oscillatory
C      coefficient functions, for which the code's output for the
C      classification information is labelled `uncertain'.  For further
C      information on the theory used by the code to generate endpoint
C      classifications and spectral category information see [2] above.
C-----------------------------------------------------------------------
C  Usage (simple explanation) -
C      The subroutine SLEDGE is called in the following manner:
C
C      SUBROUTINE SLEDGE(JOB,CONS,ENDFIN,INVEC,TOL,TYPE,EV,NUMX,XEF,EF,
C                        PDEF,T,RHO,IFLAG,STORE)
C
C      If k eigenvalues (no eigenvectors) are sought then set
C      (a) the logical 5-vector JOB to (True, False, False, False, True);
C      (b) the real 8-vector CONS to the values of A1,A1',A2,A2',B1,B2,
C          A,B for the boundary condition information.  It does not
C          matter what values are used for infinite endpoints, nor for
C          the boundary constants at a singular endpoint; the code 
C          automatically selects the Friedrichs' boundary condition at 
C          NONOSC singular endpoints, overriding user input for the 
C          boundary condition constants; for infinite endpoints the code
C          also automatically selects these constants.
C      (c) the logical 2-vector ENDFIN to (True, True) if both endpoints
C          are finite, (True, False) if A is finite but B infinite, etc.;
C      (d) the integer vector INVEC should have
C              INVEC(1) = 0 (no internal printing)
C              INVEC(2) = 0
C              INVEC(3) = k, the number of eigenvalues sought
C              INVEC(3+i) = index of ith eigenvalue sought,i = 1,...,k;
C      (e) the real 6-vector TOL should have
C              TOL(1) = absolute error tolerance desired,
C              TOL(2) = relative error tolerance desired,
C          the remaining 4 entries of TOL are ignored;
C      (f) the output estimate for the ith eigenvalue is returned
C          in EV(i), i = 1,...,k;
C      (g) the output integer k-vector IFLAG(*) should have all entries
C          zero; nonzero values indicate warnings or error returns and 
C          are explained in the detailed usage section below;
C      (h) the auxiliary vector STORE(*) should be dimensioned at least
C          155 in the calling program;
C      (i) the logical 4 by 2 vector TYPE, the real vectors XEF(1), 
C          EF(1), PDEF(1), T(1), and RHO(1) can be ignored except that
C          they need to be declared in the calling program.  The integer
C          scalar NUMX can also be ignored.
C
C      If k eigenfunctions are also desired, then follow the above 
C      pattern except make JOB(1) False and JOB(2) True.  The values of
C      TOL(3) and TOL(4) control the absolute and relative errors in 
C      each u(x); TOL(5) and TOL(6) control the absolute and relative 
C      errors in each (pu')(x).  It is usually appropriate to set TOL(5)
C      = TOL(3) = TOL(1) and TOL(6) = TOL(4) = TOL(2), but the user has 
C      the option of entering all six tolerance parameters as desired.
C      The output eigenfunction information is returned in the three 
C      real vectors X(*) for the independent variable x , EF(*) for u(x),
C      PDEF(*) for (pu')(x).  The code automatically chooses the x 
C      values; the number of values is returned in NUMX.  If you prefer
C      another choice of output points, see the detailed explanation 
C      below on usage of the code.  The values for the first requested 
C      u(x) are returned in the first NUMX locations of EF(*), those for
C      the second are in the next NUMX locations, etc.  PDEF(*) is part-
C      itioned similarly; X(*) must be dimensioned at least 31 in the 
C      calling program while EF(*) and PDEF(*) must be dimensioned at 
C      least 31*k.  The auxiliary vector STORE(*) should be dimensioned 
C      at least 420.
C
C      For other possibilities, see the detailed description which 
C      follows.
C-----------------------------------------------------------------------
C  Usage (detailed explanation) -
C
C  SUBROUTINE SLEDGE(JOB,CONS,ENDFIN,INVEC,TOL,TYPE,EV,NUMX,XEF,EF,PDEF,
C                    T,RHO,IFLAG,STORE)
C
C  Input parameters;
C     JOB(*)    = logical 5-vector,
C                 JOB(1) = .True. iff a set of eigenvalues are to be
C                                 computed but not their eigenfunctions.
C                 JOB(2) = .True. iff a set of eigenvalue and eigenfunc-
C                                 tion pairs are to be calculated.
C                 JOB(3) = .True. iff the spectral function is to be
C                                 computed over some subinterval of the
C                                 essential spectrum.
C                 JOB(4) = .True. iff the normal call to the routines for 
C                                 classification (regular/singular, etc.)
C                                 is OVERRIDDEN.  If JOB(4) is True then
C                                 TYPE(*,*) discussed below must be 
C                                 INPUT correctly!  Most users will not 
C                                 want to override the classification 
C                                 routines, but it would, of course, be 
C                                 appropriate for users experimenting with 
C                                 problems for which the coefficient 
C                                 functions do not have power-like
C                                 behavior near the singular endpoints.
C                                 Note: the code may perform poorly if 
C                                 the classification information is 
C                                 incorrect; since the cost is usually 
C                                 negligible, it is strongly recommended 
C                                 that JOB(4) be False.  The classifica-
C                                 tion is deemed sufficiently important 
C                                 for spectral density function calcul-
C                                 ations that JOB(4) is ignored when 
C                                 the input JOB(3) is True.
C                 JOB(5) = .True. iff mesh distribution is to be chosen
C                                 by SLEDGE.  If JOB(5) is True and NUMX
C                                 is zero, then the number of mesh
C                                 points is also chosen by SLEDGE; if
C                                 NUMX > 0 then NUMX mesh points will be
C                                 used.  If JOB(5) is False, then the
C                                 number (NUMX) and distribution
C                                 (XEF(*)) must be input by the user.
C                                 If JOB(3) is True and JOB(5) False
C                                 then the user must set BOTH the number
C                                 NUMX and distribution.  In this case,
C                                 NO global error estimates are made.
C     CONS(*)    = real vector of length 8, values are the boundary
C                  condition constants A1, A1', A2, A2', B1, B2, A, B.
C                  In the case of a NONOSC singular endpoint, the class-
C                  ification routine uses the Friedrichs' boundary 
C                  condition constants.  The code cannot automatically 
C                  choose a non-Friedrichs' boundary condition; however,
C                  interval truncation in the user's calling program can
C                  be used, together with many calls to SLEDGE, to 
C                  compute singular eigenvalues associated with a non-
C                  Friedrichs' boundary condition at a NONOSC endpoint
C                  (see remark 12 below).
C     ENDFIN(*) = logical 2-vector, values are
C        ENDFIN(1) = .True. iff endpoint A is finite.
C        ENDFIN(2) = .True. iff endpoint B is finite.
C     INVEC(*)  = integer vector of length 3+(number of eigenvalues
C                 desired).  This vector contains a variety of input
C                 information.
C        INVEC(1) controls the amount of internal printing: values are
C                 from 0 (no printing) to 5 (much printing).
C                 For INVEC(1) > 0 much of the output will be to a file
C                 attached to unit #21 which should be named in the
C                 user's calling program via an OPEN statement.
C                 Output for the various cases is, when INVEC(1) =
C                    0   no printing.
C                 When JOB(1) or JOB(2) is True
C                    1   initial mesh (the first 51 or fewer points),
C                        eigenvalue estimate at each level,
C                    4   the above,
C                        at each level
C                          matching point for eigenfunction shooting,
C                          X(*), EF(*), PDE(*) values,
C                    5   all the above,
C                        at each level
C                          brackets for the eigenvalue search,
C                          intermediate shooting info for the eigen-
C                          function and eigenfunction norm.
C                 When JOB(3) is True
C                    1   the actual (a,b) used at each iteration,
C                        the total number of eigenvalues computed,
C                    2   the above,
C                        switchover points to the asymptotic formulas,
C                        some intermediate Rho(t) approximations,
C                    3   all the above,
C                        initial meshes for each iteration,
C                        index of the largest EV which may be computed,
C                        various Ev and RsubN values,
C                    4   all of the above,
C                        RhoHat values at each level,
C                    5   all of the above,
C                        all Ev and RsubN values below switchover point.
C                 When JOB(4) is False
C                    2   output a description of the spectrum,
C                    3   the above plus the constants for the
C                        Friedrichs' boundary condition(s),
C                    5   all the above plus intermediate details of
C                        classification calculation.
C                 Some of the output may go to the default output device
C                 (screen or printer), but all information requested is
C                 also directed to the file attached to unit #21.
C        INVEC(2) gives the number (positive) of output values desired
C                 for the array RHO(*) (not referenced if JOB(3) is
C                 False).
C        INVEC(3) is the total number of eigenvalues to be output in
C                 EV(*).
C        INVEC(J) for J = 4, 5, ..., 3+INVEC(3) contains the indices for
C                 the eigenvalues sought. If JOB(1) and JOB(2) are
C                 False, this part of INVEC(*) is not referenced.
C     TOL(*)    = real vector of from 2 to 6 tolerances. 
C                 If JOB(1) or JOB(2) is True then
C                   TOL(1) is the absolute error tolerance for e-values,
C                   TOL(2) is the relative error tolerance for e-values,
C                   TOL(3) is the abs. error tolerance for e-functions,
C                   TOL(4) is the rel. error tolerance for e-functions,
C                   TOL(5) is the abs. error tolerance for eigenfunction
C                          derivatives,
C                   TOL(6) is the rel. error tolerance for eigenfunction
C                          derivatives.
C                   Eigenfunction tolerances need not be set if JOB(2)
C                   is False.
C                 If JOB(3) is True then
C                   TOL(1) is the absolute error tolerance,
C                   TOL(2) is the relative error tolerance;
C                   the output RHO values are NOT required to satisfy
C                   these tolerances when JOB(5) is False.
C                 All absolute error tolerances must be positive; all
C                 relative error tolerances must be at least 100 times 
C                 the unit roundoff.
C     NUMX      = integer whose value is
C                    the number of output points where each eigen-
C                    function is to be evaluated (the number of entries
C                    in XEF(*)) when JOB(2) is True,
C                 or
C                    the number of points in the initial mesh used when
C                    JOB(5) is False and NUMX>0.
C                 If JOB(5) is False, the points in XEF(*) should be
C                 chosen to have a reasonable distribution.  Since the
C                 endpoints A and B must be part of any mesh, NUMX
C                 cannot be 1 in this case.  If JOB(5) is FALSE and
C                 JOB(3) is True, then NUMX must be positive.
C     XEF(*)    = real vector of points where 
C                    eigenfunction estimates are desired (JOB(2) True)
C                 or
C                    where user's initial mesh is entered (JOB(5) False
C                    and NUMX>0).
C                 The values must satisfy
C                      A = XEF(1) < XEF(2) < ... < XEF(NUMX) = B .
C                 When JOB(2) is True the initial mesh corresponds to
C                 the set of points where eigenfunction output is
C                 desired. If JOB(2) is False and NUMX = 0, then this
C                 vector is not referenced.  When A and/or B are
C                 infinite (as indicated through ENDFIN(*)), the 
C                 entries XEF(1) and/or XEF(NUMX) are ignored; however,
C                 it is required that XEF(2) be negative when ENDFIN(1)
C                 is False, and XEF(NUMX-1) be positive when ENDFIN(2)
C                 is False (otherwise, IFLAG = -39 will result).
C     T(*)      = real vector of INVEC(2) values where the spectral
C                 function RHO(*) is desired (the existence and location
C                 of continuous spectrum can be found by first calling
C                 SLEDGE with JOB(J) False, J=1,...,4 and INVEC(1) = 1).
C                 Vector T(*) is not referenced if JOB(3) is False.  Its
C                 entries must be in increasing order.
C
C   Output parameters:
C     TYPE(*,*) = 4 by 2 logical array; column 1 carries information
C                 about endpoint A while column 2 refers to B.
C                 TYPE(1,*) = True  iff the endpoint is regular,
C                 TYPE(2,*) = True  iff it is limit circle,
C                 TYPE(3,*) = True  iff it is nonoscillatory for all EV,
C                 TYPE(4,*) = True  iff it is oscillatory for all EV,
C                 Important note: all of these must be correctly INPUT
C                 if JOB(4) is True!
C     EV(*)     = real vector containing the computed approximations to
C                 the eigenvalues whose indices are specified in
C                 INVEC(*); if JOB(1) and JOB(2) are False, then the
C                 output has no meaning.
C     NUMX      = the number of output points for eigenfunctions when
C                 input NUMX = 0, and JOB(2) or JOB(5) is True.
C     XEF(*)    = input values (if any) are changed only if JOB(2) and
C                 JOB(5) are True; in this case, the output values
C                 are chosen by the code.  If JOB(2) is False then this
C                 vector is not referenced; if JOB(2) is True and NUMX>0
C                 on input then XEF(*) should be dimensioned at least 
C                 NUMX+16 in the calling program.  If JOB(2) is True and
C                 NUMX=0 on input (so that the code chooses NUMX), then
C                 dimension XEF(*) at least 31 in the calling program.
C     EF(*)     = real vector of eigenfunction values: EF((k-1)*NUMX+i)
C                 is the estimate of u(XEF(i)) corresponding to the
C                 eigenvalue in EV(k).  If JOB(2) is False then this
C                 vector is not referenced.  Otherwise, if JOB(2) is
C                 True and NUMX>0 on input then EF(*) should be
C                 dimensioned at least NUMX*INVEC(3) in the calling
C                 program.  If JOB(2) is True and NUMX=0 on input (so
C                 that the code chooses NUMX), then dimension XEF(*)
C                 at least 31*INVEC(3) in the calling program.
C     PDEF(*)   = real vector of eigenfunction derivative values: 
C                 PDEF((k-1)*NUMX+i) is the estimate of (pu')(XEF(i))
C                 corresponding to the eigenvalue in EV(k).  If JOB(2)
C                 is False then this vector is not referenced; otherwise,
C                 it must be dimensioned as is EF(*).
C     RHO(*)    = real vector of values for the spectral density
C                 function rho(t), RHO(I) = rho(T(I)).  RHO(*) must be
C                 dimensioned at least INVEC(2); this vector is not
C                 referenced if JOB(3) is False.
C     IFLAG(*)   = integer vector carrying information about the output.
C       Declared length must be at least max(1,INVEC(3)).  For the Kth
C       requested eigenvalue (when JOB(1) or JOB(2) is true; otherwise,
C       only IFLAG(1) is used):
C       IFLAG(K) =  0, normal return, output should be reliable.
C                <  0, fatal error, calculations ceased: if
C                = -1, too many levels needed for the eigenvalue
C                      calculation; problem seems too difficult for 
C                      this algorithm at this tolerance. Are the 
C                      coefficient functions nonsmooth?  
C                = -2, too many levels needed for the eigenfunction
C                      calculation; problem seems too difficult for 
C                      this algorithm at this tolerance.  Are the 
C                      eigenfunctions ill-conditioned?
C                = -3, too many levels needed for the spectral density
C                      calculation; problem seems too difficult for 
C                      this algorithm at this tolerance.
C                = -4, the user has requested the spectral density
C                      function for a problem which has no continuous
C                      spectrum.
C                = -5, the user has requested the spectral density 
C                      function for a problem with both endpoints
C                      generating essential spectrum, i.e., both
C                      endpoints being either OSC or O-NO.  The spectral
C                      density function calculation has not been 
C                      implemented for such cases.  For spectral 
C                      category 10 (both endpoints O-NO) the spectral
C                      multiplicity is generally two, proper normal-
C                      izations for the solutions against which the 
C                      spectral functions will be normalized will 
C                      depend on how the user wants to express the 
C                      eigenfunction expansion.  Users having problems 
C                      in spectral category 10 are encouraged to supply 
C                      them to the authors, and if possible, recommend 
C                      normalizations of the two solutions to be used in
C                      writing the associated eigenfunction expansion.
C                = -6, the user has requested the spectral density 
C                      function for a problem in spectral category 2 for
C                      which a proper normalization of solution at the
C                      NONOSC endpoint is not known; for example, 
C                      problems with an irregular singular point or 
C                      infinite endpoint at one end and continuous 
C                      spectrum generated at the other.  Users with 
C                      problems of this type are encouraged to supply 
C                      them to the authors, and if possible, recommend a
C                      normalization of solution at the NONOSC endpoint
C                      which they would like to see implemented. As a
C                      rule it is best to pick a normalization which
C                      ensures that the solution is uniquely fixed and
C                      entire in the eigenvalue parameter EV for all x
C                      in the Sturm-Liouville interval; for further
C                      mathematical information on NONOSC endpoints we
C                      refer to paper [2] above.
C                = -7, problems encountered in obtaining a bracket.
C                = -8, too small a step used in the integration;
C                      TOL(*) values may be too small for this problem.
C                = -9, too small a step used in a spectral density
C                      function calculation for which the continuous
C                      spectrum is generated by a finite endpoint.  Try
C                      transforming to Liouville (or some other) form.
C                = -10, an argument to the circular trig functions is
C                       too large.  Try rerunning with a finer initial
C                       mesh, or, on singular problems, use interval
C                       truncation (see remark (12)).
C                = -15, p(x) and r(x) not positive in (A,B).
C                = -20, eigenvalues/functions were requested for a
C                       problem with an OSC singular endpoint.
C                       Interval truncation (see remark (12)) must be
C                       used on such problems.
C                = -3?, illegal input, viz.
C                  -30,  NUMX = 1 when JOB(5) is True,
C                        or NUMX = 0 when JOB(3) is True and JOB(5) is 
C                        False,
C                  -31,  B1 = B2 = 0 (at a regular endpoint),
C                  -32,  A1'*A2-A1*A2' .le. 0 when A1' or A2' nonzero,
C                  -33,  A1 = A2 = A1'= A2'= 0 (at a regular endpoint),
C                  -34,  A .ge. B (when both are finite),
C                  -35,  TOL(odd) .le. 0 ,
C                  -36,  TOL(even)  <  100*unit roundoff,
C                  -37,  INVEC(k) < 0  for some k>3 when INVEC(3)>0,
C                  -38,  INVEC(2) .le. 0 when JOB(3) is True ,
C                  -39,  XEF(*) entries out of order or not in [A,B].
C                        or XEF(2), XEF(NUMX-1) have the wrong sign in
C                           infinite interval cases,
C                        or T(*) entries are out of order.
C                >  0,  indicates some kind of warning, in this case the
C                       value may contain ANY of the following digits:
C                =  1,  failure in routine BRCKET probably due to a
C                       cluster of eigenvalues which the code cannot
C                       separate.  Calculations have continued as best
C                       as possible, but any eigenfunction results are
C                       suspect.  Try rerunning with tighter input
C                       tolerances to separate the cluster. 
C                =  2,  there is uncertainty in the classification for
C                       this problem.  Because of the limitations of the
C                       floating point arithmetic on the computer used,
C                       and the nature of the finite sampling, the
C                       routine is cannot be decisive about the 
C                       classification information at the requested
C                       tolerance.
C                =  3,  there may be some eigenvalues imbedded in the
C                       essential spectrum; using IPRINT greater than
C                       zero will result in additional output giving 
C                       the location of the approximating eigenvalues 
C                       for the step function problem.  These could be
C                       extrapolated to estimate the actual eigenvalue
C                       embedded in the essential spectrum.
C                =  5,  a change of variables was made to avoid poten-
C                       tial slow convergence; however, the global
C                       error estimates may not be as reliable.  Some
C                       experimentation using different tolerances is
C                       recommended.
C                =  6,  there were problems with eigenfunction conver-
C                       gence in a spectral density calculation; the
C                       output Rho(t) may not be accurate.
C
C   Auxiliary storage:
C     STORE(*) = real vector of auxiliary storage, must be dimensioned
C                at least
C             max(155,NUMX+16)     in general;
C               26*(NUMX+16)       for any eigenfunction calculation;
C             2400+13*INVEC(2)     for any spectral density calculation.
C-----------------------------------------------------------------------
C  SUBROUTINE INTERV(FIRST,ALPHA,BETA,CONS,ENDFIN,NFIRST,NTOTAL,
C                    IFLAG,STORE)
C
C    Input parameters:
C     FIRST      = logical; value is True if various internal variables
C                  have not yet been set.  If a prior call has been made
C                  to INTERV with FIRST True, then a little time can
C                  be saved by letting FIRST be False.
C                  IMPORTANT NOTE: setting FIRST = True will clobber any
C                  initial mesh the user has input (when NUMX > 0 or
C                  JOB(5) is False);  also, INTERV will classify the
C                  problem irregardless of what JOB(4) is set to
C                  for SLEDGE.
C      ALPHA     = real value of left end point of search interval.
C      BETA      = real value of right end point of search interval.
C      CONS(* )  = real vector of 8 input constants: A1, A1', A2, A2',
C                  B1, B2, A, B.
C      ENDFIN(*) = logical 2-vector, same meaning as in SLEDGE.
C      STORE(*)  = real vector holding initial mesh.
C
C    Output parameters:
C      NFIRST = index of first eigenvalue > ALPHA.
C      NTOTAL = total number of eigenvalues in the interval.
C      IFLAG  = integer status indicator.
C               IFLAG =   0 , normal return, output should be reliable,
C                     =  11 , there are no eigenvalues in [alpha, beta],
C                     =  12 , low confidence in NFIRST or NTOTAL or both,
C                     =  13 , BETA and/or ALPHA exceed the cutoff for 
C                             the continuous spectrum.  If only BETA
C                             is too big then NFIRST may be OK, but
C                             NTOTAL is meaningless.
C                     = -11 , ALPHA .ge. BETA,
C                     = -25 , oscillatory endpoint, output meaningless,
C                     = -3? , illegal CONS(*) values (see above comments
C                             on SLEDGE for an explanation). 
C--------------------------------------------------------------------------
C         In addition, a subroutine subprogram must be provided for the
C      coefficient functions p(x), q(x), and r(x); the form of this
C      routine is
C
C          SUBROUTINE COEFF(X,PX,QX,RX)
C          DOUBLE PRECISION X,PX,QX,RX
C               ...
C          PX = ...
C          QX = ...
C          RX = ...
C          RETURN
C          END
C
C      The subroutine name MUST be COEFF, though of course the names of
C      arguments only need follow the usual FORTRAN77 rules. X is the
C      independent variable; PX, QX, and RX are the output values of the
C      respective coefficient functions p(x), q(x), and r(x) at X.
C-----------------------------------------------------------------------
C     This is a simple sample driver for SLEDGE.
CC
CC     Declare all variables:
CC
C      INTEGER IFLAG(1),INVEC(4),NUMX, I,J,K
C      LOGICAL JOB(5),TYPE(4,2),ENDFIN(2)
C      DOUBLE PRECISION CONS(8),TOL(6),EV(1),T(3),RHO(3),STORE(2450),
C     &                 XEF(5),EF(5),PDEF(5)
CC
CC     Load the boundary condition information into CONS(*).
CC     This example has a Neumann condition at A = 1, and a
CC     singular point at B = +infinity.
CC
C      DATA CONS/0.0, 0.0, 1.0, 0.0,   0.0, 0.0,   1.0, 0.0/
C      DATA ENDFIN/.TRUE., .FALSE./
CC
CC     The eigenfunctions will be estimated at 5 points.
CC
C      DATA NUMX,XEF/5, 1.0, 1.5D0, 2.0, 4.0, 100.0/
CC
CC     Initialize the vector INVEC(*):
CC        little printing,
CC        3 output points for the density function Rho(t),
CC        estimates for the first (index 0) eigenvalue/function.
CC
C      DATA INVEC/1, 3, 1, 0/
CC
CC     Set the JOB(*) vector:
CC        estimate both eigenvalues and eigenvectors,
CC        estimate the spectral density function,
CC        classify,
CC        force the initial mesh to be the output points.
CC
C      DATA JOB/.FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE./
CC
CC     Set the tolerances:
CC
C      DATA TOL/1.D-5,1.D-4,  1.D-5,1.D-4,  1.D-5,1.D-4/
CC
CC     Initialize the 3 output points for the density function.
CC
C      DATA T/0.0, 0.5, 2.0/
CC
CC     Open file for output.
CC
C      OPEN(21,FILE = 'sample.out')
C      CALL SLEDGE(JOB,CONS,ENDFIN,INVEC,TOL,TYPE,EV,NUMX,XEF,EF,PDEF,
C     &            T,RHO,IFLAG,STORE)
CC
CC     Print results:
CC
C      DO 30 I = 1,INVEC(3)
C         WRITE (*,10) INVEC(3+I),EV(I),IFLAG(I)
C         WRITE (21,10) INVEC(3+I),EV(I),IFLAG(I)
C   10    FORMAT(' Nev =',I6,';   Ev =',D25.15,';     Flag = ',I3)
C         IF (IFLAG(I) .GT. -10) THEN
C            WRITE (*,15)
C            WRITE (21,15) 
C   15       FORMAT(13X,'x',23X,'u(x)',18X,'(pu`)(x)')
C            K = NUMX*(I-1)
C            DO 25 J = 1,NUMX
C               WRITE (21,20) XEF(J),EF(J+K),PDEF(J+K)
C   20          FORMAT(3D25.15)
C   25       CONTINUE
C         ENDIF
C   30    CONTINUE
C      WRITE (*,35)
C      WRITE (21,35)
C   35 FORMAT(/,8X,'t',21X,'Rho(t)')
C      DO 45 I = 1,INVEC(2)
C         WRITE (*,40) T(I),RHO(I)
C         WRITE (21,40) T(I),RHO(I)
C   40    FORMAT(F11.3,D32.15)
C   45 CONTINUE
C      CLOSE(21)
C      STOP
C      END
CC
C      SUBROUTINE COEFF(X,PX,QX,RX)
CC
CC     Define the coefficient functions; here a Yukawa potential.
CC
C      DOUBLE PRECISION X,PX,QX,RX, T
CC
CC     Be careful with potential over/underflows; here we assume the
CC     IEEE double precision exponent range.
CC
C      IF (X .LT. 650.0) THEN
C         T = EXP(-X)
C      ELSE
C         T = 0.0
C      ENDIF
C      PX = 1.0
C      QX = -T/X   
C      RX = PX
C      RETURN
C      END
CC
CC     End of sample driver for SLEDGE.
C-----------------------------------------------------------------------
C      General remarks:
C      (1) Two machine dependent constants must be set in a DATA
C          statement in routine START (in part 4 of the package):
C          URN   -  an estimate of the unit roundoff; infinite output
C                   values are assigned the value 1/URN.
C          UFLOW -  a number somewhat smaller than -ln(underflow level).  
C                   Values of certain variables z for which
C                             ln(abs(z)) < -under
C                   will be set to zero.
C      (2) A value of IFLAG = -1, -2, or -3 may be the result of a 
C          lack of smoothness in the coefficient functions.  In such
C          cases a user input mesh may perform better (see (4) below).
C      (3) The heuristics for generating the initial mesh distribution
C          work reasonably well over a wide range of examples, but
C          occasionally they are far from optimal.  The code's choice
C          can be over-ridden by setting JOB(5) False, setting NUMX
C          appropriately and supplying a mesh in XEF(*).
C      (4) If any of the coefficient functions p,q, or r (or their first
C          few derivatives) have finite jump discontinuities at points
C          in the interior of (A,B), then it is advantageous to have
C          these points in SLEDGE's mesh.  Currently, this can only be
C          accomplished by setting JOB(5) False and supplying an 
C          appropriate mesh using NUMX and XEF(*).
C      (5) In general, eigenvalue convergence is observed to be more
C          rapid than eigenfunction convergence; hence, it is
C          recommended that JOB(2) be False unless eigenfunction
C          information really is necessary.
C      (6) When eigenfunction output is sought, unless some knowledge
C          of the eigenfunction is known in advance, it is recommended 
C          that JOB(5) be True so that the code will attempt to choose
C          a reasonable distribution for the initial mesh points.
C      (7) Computing the spectral density function for problems having
C          continuous spectrum can be very expensive; it is recommended
C          that initially, relatively crude tolerances (0.001 or so) be
C          used to get some idea of the effort required.
C      (8) It is recommended that every problem be classified (JOB(4)
C          False) by the code before any calculation of spectral
C          quantities occurs.  Only if the user is certain as to what
C          the classification is (and describes it correctly through
C          INVEC and TYPE) should the classification option be bypassed.
C      (9) If the code does the classification of singular problems, it
C          will automatically choose the Friedrichs' boundary condition 
C          at NONOSC endpoints.  If another boundary condition is 
C          desired, the user must use interval truncation in the
C          calling program (see remark (12)).
C     (10) While all parts of the code should function on machines
C          with a fairly narrow exponent range (such as IEEE single
C          precision), it is better to have a relatively wide exponent
C          range (IEEE double precision).  The classification algorithm,
C          in particular, is far more reliable if done on a machine with
C          a fairly wide exponent range.
C     (11) Care must be taken in writing the subroutine COEFF for the
C          evaluation of p(x), q(x), and r(x) to avoid arithmetic
C          exceptions such as overflow and underflow (or trig function
C          arguments too large).  This can be especially delicate on
C          machines with a small exponent range.
C     (12) In some cases `interval truncation' is recommended.  By this
C          is meant the user should call SLEDGE several times using a 
C          sequence of regular endpoints (with appropriate boundary 
C          conditions) converging to the singular endpoint.  The eigen-
C          values of the regular problems selected by the user should be
C          arranged so as to converge to those of the desired singular
C          problem. For example, if the user wishes to compute eigen-
C          values associated with a non-Friedrichs' boundary condition 
C          for problems in spectral category 1, the user can experiment
C          with choosing a sequence of regular approximating intervals,
C          and vary the boundary conditions appropriately by means of a
C          `boundary condition function' or known solution of the 
C          equation for a real value of EV on the sequence of regular 
C          intervals until convergence of the regular eigenvalues to 
C          the desired singular one is observed.  Similarly, for 
C          problems in spectral category 3 or 5 which involve one or two
C          endpoints which are LC and OSC, the (necessarily discrete) 
C          spectrum is known to be unbounded below and above. To 
C          implement a given LC boundary condition at a singular LC 
C          endpoint one may choose a `boundary condition function' or 
C          known solution of the equation for a real value of EV and 
C          make use of it on a sequence of regular approximating 
C          intervals to vary the boundary condition on successive calls
C          to SLEDGE for the sequence of regular intervals until 
C          convergence to the desired singular eigenvalue is observed. 
C          At present these methods are highly experimental and problem-
C          dependent as good heuristics for the choice of the rate of 
C          convergence of the regular intervals to the singular one 
C          which work well over a wide class of problems are not known.
C          (The only case in which SLEDGE automatically selects regular
C          approximating subintervals is for spectral density function
C          calculations for problems in spectral category 2; but
C          here the singular endpoint is of LP type, so no singular
C          boundary condition is required to be implemented.)
C     (13) Problems of slow convergence can sometimes be avoided by a
C          judicious change of either dependent or independent variable
C          (or both).
C     (14) If the Liouville normal form potential Q(t) has a minimum
C          far from zero, then the heuristics for generating the initial
C          mesh may well miss it.  In this case, it is advisable to
C          shift the independent variable.
C     (15) The determination of the total number of eigenvalues is the
C          most difficult part of the classification process.  When the
C          theory provides this number, of course, there is no problem;
C          otherwise, it should be viewed with some skepticism.  A more
C          reliable count of the eigenvalues below the cutoff point of
C          the essential spectrum can be gained (at some expense) by
C          trying to compute many eigenvalues near that point.
C///////////////////////////////////////////////////////////////////////
      SUBROUTINE SLEDGE(JOB,CONS,ENDFIN,INVEC,TOL,TYPE,EV,NUMX,XEF,EF,
     &                  PDEF,T,RHO,IFLAG,STORE)
      INTEGER INVEC(*),NUMX,IFLAG(*)
      LOGICAL JOB(*),ENDFIN(*),TYPE(4,*)
      DOUBLE PRECISION CONS(*),TOL(*),EV(*),XEF(*),EF(*),PDEF(*),T(*),
     &                 RHO(*),STORE(*)
C
C     This is the interface routine between the user and other routines
C     which carry out most of the actual calculations.
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &        KCLASS(2),I,IBASE,IEV,IPRINT,J,JTOL,K,KCL1,KCL2,LASTEV,
     &        MAXITS,MAXT,MU1,MU2,NEV,NEXTRP,NUMEV,NUMT
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &        CP(2),CR(2),CUTOFF,D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     &        ETA(2,2),PNU(2),
     &        AA,ALPHA,BB,CEV(2),DENS,DENSHI,DENSLO,DENSOP,ENDFAC,
     &        ENDI(5),ERROR,FZ,HMIN,RHOTOL,SGN,TOL1,TOLMAX,XTOL,ZETA,
     &        ZETAI(5),
     &        ZERO,HALF,ONE,TWO,FOUR
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2),
     &        AAFIN,BBFIN,CSPEC(2),DOMESH,DONE,EDONE,JOBST(3),
     &        LBASE,LMESH,LPLC(2),OSCILL
      COMMON /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, HALF = 0.5D0, ONE = 1.0, TWO = 2.0, 
     &           FOUR = 4.0, TOLMAX = 1.D-4)
      DATA DENSLO,DENSOP,DENSHI/4.0, 6.0, 12.0/
      DATA ENDI/12.0, 20.0, 85.0, 240.0, 500.0/
      DATA ZETAI/2.2, 2.0, 1.5, 1.4, 1.3/
C
C     Initialize.
C
      AFIN = ENDFIN(1)
      BFIN = ENDFIN(2)
      IPRINT = INVEC(1)
      NUMT = INVEC(2)
      NEV = INVEC(3)
      LNF = .FALSE.
      DOMESH = .TRUE.
      LMESH = .FALSE.
      FLAG = 0
      IFLAG(1) = 0
      IBASE = 1
      LBASE = .FALSE.
      DO 5 I = 1,6
         LFLAG(I) = .FALSE.
    5 CONTINUE
      TOL1 = MIN(TOL(1)+TOL(2),TOLMAX)
      JOBST(1) = JOB(2)
      IF ((NUMX .GT. 0) .AND. (.NOT. JOB(5))) THEN
         JOBST(2) = .TRUE.
      ELSE
         JOBST(2) = .FALSE.
      ENDIF
      JOBST(3) = JOB(3)
      IF ((.NOT. JOB(1)) .AND. (.NOT. JOB(2))) NEV = 0
      CALL START(JOBST,CONS,TOL,NEV,INVEC(4),NUMX,XEF,NUMT,T,
     &           NEXTRP,STORE)
      IF (JOB(4)) THEN
         ALPHA = A2*A1P-A1*A2P
         IF ((A1P .NE. ZERO) .OR. (A2P .NE. ZERO)) THEN
            IF (ALPHA .LE. ZERO) FLAG = -32
         ELSE
            IF ((A1 .EQ. ZERO) .AND. (A2 .EQ. ZERO)) FLAG = -33
         ENDIF
         IF ((B1 .EQ. ZERO) .AND. (B2 .EQ. ZERO)) FLAG = -31
      ENDIF
      IF (FLAG .LT. 0) GOTO 120
      IF (JOB(1) .OR. JOB(2)) THEN
         DO 10 K = 1,NEV
            EV(K) = ZERO
   10    CONTINUE
      ENDIF
      IF (JOB(3)) THEN
         DO 15 K = 1,NUMT
            RHO(K) = ZERO
   15    CONTINUE
      ENDIF
      IF ((.NOT. JOB(4)) .OR. JOB(3)) THEN
         CALL CLASS(IPRINT,TOL1,JOBST(2),CSPEC,CEV,LASTEV,LPLC,STORE,
     &              JOB(5),HMIN,DOMESH)
         ALPHA = A2*A1P-A1*A2P
         IF ((A1P .NE. ZERO) .OR. (A2P .NE. ZERO)) THEN
            IF (ALPHA .LE. ZERO) FLAG = -32
         ELSE
            IF ((A1 .EQ. ZERO) .AND. (A2 .EQ. ZERO)) FLAG = -33
         ENDIF
         IF ((B1 .EQ. ZERO) .AND. (B2 .EQ. ZERO)) FLAG = -31
         IF (FLAG .LT. 0) GOTO 120
         DO 20 K = 1,2
            TYPE(1,K) = REG(K)
            TYPE(2,K) = LC(K)
            TYPE(3,K) = .NOT. OSC(K)
            TYPE(4,K) = OSC(K)
            IF (CSPEC(K)) THEN
               TYPE(3,K) = .FALSE.
               TYPE(4,K) = .FALSE.
            ENDIF
   20    CONTINUE
         IF (IPRINT .GT. 2) 
     &   CALL DSCRIP(LC,LPLC,TYPE,REG,CSPEC,CEV,CUTOFF,LASTEV,
     &               A1,A1P,A2,A2P,B1,B2)
      ELSE
         LNF = .FALSE.
         KCLASS(1) = 0
         KCLASS(2) = 0
         DO 25 K = 1,2
            REG(K) = TYPE(1,K)
            LC(K) = TYPE(2,K)
            OSC(K) = TYPE(4,K)
            CSPEC(K) = .NOT. (TYPE(3,K) .OR. TYPE(4,K))
   25    CONTINUE
         IF (.NOT. AFIN) STORE(1) = -99999.0
         IF (.NOT. BFIN) STORE(NXINIT) = 99999.0
      ENDIF
C
C     Use NSGNF to hold the sign of F when EV is large negative.
C
      SGN = A2P*B2
      IF (SGN .NE. ZERO) THEN
         NSGNF = SIGN(ONE,SGN)
      ELSE
         SGN = A1P*B2+A2P*B1
         IF (SGN .NE. ZERO) THEN
            NSGNF = SIGN(ONE,SGN)
         ELSE
            SGN = A1P*B1+A2*B2
            IF (SGN .NE. ZERO) THEN
               NSGNF = SIGN(ONE,SGN)
            ELSE
               SGN = A1*B2+A2*B1
               IF (SGN .NE. ZERO) THEN
                  NSGNF = SIGN(ONE,SGN)
               ELSE
                  NSGNF = SIGN(ONE,A1*B1)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      OSCILL = ((.NOT. TYPE(1,1)) .AND. TYPE(4,1)) .OR.
     &         ((.NOT. TYPE(1,2)) .AND. TYPE(4,2))
      TOL1 = TOL(1)+TOL(2)
      IF (JOB(1) .OR. JOB(2)) THEN
C
C        Set up approximating regular problems for eigenvalues.
C
         IF (OSCILL) THEN
            IF (IPRINT .GE. 1) THEN
               WRITE (*,30)
               WRITE (21,30)
   30          FORMAT(' This problem is oscillatory, you must use ',
     &                'interval truncation.')
            ENDIF
            FLAG = -20
            GOTO 120
         ENDIF
         IF (DOMESH) THEN
C
C           Calculate the initial mesh.
C
            K = NXINIT+16
            CALL MESH(JOB(5),-1,STORE,STORE(K),STORE(2*K+1),
     &                STORE(3*K+1),STORE(4*K+1),TOL1,HMIN)
            IF (FLAG .LT. 0) GOTO 120
         ENDIF
         IF (((KCLASS(1) .EQ. 3) .OR. (KCLASS(2) .EQ. 3)) .AND. JOB(5))
     &      LMESH = .TRUE.
         IF ((.NOT. LMESH) .AND. (IPRINT .GE. 1)) THEN
            WRITE (*,35) (STORE(I),I=1,NXINIT)
            WRITE (21,35) (STORE(I),I=1,NXINIT)
  35        FORMAT(' Level 0 mesh:',/,(5G15.6))
         ENDIF
         IF (JOB(5)) NUMX = NXINIT
C
C        Set MAXLVL, the maximum number of levels (mesh bisections).
C
C        IMPORTANT NOTE: the size of various fixed arrays in this
C        package depends on the value of MAXLVL in this FORTRAN77
C        implementation.  If MAXLVL is increased, then more storage
C        may have to be allocated to these arrays.  In particular,
C        check RATIO(*), R(*,*), and W(*,*) in EXTRAP; EVEXT(*)
C        in REGULR.
C
         MAXLVL = 10
C

         DO 45 K = 1,NEV
            EV(K) = ZERO
            IFLAG(K) = 0
            FLAG = 0
            CALL REGULR(JOB(2),LMESH,TOL,INVEC(3+K),EV(K),IPRINT,
     &                  NEXTRP,XEF,EF(1+NUMX*(K-1)),
     &                  PDEF(1+NUMX*(K-1)),HMIN,STORE)
            IF ((CSPEC(1) .OR. CSPEC(2)) .AND. (IPRINT .GE. 1) .AND.
     &          (.NOT. JOB(4)) .AND. (FLAG .GT. -5)) THEN
               IF ((EV(K) .GE. CUTOFF) .OR. ((LASTEV .NE. -5) .AND.  
     &             (INVEC(3+K) .GE. LASTEV))) THEN
                  WRITE (*,40) INVEC(3+K)
                  WRITE (21,40) INVEC(3+K)
   40             FORMAT(' WARNING: Requested eigenvalue ',I6,
     &                   ' may not be below the continuous spectrum.')
               ENDIF
            ENDIF
            IF (LFLAG(1)) THEN
               IFLAG(K) = IFLAG(K)+IBASE
               LFLAG(1) = .FALSE.
               LBASE = .TRUE.
            ENDIF
            IF (FLAG .LT. 0) IFLAG(K) = FLAG
   45    CONTINUE
         IF (LBASE) THEN
            IBASE = 10*IBASE
            LBASE = .FALSE.
         ENDIF
      ENDIF
      IF (JOB(3)) THEN
         IF (CSPEC(1) .AND. CSPEC(2)) THEN
            IFLAG(1) = -5
            IF (IPRINT .GT. 0) WRITE (*,50)
   50       FORMAT(' This problem has continuous spectrum generated by',
     &' both endpoints.  The',/,' calculation of the spectral density',
     &' function has not yet been implemented',/,' for such cases.',/)
            GOTO 120
         ENDIF
         IF (.NOT. (CSPEC(1) .OR. CSPEC(2))) THEN
            IFLAG(1) = -4
            IF (IPRINT .GT. 0) WRITE (*,55)
   55       FORMAT(' This problem has no continuous spectrum.')
            GOTO 120
         ENDIF
         IF ((CSPEC(1) .AND. ((KCLASS(2) .EQ. 5) .OR.
     &       (KCLASS(2) .EQ. 9))) .OR. (CSPEC(2) .AND.
     &       ((KCLASS(1) .EQ. 5) .OR. (KCLASS(1) .EQ. 9)))) THEN
            IFLAG(1) = -6
            IF (IPRINT .GT. 0) WRITE (*,60)
   60       FORMAT(' The normalization of the spectral density function'
     &            ,' is unknown for this problem.')
            GOTO 120
         ENDIF
         IF ((CSPEC(1) .AND. (.NOT. BFIN)) .OR. (CSPEC(2) .AND.
     &       (.NOT. AFIN))) THEN
            IFLAG(1) = -6
            IF (IPRINT .GT. 0) WRITE (*,60)
            GOTO 120
         ENDIF
         IF (OSCILL) THEN
            FLAG = -25
            IF (IPRINT .GT. 0) THEN
               WRITE (*,30)
               WRITE (21,30)
            ENDIF
            GOTO 120
         ENDIF
         XTOL = -LOG10(MAX(TOL(1),TOL(2)))
         JTOL = XTOL-HALF
         JTOL = MIN(MAX(JTOL,1),5)
         DENSOP = 3*JTOL
         MAXITS = (15-JTOL)/3
C
C        Set Maxlvl for the density function calculation; see above
C        "IMPORTANT NOTE" if this is to be increased.
C
         MAXLVL = (7+JTOL)/2
         AAFIN = AFIN
         AA = A
         KCL1 = KCLASS(1)
         BBFIN = BFIN
         BB = B
         KCL2 = KCLASS(2)
         IF (JOB(5)) THEN
C
C           Use interval truncation in this oscillatory regime.
C
            OSCILL = .FALSE.
            IF ((.NOT. JOB(4)) .AND. ((KCLASS(1) .EQ. 1) .OR.
     &          (KCLASS(2) .EQ. 1))) OSCILL = .TRUE.
            IF (.NOT. OSCILL) THEN
               NXINIT = 4*JTOL+5
               ENDFAC = ENDI(JTOL)
            ELSE
               NXINIT = 24*JTOL+36
               ENDFAC = 48.0
            ENDIF
            IF (CSPEC(1)) THEN
               IF (AFIN) THEN
                  KCLASS(1) = 7
                  ENDFAC = 4.0*ENDFAC
                  IF (BFIN) THEN
                     A = AA+(BB-AA)/ENDFAC
                  ELSE
                     A = AA+ABS(AA)/ENDFAC
                  ENDIF 
               ELSE
                  AFIN = .TRUE.
                  KCLASS(1) = 0
                  IF (BFIN) THEN
                     A = -ENDFAC-MIN(-B,ZERO)
                  ELSE
                     A = -ENDFAC
                  ENDIF
               ENDIF
            ELSE
               IF (BFIN) THEN
                  KCLASS(2) = 7
                  ENDFAC = 4.0*ENDFAC
                  IF (AFIN) THEN
                     B = BB-(BB-AA)/ENDFAC
                  ELSE
                     B = BB-ABS(BB)/ENDFAC
                  ENDIF
               ELSE
                  BFIN = .TRUE.
                  KCLASS(2) = 0
                  IF (AFIN) THEN
                     B = ENDFAC+MAX(A,ZERO)
                  ELSE
                     B = ENDFAC
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            IF (CSPEC(1)) AFIN = .TRUE.
            IF (CSPEC(2)) BFIN = .TRUE.
            IF (NUMX .EQ. 0) THEN
               IFLAG(1) = -30
               RETURN
            ENDIF
         ENDIF
         MAXT = NUMT
C
C        Loop over the choices of intervals.
C
         NUMEV = 0
         DO 105 K = 1,MAXITS
            STORE(1) = A
            STORE(NXINIT) = B
            LFLAG(3) = .FALSE.
            FLAG = 0
            IF (IPRINT .GE. 1) THEN
               WRITE (*,65) K
               WRITE (21,65) K
   65          FORMAT(60('-'),/,' Iteration ',I2)
               WRITE (21,70) A,B,NXINIT
   70          FORMAT(/,' For a, b =',2F15.8,/,' Nxinit = ',I4,/)
            ENDIF
            IF (JOB(5)) THEN
               I = NXINIT+16
               CALL MESH(.TRUE.,-1,STORE,STORE(I+1),STORE(2*I+1),
     &                   STORE(3*I+1),STORE(4*I+1),TOL1,HMIN)
            ENDIF
            IF (IPRINT .GE. 3) THEN
               WRITE (*,75) (STORE(I),I=1,NXINIT)
               WRITE (21,75) (STORE(I),I=1,NXINIT)
   75          FORMAT(' Level 0 mesh:',/,(5G15.6))
            ENDIF
            CALL DENSEF(TOL,CSPEC,IPRINT,K,NEXTRP,MAXT,T,RHO,IEV,
     &                  HMIN,NUMEV,STORE)
            IF (FLAG .EQ. -3) THEN
               LFLAG(6) = .FALSE.
               FLAG = 0
            ENDIF
            IF (FLAG .LT. 0) GOTO 120
            IF (.NOT. JOB(5)) GOTO 110
            IF (K .GT. 1) THEN
               DONE = .TRUE.
               J = MAXT
               DO 80 I = 1,J
                  RHOTOL = TWO*ZETA*MAX(TOL(1),TOL(2)*RHO(I))
                  ERROR = RHO(I)-STORE(2320+(MAXLVL+2)*NUMT+I)
                  IF (ABS(ERROR) .LE. RHOTOL) THEN
                     EDONE = .TRUE.
                  ELSE
                     EDONE = .FALSE.
                     MAXT = I
                  ENDIF
                  DONE = DONE .AND. EDONE
   80          CONTINUE
               IF (DONE) GOTO 110
            ENDIF
            IF (IPRINT .GE. 2) THEN
               WRITE (*,85)
               WRITE (21,85)
   85          FORMAT(9X,'t',15X,'Truncated Rho(t)')
               DO 95 I = 1,NUMT
                  WRITE (*,90) T(I),RHO(I)
                  WRITE (21,90) T(I),RHO(I)
   90             FORMAT(F12.4,D31.15)
   95          CONTINUE
            ENDIF
            DO 100 I = 1,MAXT
               STORE(2320+(MAXLVL+2)*NUMT+I) = RHO(I)
  100       CONTINUE
            COUNTZ = .TRUE.
            CALL SHOOT(CUTOFF,STORE,MU1,FZ)
            CALL SHOOT(T(NUMT),STORE,MU2,FZ)
            COUNTZ = .FALSE.
            IF (T(NUMT) .GT. CUTOFF) THEN
               DENS = (MU2-MU1)/(K*(T(NUMT)-CUTOFF))
            ELSE
               DENS = DENSOP
            ENDIF
            IF (.NOT. OSCILL) THEN 
               NXINIT = NXINIT+10
               ZETA = ZETAI(JTOL)
            ELSE
               ZETA = 2.0
               NXINIT = ZETA*NXINIT
            ENDIF
            IF (CSPEC(1)) THEN
               IF (AAFIN) THEN
                  ENDFAC = 5.0*ZETA
                  IF (DENS .LT. DENSLO) ENDFAC = 75.0
                  IF ((DENS .GT. DENSHI) .AND. (.NOT. OSCILL))
     &               ENDFAC = 8.0
                  A = AA+(A-AA)/ENDFAC
                  IF ((AA-A)**2 .LT. U) THEN
                     FLAG = -9
                     GOTO 110
                  ENDIF
               ELSE
                  IF (DENS .LT. DENSLO) ZETA = 2.0
                  IF ((DENS .GT. DENSHI) .AND. (.NOT. OSCILL))ZETA = 1.4
                  A = ZETA*A
               ENDIF
            ENDIF
            IF (CSPEC(2)) THEN
               IF (BBFIN) THEN
                  ENDFAC = 5.0*ZETA
                  IF (DENS .LT. DENSLO) ENDFAC = 75.0
                  IF ((DENS .GT. DENSHI) .AND. (.NOT. OSCILL))
     &               ENDFAC = 8.0
                  B = BB-(BB-B)/ENDFAC
                  IF ((B-BB)**2 .LT. U) THEN
                     FLAG = -9
                     GOTO 110
                  ENDIF
               ELSE
                  IF (DENS .LT. DENSLO) ZETA = 2.0
                  IF ((DENS .GT. DENSHI) .AND. (.NOT. OSCILL))ZETA = 1.4
                  B = ZETA*B
               ENDIF
            ENDIF
            IF (MOD(NXINIT,2) .EQ. 0) NXINIT = NXINIT+1
            NXINIT = MIN(464,NXINIT)
  105    CONTINUE
         FLAG = -3
  110    IF (IPRINT .GE. 1) WRITE (21,115) NUMEV
  115    FORMAT(' The total number of eigenvalues computed was ',I10)
         IF (CSPEC(1)) THEN
            IF (AAFIN) THEN
               A = AA
               KCLASS(1) = KCL1
            ELSE
               AFIN = .FALSE.
            ENDIF
         ENDIF
         IF (CSPEC(2)) THEN
            IF (BBFIN) THEN
               B = BB
               KCLASS(2) = KCL2
            ELSE
               BFIN = .FALSE.
            ENDIF
         ENDIF
      ENDIF
C
C     Set fatal output flags.
C
  120 IF (FLAG .LT. -9) THEN
         DO 125 K = 1,MAX(NEV,1)
            IFLAG(K) = FLAG
  125    CONTINUE
         RETURN
      ELSE
         IF ((FLAG .LT. 0) .AND. (.NOT. (JOB(1) .OR. JOB(2))))
     &      IFLAG(1) = FLAG
      ENDIF
C
C     Set warning flags.
C
      DO 135 I = 2,5
         DO 130 K = 1,MAX(NEV,1)
            IF (LFLAG(I) .AND. (IFLAG(K) .GE. 0)) IFLAG(K) = 
     &                                            IFLAG(K)+I*IBASE
  130    CONTINUE
         IF (LFLAG(I)) IBASE = 10*IBASE
  135 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE INTERV(FIRST,ALPHA,BETA,CONS,ENDFIN,NFIRST,NTOTAL,
     &                  IFLAG,X)
      LOGICAL FIRST,ENDFIN(*)
      INTEGER NFIRST,NTOTAL,IFLAG
      DOUBLE PRECISION ALPHA,BETA,CONS(*),X(*)
***********************************************************************
*                                                                     *
*     INTERV calculates the indices of eigenvalues found in a         *
*     specified interval.                                             *
*                                                                     *
***********************************************************************
C     Local variables:
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &  IDUMMY(1),I1,I2,I3,J1,J2,J3,K,LASTEV,LEVEL0,MU,NEXTRP,NLAST,NUMX
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &        CEV(2),CUTOFF,HMIN,TOL(6),V
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2),
     &        CSPEC(2),DOMESH,JOBST(3),LPLC(2)
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      SAVE CUTOFF
C
      NFIRST = -5
      NLAST = -5
      IFLAG = 0
      IF (ALPHA .GE. BETA) THEN
         IFLAG = -9
         RETURN
      ENDIF
      TOL(1) = 0.000001
      TOL(2) = 0.000001
      IF (FIRST) THEN
         JOBST(1) = .FALSE.
         JOBST(2) = .FALSE.
         JOBST(3) = .FALSE.
         AFIN = ENDFIN(1)
         BFIN = ENDFIN(2)
         MAXLVL = 10
         NUMX = 0
         CALL START(JOBST,CONS,TOL,0,IDUMMY,NUMX,X,0,X,NEXTRP,X)
         CALL CLASS(0,TOL(1),JOBST(2),CSPEC,CEV,LASTEV,LPLC,X,
     &              .TRUE.,HMIN,DOMESH)
         CUTOFF = MIN(CEV(1),CEV(2))
         IF (FLAG .LT. 0) THEN
            IFLAG = FLAG
            RETURN
         ENDIF
         IF (OSC(1) .OR. OSC(2)) THEN
            IFLAG = -25
            RETURN
         ENDIF
         IF (DOMESH) THEN
            K = NUMX+16
            CALL MESH(.TRUE.,-1,X,X(K),X(2*K+1),X(3*K+1),X(4*K+1),
     &                TOL(1),HMIN)
         ENDIF
      ENDIF
      IF (FLAG .LT. 0) THEN
         IFLAG = FLAG
         RETURN
      ENDIF
      LEVEL0 = LEVEL
      COUNTZ = .TRUE.
      LEVEL = 3
      CALL SHOOT(ALPHA,X,MU,V)
      I1 = MU
      CALL SHOOT(BETA,X,MU,V)
      J1 = MU
      LEVEL = LEVEL+1
      CALL SHOOT(ALPHA,X,MU,V)
      I2 = MU
      CALL SHOOT(BETA,X,MU,V)
      J2 = MU
   10 LEVEL = LEVEL+1
      IF (NFIRST .EQ. -5) THEN
         CALL SHOOT(ALPHA,X,MU,V)
         I3 = MU
         IF ((I1 .EQ. I2) .AND. (I2 .EQ. I3)) THEN
            NFIRST = I1
            GOTO 15
         ENDIF
         I1 = I2
         I2 = I3
      ENDIF
   15 IF (NLAST .EQ. -5) THEN
         CALL SHOOT(BETA,X,MU,V)
         J3 = MU
         IF ((J1 .EQ. J2) .AND. (J2 .EQ. J3)) THEN
            NLAST = J1-1
            IF (NFIRST .NE. -5) GOTO 20
         ENDIF
         J1 = J2
         J2 = J3
      ENDIF
      IF (LEVEL .LT. MAXLVL) GOTO 10
   20 IF (NFIRST .EQ. -5) THEN
         NFIRST = I3
         IFLAG = 12
      ENDIF
      IF (NLAST .EQ. -5) THEN
         NLAST = J3
         IFLAG = 12
      ENDIF
      NTOTAL = NLAST+1-NFIRST
      IF (NTOTAL .EQ. 0) IFLAG = 11
      IF (BETA .GT. CUTOFF) IFLAG = 13
      COUNTZ = .FALSE.
      LEVEL = LEVEL0
      RETURN
      END

C///////////////////////////////////////////////////////////////////////
      SUBROUTINE AITKEN(XLIM,TOL,N,X,ERROR)
      DOUBLE PRECISION XLIM,TOL,X(*),ERROR
      INTEGER N
C
C     Use Aitken's algorithm to accelerate convergence of the sequence
C     in X(*).
C
      DOUBLE PRECISION DENOM,XOLD,ZERO,ONE,TWO
      INTEGER I
      PARAMETER (ZERO = 0.0, ONE = 1.0, TWO = 2.0)
C
      IF (N .LE. 2) THEN
         XLIM = X(N)
         ERROR = ZERO
         RETURN
      ENDIF
      XOLD = 1.D30
      DO 10 I = 1,N-2
         DENOM = X(I+2)-TWO*X(I+1)+X(I)
         IF (DENOM .NE. ZERO) THEN
            XLIM = X(I)-(X(I+1)-X(I))**2/DENOM
            ERROR = XLIM-XOLD
         ELSE
            ERROR = X(I+2)-X(I+1)
            XLIM = X(I+2)
         ENDIF
         IF (ABS(ERROR) .LT. MAX(ONE,ABS(XLIM))*TOL) RETURN
         XOLD = XLIM
   10 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION ASYMEV(NEV,QINT,RPINT,ALPHA1,ALPHA2,
     &                                 BETA1,BETA2)
      INTEGER NEV
      DOUBLE PRECISION QINT,RPINT,ALPHA1,ALPHA2,BETA1,BETA2	
C
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2)
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &                 FNEV,
     &                 ZERO,HALF,ONE,TWO,PI
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, HALF = 0.5D0, ONE = 1.0, TWO = 2.0,
     &           PI = 3.14159265358979324D0)
C
C     Evaluate the asymptotic formula for eigenvalue NEV.  
C        Note: not all cases have been implemented yet.
C
      ASYMEV = -999999.0
      FNEV = NEV
      IF (REG(1)) THEN
         IF ((A1P .NE. ZERO) .OR. (A2P .NE. ZERO)) THEN
            IF (A2P .NE. ZERO) THEN
               ASYMEV = (TWO*A1P/A2P+QINT)/RPINT
               IF (B2 .NE. ZERO) THEN
C                 Case 1
                  ASYMEV = ASYMEV+TWO*B1/(B2*RPINT)+((FNEV-ONE)*PI/
     &                                               RPINT)**2
               ELSE
C                 Case 2
                  ASYMEV = ASYMEV+((FNEV-HALF)*PI/RPINT)**2
               ENDIF
            ELSE
               ASYMEV = (TWO*A2/A1P+QINT)/RPINT
               IF (B2 .NE. ZERO) THEN
C                 Case 3
                  ASYMEV = ASYMEV+TWO*B1/(B2*RPINT)+((FNEV-HALF)*PI/
     &                                                RPINT)**2
               ELSE
C                 Case 4
                  ASYMEV = ASYMEV+(FNEV*PI/RPINT)**2
               ENDIF
            ENDIF
         ELSE
            IF (A2 .NE. ZERO) THEN
               IF (B2 .NE. ZERO) THEN
C                 Case 1
                  ASYMEV = (FNEV*PI/RPINT)**2+(TWO*(BETA1/BETA2+
     & 	                              ALPHA1/ALPHA2)+QINT)/RPINT
               ELSE
C                 Case 2   (Dirichlet at B)
                  ASYMEV = ((FNEV+HALF)*PI/RPINT)**2+(TWO*ALPHA1/ALPHA2
     &                     +QINT)/RPINT
               ENDIF
            ELSE
               IF (B2 .NE. ZERO) THEN
C                 Case 3   (Dirichlet at A)
                  ASYMEV = ((FNEV+HALF)*PI/RPINT)**2+(TWO*BETA1/BETA2
     &                     +QINT)/RPINT
               ELSE
C                 Case 4   (Dirichlet at A and at B)
                  ASYMEV = ((FNEV+ONE)*PI/RPINT)**2+QINT/RPINT
               ENDIF
            ENDIF
         ENDIF
         RETURN
      ENDIF
      IF (REG(2)) THEN
         IF (B2 .NE. ZERO) THEN
            IF (A2 .NE. ZERO) THEN
C              Case 1
               ASYMEV = (FNEV*PI/RPINT)**2+(TWO*(ALPHA1/ALPHA2+
     &                             BETA1/BETA2)+QINT)/RPINT
            ELSE
C              Case 2   (Dirichlet at A)
               ASYMEV = ((FNEV+HALF)*PI/RPINT)**2+(TWO*BETA1/BETA2
     &                  +QINT)/RPINT
            ENDIF
         ELSE
            IF (B2 .NE. ZERO) THEN
C              Case 3   (Dirichlet at B)
               ASYMEV = ((FNEV+HALF)*PI/RPINT)**2+(TWO*ALPHA1/ALPHA2
     &                  +QINT)/RPINT
            ELSE
C              Case 4   (Dirichlet at A and at B)
               ASYMEV = ((FNEV+ONE)*PI/RPINT)**2+QINT/RPINT
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
C----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION ASYMR(NEV,RPINT,RPATA,RPATB,SCALE)
      INTEGER NEV
      DOUBLE PRECISION RPINT,RPATA,RPATB,SCALE
C
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2)
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &                 FNEV,
     &                 ZERO,HALF,ONE,TWO,PI
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, HALF = 0.5D0, ONE = 1.0, TWO = 2.0,
     &           PI = 3.14159265358979324D0)
C
C     Evaluate the asymptotic formula for RsubNEV.
C        Note: not all cases have been implemented yet.  See the note
C        above in ASYMEV.
C
      ASYMR = ZERO
      FNEV = MAX(NEV,2)
      IF (REG(1)) THEN
         IF ((A1P .NE. ZERO) .OR. (A2P .NE. ZERO)) THEN
            IF (A2P .NE. ZERO) THEN
               IF (B2 .NE. ZERO) THEN
                  ASYMR = TWO*RPINT**3/(RPATA*A2P**2*((FNEV-ONE)*PI)**4)
               ELSE
                 ASYMR = TWO*RPINT**3/(RPATA*A2P**2*((FNEV-HALF)*PI)**4)
               ENDIF
            ELSE
               IF (B2 .NE. ZERO) THEN
                  ASYMR = RPATA*TWO*RPINT/(A1P*(FNEV-HALF)*PI)**2
               ELSE
                  ASYMR = RPATA*TWO*RPINT/(A1P*FNEV*PI)**2
               ENDIF
            ENDIF
         ELSE
            IF (A2 .NE. ZERO) THEN
               ASYMR = TWO/(RPATA*A2*A2*RPINT)
            ELSE
               IF (B2 .NE. ZERO) THEN
                  ASYMR = TWO*RPATA*((FNEV+HALF)*PI/A1)**2/RPINT**3
               ELSE
                  ASYMR = TWO*RPATA*((FNEV+ONE)*PI/A1)**2/RPINT**3
               ENDIF
            ENDIF
         ENDIF
         RETURN
      ENDIF
      IF (REG(2)) THEN
         IF (A2 .NE. ZERO) THEN
            ASYMR = TWO/(RPATB*B2*B2*RPINT)
         ELSE
            IF (B2 .NE. ZERO) THEN
               ASYMR = TWO*RPATB*((FNEV+HALF)*PI/B1)**2/RPINT**3
            ELSE
               ASYMR = TWO*RPATB*((FNEV+ONE)*PI/B1)**2/RPINT**3
            ENDIF
         ENDIF
      ENDIF
      ASYMR = ASYMR*SCALE**2
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE BRCKET(N,EVLOW,EVHIGH,FLOW,FHIGH,ABSERR,RELERR,X)
      INTEGER N
      DOUBLE PRECISION EVLOW,EVHIGH,FLOW,FHIGH,ABSERR,RELERR,X(*)
C
C     Find values for EVLOW and EVHIGH which bracket the Nth eigenvalue;
C     in particular,
C           EV(N-1) < EVLOW < EV(N) < EVHIGH < EV(N+1)  .
C     It is assumed that if U(X,LAMBDA) has NZ zeros in (A,B) then
C           EV(MU-1) < LAMBDA < EV(MU)
C     where MU is a function of NZ, LAMBDA, and the constants in the
C     boundary conditions.  The value of MU for a given LAMBDA is
C     returned by a call to subprogram SHOOT.
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &        K,MU
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &                 DIFF,EV,EVSIGN,FEV,
     &                 ZERO,TWO
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2),
     &        LOW,HIGH
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, TWO = 2.0)
C
C     Set COUNTZ so that zeros are counted in SHOOT.
C
      COUNTZ = .TRUE.
      EVSIGN = NSGNF
C
C     SHOOT with Ev = Evlow should return FEV having sign EVSIGN.
C
      IF (N .NE. 2*(N/2)) EVSIGN = -EVSIGN
      LOW = .FALSE.
      HIGH = .FALSE.
C
C     Make EVLOW a lower bound for EV(N).
C
      EV = EVLOW
      DIFF = ABS(EVHIGH-EVLOW)
      IF (DIFF .EQ. ZERO) DIFF = ABSERR+RELERR
   10 CALL SHOOT(EV,X,MU,FEV)
      IF (FLAG .LT. 0) THEN
         COUNTZ = .FALSE.
         RETURN
      ENDIF
      IF (MU .GT. N) THEN
         EVHIGH = EV
         FHIGH = FEV
         EV = EV-DIFF
         DIFF = TWO*DIFF
         IF ((MU .EQ. N+1) .AND. (EVSIGN*FEV .LE. ZERO)) HIGH = .TRUE.
         GOTO 10
      ELSE
         EVLOW = EV
         FLOW = FEV
         IF ((MU .EQ. N) .AND. (EVSIGN*FEV .GE. ZERO)) LOW = .TRUE.
      ENDIF
C
C     Make EVHIGH an upper bound for EV(N).
C
      IF (.NOT. HIGH) THEN
         EV = EVHIGH
         DIFF = ABS(EVHIGH-EVLOW)
   20    CALL SHOOT(EV,X,MU,FEV)
         IF (FLAG .LT. 0) RETURN
         K = NSGNF*(-1)**MOD(MU,2)
         IF (K*FEV .LT. ZERO) THEN
            EV = EV+DIFF
            DIFF = TWO*DIFF
            GOTO 20
         ELSE
            IF (MU .LE. N) THEN
               EVLOW = EV
               FLOW = FEV
               EV = EV+DIFF
               DIFF = TWO*DIFF
               GOTO 20
            ELSE
               EVHIGH = EV
               FHIGH = FEV
               IF (MU .EQ. N+1) HIGH = .TRUE.
            ENDIF
         ENDIF
      ENDIF
C
C     Refine the interval [EVLOW,EVHIGH] to include only the Nth
C     eigenvalue.
C
   30 IF ((.NOT. LOW) .OR. (.NOT. HIGH)) THEN
         DIFF = EVHIGH-EVLOW
         EV = EVLOW+DIFF/TWO
C
C        Check for a cluster of eigenvalues within user's tolerance.
C
         IF  (TWO*DIFF .LT. MAX(ABSERR,RELERR*(MAX(ABS(EVLOW),
     &                          ABS(EVHIGH))))) THEN
            LFLAG(1) = .TRUE.
            COUNTZ = .FALSE.
            RETURN
         ENDIF
         CALL SHOOT(EV,X,MU,FEV)
         IF (FLAG .LT. 0) THEN
            COUNTZ = .FALSE.
            RETURN
         ENDIF
C
C        Update EVLOW and EVHIGH.
C
         IF (MU .EQ. N) THEN
             EVLOW = EV
             LOW = .TRUE.
             FLOW = FEV
         ELSE
            IF (MU .EQ. N+1) THEN
               EVHIGH = EV
               HIGH = .TRUE.
               FHIGH = FEV
            ELSE
               IF (MU .LT. N) THEN
                  EVLOW = EV
                  FLOW = FEV
               ELSE
                  EVHIGH = EV
                  FHIGH = FEV
               ENDIF
            ENDIF
         ENDIF
         GOTO 30
      ENDIF
      COUNTZ = .FALSE.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE CLASS(IPRINT,TOL,JOB,CSPEC,CEV,LASTEV,LPLC,X,
     &                 JMESH,HMIN,DOMESH)
      INTEGER IPRINT,LASTEV
      LOGICAL JOB,CSPEC(*),LPLC(*),JMESH,DOMESH
      DOUBLE PRECISION TOL,CEV(*),X(*),HMIN
C
C     This routine classifies the Sturm-Liouville problem.  Note:
C     (1) any computational algorithm must be based on a finite
C     amount of information; hence, there will always be cases that
C     any algorithm misclassifies.  In addition, some problems are
C     inherently ill-conditioned, in that a small change in the 
C     coefficients can produce a totally different classification.
C     (2)  The maximum number of points sampled for singular problems
C     is given by the variable KMAX.  By increasing this number, the
C     reliability of the classification may increase; however, the
C     computing time may also increase.  The values we have chosen
C     seem to be a reasonable balance for most problems.
C     (3) The algorithms apply standard theorems involving limits of
C     the Liouville normal form potential.  When this is not available,
C     each coefficient function is approximated by a power function
C     (c*x^r) and classified according to the properties of the
C     resulting Liouville approximation.
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &        KCLASS(2),IFLAG,J,K,KMAX,KUSED(2),LAST(3),M
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &        CP(2),CR(2),CUTOFF,D(4,2),EMU(2),EP(2),EQLNF(2),
     &        ER(2),ETA(2,2),PNU(2),
     &        BASE,BC(2),END,EV,FEV,OVER,PZ(40,2),QZ(40,2),
     &        RZ(40,2),S,SGN,Y(40),Z(40,2),
     &        ZERO,TENTH,ONE,TWO,EIGHT
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2),
     &        ENDFIN(2)
      COMMON /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, TENTH = 0.1D0, ONE = 1.0, TWO = 2.0,
     &           EIGHT = 8.0)
C
C     Sample the coefficient functions; determine if the problem as
C     given is in Liouville normal form.
C
      LNF = .TRUE.
      DOMESH = .TRUE.
      DO 40 J = 1,2
         IF (J .EQ .1) THEN
            END = A
            ENDFIN(1) = AFIN
            SGN = ONE
         ELSE
            END = B
            ENDFIN(2) = BFIN
            SGN = -ONE
         ENDIF
         K = TENTH-LOG10(TOL)
         KMAX = MIN(MAX(4*K,10),40)
         M = 0
         BASE = ONE/MIN(MAX(K,4),8)
         IF (ENDFIN(J)) THEN
            IF (END .NE. ZERO) KMAX = MIN(-INT(LOG10(U))-1,KMAX)
         ELSE
            KMAX = MIN(KMAX,20)
            BASE = EIGHT
         ENDIF
         DO 10 K = 1,KMAX
            Z(K,J) = BASE**K
  10     CONTINUE
         OVER = UNDER/TWO
         DO 30 K = 1,KMAX
            IF (ENDFIN(J)) THEN
               S = END+SGN*Z(K,J)
            ELSE
               S = -SGN*Z(K,J)
            ENDIF
            CALL COEFF(S,PZ(K,J),QZ(K,J),RZ(K,J))
            NCOEFF = NCOEFF+1
            IF ((PZ(K,J) .LE. ZERO) .OR. (RZ(K,J) .LE. ZERO)) THEN
               FLAG = -15
               RETURN
            ENDIF
            IF (LNF .AND. (K .GT. 1)) THEN
               IF (PZ(K,J) .NE. PZ(K-1,J)) LNF = .FALSE.
               IF (RZ(K,J) .NE. RZ(K-1,J)) LNF = .FALSE.
            ENDIF
            S = LOG(PZ(K,J))
            IF (ABS(S) .GT. OVER) M = 1
            S = LOG(ONE+ABS(QZ(K,J)))
            IF (ABS(S) .GT. OVER) M = 1
            S = LOG(RZ(K,J))
            IF (ABS(S) .GT. OVER) M = 1
            IF (M .NE. 0) GOTO 35
   30    CONTINUE
         K = KMAX
   35    KUSED(J) = K-1
   40 CONTINUE
      DO 50 J = 1,2
         CALL CLSEND(Z(1,J),PZ(1,J),QZ(1,J),RZ(1,J),KUSED(J),IPRINT,
     &        END,ENDFIN(J),J,TOL,CEV(J),CSPEC(J),BC,Y,LPLC,IFLAG)
         IF (.NOT. JOB) THEN
            IF ((.NOT. AFIN) .AND. (J .EQ. 1)) X(1) = -END
            IF ((.NOT. BFIN) .AND. (J .EQ. 2)) X(NXINIT) = END
         ENDIF
         IF ((CSPEC(J) .OR. (.NOT. OSC(J))) .AND. (.NOT. REG(J))) THEN
            IF (J .EQ. 1) THEN
               A1 = BC(1)
               A1P = ZERO
               A2 = BC(2)
               A2P = ZERO
            ELSE
               B1 = BC(1)
               B2 = BC(2)
            ENDIF
         ENDIF
         IF (IFLAG .EQ. 1) LFLAG(2) = .TRUE.
   50 CONTINUE
      CUTOFF = MIN(CEV(1),CEV(2))
C
C     Find the number of eigenvalues below the start of the
C     continuous spectrum.
C
      LASTEV = -5
      IF ((.NOT. OSC(1)) .AND. (.NOT. LC(2)) .AND. OSC(2)) LASTEV = 0
      IF ((.NOT. OSC(2)) .AND. (.NOT. LC(1)) .AND. OSC(1)) LASTEV = 0
      IF (OSC(1) .AND. (.NOT. LC(1)) .AND. OSC(2) .AND. LC(2))LASTEV = 0 
      IF (OSC(2) .AND. (.NOT. LC(2)) .AND. OSC(1) .AND. LC(1))LASTEV = 0 
      IF ((CSPEC(1) .AND. OSC(2)) .OR. (CSPEC(2) .AND. OSC(1))) LASTEV=0
      IF ((CSPEC(1) .AND. (.NOT. OSC(2))) .OR.
     &    (CSPEC(2) .AND. (.NOT. OSC(1))) .OR.
     &    (CSPEC(1) .AND. CSPEC(2))) THEN
         K = NXINIT+16
         CALL MESH(JMESH,-1,X,X(K),X(2*K+1),X(3*K+1),X(4*K+1),TOL,HMIN)
         DOMESH = .FALSE.
         COUNTZ = .TRUE.
         DO 75 J = 1,2
            LEVEL = 3*J
            EV = CUTOFF
            CALL SHOOT(EV,X,LAST(J),FEV)
            IF (FLAG .LT. 0) RETURN
            IF (IPRINT .GE. 5) WRITE (21,70) LEVEL,LAST(J)
   70       FORMAT(' When level = ',I2,', Ev index at cutoff is ',I12)
   75    CONTINUE
         COUNTZ = .FALSE.         
         IF (LAST(1) .GE. LAST(2)) THEN
            LASTEV = LAST(2)
            IF (LAST(1) .NE. LAST(2)) THEN
               LFLAG(2) = .TRUE.
               IF (IPRINT .GE. 3) WRITE (21,80)
   80          FORMAT(' The eigenvalue count is uncertain.')
               IF (LAST(1) .GT. 2*LAST(2)) LASTEV = 0
            ENDIF
         ELSE
            LASTEV = -5
         ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE CLSEND(Z,PZ,QZ,RZ,KMAX,IPRINT,END,ENDFIN,IEND,TOL,CEV,
     &                  CSPEC,BC,Y,LPLC,IFLAG)
      INTEGER KMAX,IPRINT,IEND,IFLAG
      LOGICAL ENDFIN,CSPEC,LPLC(*)
      DOUBLE PRECISION Z(*),PZ(*),QZ(*),RZ(*),END,TOL,CEV,BC(*),Y(*)
C
C     Iflag = 0    if reasonably certain of the classification;
C           = 1    if not sure.
C
C     Information about the nature of the problem at singular point
C     IEND is passed through the variable KCLASS(IEND):
C         KCLASS(*) = 0    normal;
C                   = 1    oscillatory coefficient function;
C                   = 2    regular, but 1/p, q, or r unbounded;
C                   = 3    infinite endpoint, Eqlnf = -1 ;
C                   = 4    finite singular endpoint, Tau unbounded,
C                          (not 8-10);
C                   = 5    not "hard", irregular; 
C                   = 6    "hard" irregular with Eta(1) < 0;
C                   = 7    finite end which generates Cspectrum;
C                   = 8    Q is unbounded (< 1/t^2) near a nonoscill-
C                          atory finite end;
C                   = 9    Q is unbounded (like 1/t^2) near a nonosc-
C                          illatory finite end;
C                   = 10   "hard", irregular, Eta(1) > 0.
C                    Note: "hard" means Tau goes to +infinity at a 
C                          finite nonoscillatory endpoint.
C         REG(*) = .True.  iff endpoint is regular.
C         LC(*)  = .True.  iff endpoint is limit circle.
C         OSC(*) = .True.  iff endpoint is oscillatory for all Ev.
C         CSPEC  = .True.  iff endpoint generates continuous spectrum.
C         LPLC(*)= .True.  iff theory yields Lp/Lc classification.
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &        KCLASS(2),I,IQLNF,K
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2),
     &        EX,EXACT,IRREG,POSC,QOSC,ROSC
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &        CPT(2),CRT(2),CUTOFF,D(4,2),EMU(2),EPT(2),EQLNF(2),
     &        ERT(2),ETA(2,2),PNU(2),
     &        CP,CQ,CR,C1,C2,C3,DELTA,EP,EQ,ER,GAMMA,SGN,TOL4,ZZ,
     &        ZERO,QUART,HALF,QUART3,ONE,TWO,FOUR
      COMMON /SLCLSS/CPT,CRT,CUTOFF,D,EMU,EPT,EQLNF,ERT,ETA,PNU,KCLASS
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, QUART = 0.25D0, HALF = 0.5D0, 
     &           QUART3 = 0.75D0, ONE = 1.0, TWO = 2.0, FOUR = 4.0)
C
      IFLAG = 0
      IF (IPRINT .GE. 3) THEN
         IF (IEND .EQ. 1) WRITE (21,5)
    5    FORMAT(/,' For endpoint A')
         IF (IEND .EQ. 2) WRITE (21,6)
    6    FORMAT(/,' For endpoint B:')
         WRITE (21,7) KMAX
    7    FORMAT('  Kmax = ',I3)
      ENDIF
      CSPEC = .FALSE.
      IRREG = .FALSE.
      LPLC(IEND) = .TRUE.
      KCLASS(IEND) = 0
      PNU(IEND) = ZERO
      CEV = ONE/U
      EX = .TRUE.
      C1 = ZERO
      C2 = ZERO
      TOL4 = FOUR*TOL
C
C     Seek monomial approximations to each coefficient function.
C
      CALL POWER(Z,QZ,KMAX,TOL,IPRINT,EQ,CQ,QOSC,EXACT,Y,IFLAG)
      IF (ABS(CQ) .LE. TOL) CQ = ZERO
      IF (CQ .EQ. ZERO) EQ = ZERO
      IF (LNF) THEN
         EP = ZERO
         CP = PZ(1)
         ER = ZERO
         CR = RZ(1)
         POSC = .FALSE.
         ROSC = .FALSE.
      ELSE
         CALL POWER(Z,PZ,KMAX,TOL,IPRINT,EP,CP,POSC,EXACT,Y,IFLAG)
         IF (ABS(CP) .LE. TOL) EP = ZERO
         EX = EX .AND. EXACT
         CALL POWER(Z,RZ,KMAX,TOL,IPRINT,ER,CR,ROSC,EXACT,Y,IFLAG)
         IF (ABS(CR) .LE. TOL) ER = ZERO
      ENDIF
      IF (POSC .OR. ROSC) THEN
         IF (ENDFIN) THEN
            REG(IEND) = .TRUE.
         ELSE
            IFLAG = 1
            IF (IPRINT .GE. 3) WRITE (21,10)
   10       FORMAT('  WARNING: p(x) or r(x) is not well-approximated by'
     &        ,' a power potential.',/,'  Classification is uncertain.')
            REG(IEND) = .FALSE.
            KCLASS(IEND) = 1
         ENDIF
         LC(IEND) = .TRUE.
         OSC(IEND) = .FALSE.
      ENDIF
      IF (QOSC) THEN
         IFLAG = 1
         IF (IPRINT .GE. 3) WRITE (21,20)
   20    FORMAT('  WARNING: q(x) is not well-approximated by a power ',
     &          'potential.',/,'  Classification is uncertain.')
         IF (ENDFIN) THEN
            REG(IEND) = .TRUE.
            LC(IEND) = .TRUE.
            OSC(IEND) = .FALSE.
         ELSE
            KCLASS(IEND) = 1
            REG(IEND) = .FALSE.
            LC(IEND) = .FALSE.
            CSPEC = .TRUE.
            OSC(IEND) = .FALSE.
            BC(1) = ONE
            BC(2) = ZERO
            CEV = QZ(KMAX-1)
            K = 40
            DELTA = (Z(KMAX)-Z(KMAX-1))/(K+1)
            DO 30 I = 0,K
               ZZ = Z(KMAX)-I*DELTA
               CALL COEFF(ZZ,CP,CQ,CR)
               NCOEFF = NCOEFF+1
               CEV = MIN(CEV,CQ)
   30       CONTINUE
            IF (ABS(CEV) .LT. TOL4) CEV = ZERO
            IF (U*ABS(CEV) .GE. ONE) CSPEC = .FALSE.
         ENDIF
         EQLNF(IEND) = ZERO
      ENDIF
      IF (IPRINT .GE. 3) THEN
         WRITE (21,40) CP,EP,CQ,EQ,CR,ER
   40    FORMAT('  Cp, Ep; Cq, Eq; Cr, Er =',/,3(2D25.12,/))
      ENDIF
      CPT(IEND) = CP
      EPT(IEND) = EP
C
C     Analyze this endpoint.
C
      IF ((EP .LT. ONE) .AND. (EQ .GT. -ONE) .AND. 
     &    (ER .GT. -ONE) .AND. ENDFIN) THEN
         REG(IEND) = .TRUE.
         LC(IEND) = .TRUE.
         OSC(IEND) = .FALSE.
         CSPEC = .FALSE.
         IF ((EP .GT. ZERO) .OR. (EQ .LT. ZERO) .OR. (ER .LT. ZERO))
     &      KCLASS(IEND) = 2
         RETURN
      ENDIF
      REG(IEND) = .FALSE.
      ETA(1,IEND) = HALF*(ER-EP+TWO)
      IF (ABS(ETA(1,IEND)) .LE. TOL) ETA(1,IEND) = ZERO
      IF (ETA(1,IEND) .NE. ZERO) THEN
         EQLNF(IEND) = (EQ-ER)/ETA(1,IEND)
         IQLNF = EQLNF(IEND)+SIGN(HALF,EQLNF(IEND))
         IF (ABS(IQLNF-EQLNF(IEND)) .LT. TOL4) EQLNF(IEND) = IQLNF
         C1 = (CQ/CR)*(ABS(ETA(1,IEND))*SQRT(CP/CR))**EQLNF(IEND)
         IF (C1 .EQ. ZERO) EQLNF(IEND) = ZERO
         C2 = (EP+ER)/(FOUR*ETA(1,IEND))
         C2 = C2*(C2-ONE)
         IF (IPRINT .GE. 5) THEN
           WRITE (21,50) C1,C2,EQLNF(IEND),ETA(1,IEND)
   50      FORMAT('  C1, C2      =',2D20.10,/,'  Eqlnf, Eta1 =',2D20.10)
         ENDIF
      ELSE
         C3 = (CP/CR)*(QUART*(EP+ER))**2
         IF (IPRINT .GE. 5) THEN
            WRITE (21,60) C3
   60       FORMAT('  C3 = ',D19.10)
         ENDIF
      ENDIF
      IF (.NOT. ENDFIN) THEN
C
C        Make an initial estimate for "infinity" (used in MESH).
C
         IF ((ER .GT. EQ) .OR. (CQ .EQ. ZERO)) THEN
            IF (ER .NE. EP) THEN
               GAMMA = ER-EP
               DELTA = CR/CP
            ELSE
               GAMMA = EQ-EP
               DELTA = ABS(CQ)/CP
               IF (DELTA .EQ. ZERO) DELTA = ONE
            ENDIF
         ELSE
            IF (EQ .GT. ER) THEN
               GAMMA = EQ-EP
               DELTA = ABS(CQ)/CP
            ELSE
               IF (ER .GT. EP) THEN
                  GAMMA = ER-EP
                  DELTA = CR/CP
               ELSE
                  GAMMA = ZERO
                  IF (ER .EQ. EP) THEN
                     DELTA = ABS(CR-CQ)/CP
                  ELSE
                     DELTA = ABS(CQ)/CP
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (GAMMA .GT. HALF) THEN
            IF (GAMMA .LT. TWO) THEN
               END = 80.0
            ELSE
               END = MIN(MAX(64.0/((TWO*GAMMA-3.0)*DELTA**(ONE/
     &                       (GAMMA+TWO))),ONE),80.D0)
               IF (GAMMA .GT. 24.0) END = 12.0
            ENDIF
         ELSE
            IF (GAMMA .LT. -HALF) THEN
               END = MAX(MIN(600.0*DELTA**(ONE/GAMMA)*5.0**GAMMA,
     &                       120.0D0),ONE)
            ELSE
               END = 12.0
            ENDIF
            IF ((GAMMA .EQ. ZERO) .AND. (CQ .NE. ZERO)) END = 40.0
         ENDIF
      ENDIF
C
C     Test for finite irregular singular points.
C
      IF (ENDFIN) THEN
         SGN = ONE
         I = ER-EP+SIGN(HALF,ER-EP)
         K = EQ-EP+SIGN(HALF,EQ-EP)
         IRREG = .TRUE.
         IF (CQ .EQ. ZERO) THEN
            IF ((I .GE. -2) .AND. (ABS(ER-EP-I) .LE. TOL4))
     &      IRREG = .FALSE.
         ELSE
            IF ((ER .LE. EQ) .AND. (I .GE. -2) .AND.
     &          (ABS(ER-EP-I) .LE. TOL4)) IRREG = .FALSE.
            IF ((ER .GT. EQ) .AND. (K .GE. -2) .AND.
     &          (ABS(EQ-EP-K) .LE. TOL4)) IRREG = .FALSE.
         ENDIF
         EMU(IEND) = HALF*(ONE-EP)
         IF (IRREG) THEN
            IF (IPRINT .GE. 3) WRITE (21,70)
   70       FORMAT('  This is an irregular singular point.')
            KCLASS(IEND) = 5
         ELSE
C
C           Compute the principal Frobenius root.
C
            IF (ETA(1,IEND) .NE. ZERO) THEN
               IF ((CQ .NE. ZERO) .AND. (ER .GT. EQ) .AND. (K .EQ. -2))
     &            THEN
                  PNU(IEND) = EMU(IEND)**2+CQ/CP
                  IF (ABS(PNU(IEND)) .LE. TOL4) PNU(IEND) = ZERO
                  IF (PNU(IEND) .GE. ZERO) THEN
                     PNU(IEND) = EMU(IEND)+SQRT(PNU(IEND))
                  ELSE
                     PNU(IEND) = -EP
                  ENDIF
               ELSE
                  PNU(IEND) = MAX(ONE-EP,ZERO)
               ENDIF
               IF (PNU(IEND) .GT. -EP) THEN
                  IF (IPRINT .GE. 5) WRITE (21,75) PNU(IEND)
   75             FORMAT('  The principal Frobenius root is ',E20.8)
               ENDIF
            ENDIF
         ENDIF
      ELSE
         SGN = -ONE
      ENDIF
      IF (SGN*ETA(1,IEND) .GT. ZERO) THEN
C
C        Carry out the Case 1 tests.
C
         K = 0
         IF (EQLNF(IEND) .LT. -TWO) THEN
            IF (CQ .LT. ZERO) K = 1
            IF (CQ .GT. ZERO) K = -1
         ENDIF
         IF (EQLNF(IEND) .EQ. -TWO) THEN
            IF (ABS(C1+C2+QUART) .LE. TOL4) THEN
               IF (IPRINT .GE. 3) WRITE (21,80)
   80          FORMAT('  WARNING: borderline nonoscillatory/oscillato',
     &                'ry classification.')
               K = -1
               IFLAG = 1
            ELSE
               IF (C1+C2 .LT. -QUART-TOL4) K = 1
               IF (C1+C2 .GT. -QUART) K = -1
            ENDIF
         ENDIF
         IF (EQLNF(IEND) .GT. -TWO) THEN
            IF (ABS(C2+QUART) .LE. TOL4) THEN
               C2 = -QUART   
               IF (IPRINT .GE. 3) WRITE (21,80)
               IFLAG = 1
            ENDIF
            IF (C2 .GE. -QUART) K = -1
         ENDIF
         IF (K .EQ. 1) THEN
            OSC(IEND) = .TRUE.
         ELSE
            IF (K .EQ. -1) THEN
               OSC(IEND) = .FALSE.
            ELSE
               IF (IPRINT .GE. 3) WRITE (21,85)
   85          FORMAT('  NO INFORMATION on osc/nonosc class.')
            ENDIF
         ENDIF
         K = 0
         IF (EQLNF(IEND) .LT. -TWO) THEN
            IF (CQ .GT. ZERO) K = -1
            IF (CQ .LT. ZERO) K = 1
         ENDIF
         IF (EQLNF(IEND) .EQ. -TWO) THEN
            IF (ABS(C1+C2-QUART3) .LE. TOL4) THEN
               K = -1
               IF (IPRINT .GE. 3) WRITE (21,90)
   90          FORMAT('  WARNING: borderline Lc/Lp classification.')
               IFLAG = 1
            ENDIF
            IF (C1+C2 .GE. QUART3) K = -1
            IF (ABS(C1+C2) .LT. QUART3-TOL4) K = 1
            IF (C1+C2 .LT. -TOL4) K = 1
         ENDIF
         IF (EQLNF(IEND) .GT. -TWO) THEN
            IF (ABS(C2-QUART3) .LE. TOL4) THEN
               K = -1
               IF (IPRINT .GE. 3) WRITE (21,90)
               IFLAG = 1
            ENDIF
            IF (C2 .GE. QUART3) K = -1
            IF (ABS(C2) .LT. QUART3-TOL4) K = 1
            IF (C2 .LT. -TOL4) K = 1
         ENDIF
         IF (K .EQ. 1) THEN
            LC(IEND) = .TRUE.
         ELSE
            IF (K .EQ. -1) THEN
               LC(IEND) = .FALSE.
            ELSE
               WRITE (21,95)
   95          FORMAT('  NO INFORMATION on Lp/Lc class.')
            ENDIF
         ENDIF
      ENDIF
      IF (SGN*ETA(1,IEND) .LT. ZERO) THEN
C
C        Carry out the Case 2 tests.
C
         K = 0
         IF ((EQLNF(IEND) .GT. ZERO) .AND. (CQ .LT. ZERO)) K = 1
         IF ((EQLNF(IEND) .GT. ZERO) .AND. (CQ .GT. ZERO)) K = -1
         IF (EQLNF(IEND) .EQ. ZERO) THEN
            K = -1
            CEV = CQ/CR
            CSPEC = .TRUE.
            IF (U*ABS(CEV) .GE. ONE) CSPEC = .FALSE.
         ENDIF
         IF (EQLNF(IEND) .LT. ZERO) THEN
            K = -1
            CEV = ZERO
            CSPEC = .TRUE.
         ENDIF
         IF (K .EQ. 1) THEN
            OSC(IEND) = .TRUE.
         ELSE
            IF (K .EQ. -1) THEN
               OSC(IEND) = .FALSE.
            ELSE
               WRITE (21,100)
  100          FORMAT('  NO INFORMATION on Osc/Nonosc class.')
            ENDIF
         ENDIF
         K = 0
         IF ((EQLNF(IEND) .GT. TWO) .AND. (CQ .GT. ZERO)) K = -1
         IF (EQLNF(IEND) .LE. TWO) K = -1
         IF ((EQLNF(IEND) .GT. TWO) .AND. (CQ .LT. ZERO)) K = 1
         IF (K .EQ. 1) THEN
            LC(IEND) = .TRUE.
         ELSE
            IF (K .EQ. -1) THEN
               LC(IEND) = .FALSE.
            ELSE
               WRITE (21,105)
  105          FORMAT('  NO INFORMATION on Lp/Lc class.')
            ENDIF
         ENDIF
      ENDIF
      IF (ETA(1,IEND) .EQ. ZERO) THEN
C
C        Carry out the Case 3 and 4 tests. 
C
         IF ((SGN*(EQ-ER) .LT. ZERO) .AND. (CQ .LT. ZERO)) THEN
            OSC(IEND) = .TRUE.
            LC(IEND) = .TRUE.
         ENDIF
         IF ((SGN*(EQ-ER) .LT. ZERO) .AND. (CQ .GT. ZERO)) THEN
            OSC(IEND) = .FALSE.
            LC(IEND) = .FALSE.
         ENDIF 
         IF (EQ .EQ. ER) THEN
            OSC(IEND) = .FALSE.
            LC(IEND) = .FALSE.
            CEV = CQ/CR+C3
            CSPEC = .TRUE.
            IF (U*ABS(CEV) .GE. ONE) CSPEC = .FALSE.
         ENDIF
         IF ((SGN*(EQ-ER) .GT. ZERO) .OR. (CQ .EQ. ZERO)) THEN
            OSC(IEND) = .FALSE.
            LC(IEND) = .FALSE.
            CEV = C3
            CSPEC = .TRUE.
            IF (U*ABS(CEV) .GE. ONE) CSPEC = .FALSE.
         ENDIF
      ENDIF
      IF (ABS(CEV) .LE. TOL4) CEV = ZERO
C
C     Calculate the Friedrichs boundary condition (if appropriate).
C
      IF (CSPEC) THEN
         BC(1) = ONE
         BC(2) = ZERO
      ENDIF
      IF ((.NOT. CSPEC) .AND. (.NOT. OSC(IEND))) THEN
         IF ((SGN*(ER+EP) .GT. ZERO) .AND. (SGN*(EQ+EP) .GT. ZERO)) THEN
            BC(1) = ZERO
            BC(2) = ONE
         ELSE
            IF ((SGN*(ER+EP) .GT. ZERO) .AND. (EQ+EP .EQ. ZERO)) THEN
               BC(1) = SQRT(CP*ABS(CQ))
               BC(2) = ONE
               IF (BC(1) .GT. ONE) THEN
                  BC(2) = ONE/BC(1)
                  BC(1) = ONE
               ENDIF
            ELSE
               IF ((SGN*(ER+EP) .LT. ZERO) .OR.
     &             (SGN*(EQ+EP) .LT. ZERO)) THEN
                  BC(1) = ONE
                  BC(2) = ZERO
               ENDIF
            ENDIF
         ENDIF            
      ENDIF
      IF (.NOT. OSC(IEND)) THEN
         IF ((.NOT. ENDFIN) .AND. (EQLNF(IEND) .EQ. -ONE)) 
     &      KCLASS(IEND) = 3
         I = ER-EP
         IF (CQ .NE. ZERO) THEN
            K = EQ-EP
            IF (CQ .GT. ZERO) THEN
               IF (K .LT. I) I = 0
            ELSE
               I = MIN(I,K)
            ENDIF
         ENDIF
         IF (ENDFIN .AND. (I .LT. 0)) THEN
            IF (IRREG) KCLASS(IEND) = 6
            IF (ETA(1,IEND) .GT. ZERO) THEN
C
C              Transform some nonoscillatory problems for which Tau
C              is unbounded near a finite singular endpoint.
C
               IF (IRREG) KCLASS(IEND) = 10
               CPT(IEND) = CP
               CRT(IEND) = CR
               EPT(IEND) = EP
               ERT(IEND) = ER
               EMU(IEND) = ZERO
               D(1,IEND) = (ETA(1,IEND)*SQRT(CP/CR))**
     &                     (ONE/ETA(1,IEND))
               D(2,IEND) = ONE/SQRT(SQRT(CP*CR*D(1,IEND)**(EP+ER)))
               IF ((EQLNF(IEND) .EQ. -TWO) .OR. (C1 .EQ. ZERO)) THEN
                  IF (.NOT. IRREG) KCLASS(IEND) = 9
                  EMU(IEND) = ABS(QUART+C1+C2)
                  IF (EMU(IEND) .LT. TOL4) THEN
                     EMU(IEND) = HALF
                  ELSE
                     EMU(IEND) = HALF+SQRT(EMU(IEND))
                  ENDIF
               ELSE
                  IF (.NOT. IRREG) KCLASS(IEND) = 8
               ENDIF
               ETA(2,IEND) = EMU(IEND)-QUART*(EP+ER)/ETA(1,IEND)
               IF ((KCLASS(IEND) .EQ. 10) .AND. (EMU(IEND) .EQ. ZERO))
     &         THEN
                  ETA(2,IEND) = HALF*(ONE-EP)/ETA(1,IEND)
                  EMU(IEND) = ETA(2,IEND)+QUART*(EP+ER)/ETA(1,IEND)
               ENDIF
               D(3,IEND) = ETA(2,IEND)*(ETA(2,IEND)+(EP-ONE)/
     &                                  ETA(1,IEND))
               D(4,IEND) = D(3,IEND)
               IF (EQLNF(IEND) .EQ. -TWO) D(4,IEND) = D(4,IEND)-C1
               IF (ABS(D(4,IEND)) .LE. TOL4) D(4,IEND) = ZERO
               D(4,IEND) = SQRT(ABS(D(4,IEND)))
               IF (IPRINT .GE. 5) THEN
                  WRITE (21,110) EMU(IEND),ETA(2,IEND)
                  WRITE (21,115) D(3,IEND),D(4,IEND)
  110             FORMAT('  Mu =',D20.6,';  Eta2 =',D20.6)
  115             FORMAT('  D3 =',D20.6,';    D4 =',D20.6)
               ENDIF
            ENDIF
            IF (KCLASS(IEND) .GE. 9) THEN
               IF (.NOT. EX) LFLAG(5) = .TRUE.
               IF (IPRINT .GE. 3) THEN
                  WRITE (21,120)
  120             FORMAT('  This problem has unbounded ',
     &                   '[Ev*r(x)-q(x)]/p(x).')
                  IF ((EMU(IEND) .GT. ZERO) .OR. (EP+ER .NE. ZERO))
     &               WRITE (21,125)
  125             FORMAT('  A change of variables will be used near',
     &                   ' this endpoint.')
               ENDIF
            ENDIF
         ENDIF
         IF (ENDFIN .AND. (.NOT. REG(IEND)) .AND. (KCLASS(IEND)
     &       .EQ. 0)) KCLASS(IEND) = 4
      ENDIF
      IF ((POSC .OR. QOSC .OR. ROSC) .AND. (.NOT. ENDFIN)) END = 99.0
      IF (IPRINT .GE. 5) THEN
         WRITE (21,130) KCLASS(IEND)
  130    FORMAT('  Classification type (KCLASS) is: ',I2)
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DENSEF(TOL,CSPEC,IPRINT,ITER,NEXTRP,NUMT,T,RHO,NEV,
     &                  HMIN,NUMEV,STORE)
      INTEGER IPRINT,ITER,NEXTRP,NUMT,NEV,NUMEV
      LOGICAL CSPEC(2)
      DOUBLE PRECISION TOL(*),T(*),RHO(*),HMIN,STORE(*)
***********************************************************************
*                                                                     *
*     This routine computes the spectral density function rho(t).     *
*                                                                     *
***********************************************************************
C
C    Input parameters:
C      TOL(*) as in SLEDGE.
C      CSPEC(*)  = logical 2-vector; CSPEC(i) = .true. iff endpoint i
C                  (1 = A, 2 = B) generates continuous spectrum.
C      IPRINT    = integer controlling printing.
C      ITER      = iteration from SLEDGE.
C      NEXTRP    = integer giving maximum no. of extrapolations.
C      NUMT      = integer equalling number of T(*) points.
C      T(*)      = real vector of abcissae for spectral function rho(t). 
C      HMIN      = minimum stepsize in Level 0 mesh.
C
C    Output parameters:
C      RHO(*)    = real vector of values for spectral density function
C                  rho(t),  RHO(I) = rho(T(I)).
C      NEV       = integer pointer to eigenvalue.  On a normal return
C                  (FLAG = 0) this is set to the index of the last 
C                  eigenvalue computed; if FLAG is not zero, then NEV
C                  gives the index of the eigenvalue where the problem
C                  occurred.
C      NUMEV     = cumulative number of eigenvalues computed.
C
C    Auxiliary storage:
C      STORE(*) = real vector of auxiliary storage, must be dimensioned
C                 at least 5*Nxinit+(Maxlvl+2)*NUMT.  The value of
C                 Nxinit is either the input NUMX or Maxint.  Currently,
C                 Maxlvl = 8 and Maxint = 235.
C            1    ->     Nxinit         vector of mesh points X(*),
C       Nxinit+1  ->    5*Nxinit        intermediate RsubN calculations,
C     5*Nxinit+1  -> 5*Nxinit+10*NUMT   intermediate RHO values.
C-----------------------------------------------------------------------
C     The definition of a spectral density function assumes a certain
C     normalization on the eigenfunctions.  For the case when x = b
C     generates the continuous spectrum, the normalization used here 
C     (and in routine GETRN below) is:
C     (1) when x = a is regular  u(a) = (A2-A2P*Ev)/SCALE
C                            (pu')(a) = (A1-A1P*Ev)/SCALE,
C         with SCALE = sqrt(A1**2+A2**2) when A1' = A2' = 0, and
C              SCALE = sqrt(ALPHA) otherwise. 
C     (2) When x = a is a regular singular point then u(x) is taken to 
C         be asymptotic to the principal Frobenius solution, i.e.,
C         near x = a   u(x) ~ (x-a)**Nu    with Nu the larger
C         root of the indicial equation.
C     Analogous normalizations hold at x = b when the endpoint x = a
C     generates the continuous spectrum.
C----------------------------------------------------------------------
C     Local variables:
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,     
     &     KCLASS(2),I,IASYMP,J,JS,KLVL,KRHO,LPRINT,MAXNEV,NRHO,NSAVE
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2),
     &        DONE,EFIN(2,2),JUMP,RDONE
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &       CP(2),CR(2),CUTOFF,D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     &       ETAT(2,2),PNU(2),
     &       ABSERR,ALPHA,ALPHA1,ALPHA2,ASYMEV,ASYMR,AEV,ARN,AVGRHO,
     &       BETA1,BETA2,BIG,DELTA,DENOM,DX,EFNORM,ERROR,ETA,ETAOLD,EV,
     &       EVHAT,EVHIGH,EVLOW,EVOLD(200),EVSAVE,FLOW,FHIGH,H,HALFH,
     &       PDU,PSRHO,PX,QINT,QLNF,QX,RELERR,RHOSUM,ROLD,RPATA,RPATB,
     &       RPINT,RSAVE(5),RSUBN,RX,SCALE,SQRTRP,TENU,TOL1,TOL2,
     &       TOLMAX,TOLMIN,UX,XLEFT,XTOL,Z,ZABS,ZREL,
     &       ZERO,HALF,ONE,TWO,THREE,FOUR,SIX,EIGHT,TEN
      COMMON /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETAT,PNU,KCLASS
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, HALF = 0.5D0, ONE = 1.0,
     &           TWO = 2.0, THREE = 3.0, FOUR = 4.0, SIX = 6.0,
     &           EIGHT = 8.0, TEN = 10.0, TOLMIN = 5.D-3)
C
C     Initialization:
C
      NSAVE = 200
      JS = 5*NXINIT
      TENU = TEN*U
      BIG = ONE/U
      AVGRHO = BIG
      LPRINT = MIN(IPRINT,3)
      MAXNEV = 0
      IF (REG(1)) THEN
         DX = SQRT(U)*MAX(ONE,ABS(A))
         CALL COEFF(A+DX,PX,QX,RX)
         IF (FLAG .LT. 0) RETURN
         RPATA = RX*PX
         ALPHA1 = ONE/RX
         CALL COEFF(A+TWO*DX,PX,QX,RX)
         ALPHA1 = ALPHA1*(RPATA-PX*RX)/(DX*FOUR)
         RPATA = SQRT(RPATA)
         ALPHA2 = SQRT(RPATA)
         ALPHA1 = (A1+A2*ALPHA1)/ALPHA2
         ALPHA2 = ALPHA2*A2
         NCOEFF = NCOEFF+2
      ELSE
         ALPHA1 = ZERO
         ALPHA2 = ONE
         RPATA = ONE
      ENDIF
      IF (REG(2)) THEN
         DX = SQRT(U)*MAX(ONE,ABS(B))
         CALL COEFF(B-DX,PX,QX,RX)
         IF (FLAG .LT. 0) RETURN
         RPATB = RX*PX
         BETA1 = ONE/RX
         CALL COEFF(B-TWO*DX,PX,QX,RX)
         BETA1 = BETA1*(RPATB-PX*RX)/(DX*FOUR)
         RPATB = SQRT(RPATB)
         BETA2 = SQRT(RPATB)
         BETA1 = (B1-B2*BETA1)/BETA2
         BETA2 = BETA2*B2
         NCOEFF = NCOEFF+2
      ELSE
         BETA1 = ZERO
         BETA2 = ONE
      ENDIF
      ALPHA = A1P*A2-A1*A2P
      IF (CSPEC(2)) THEN
         IF (ALPHA .EQ. ZERO) THEN
            SCALE = SQRT(A1**2+A2**2)
         ELSE
            SCALE = SQRT(ALPHA)
         ENDIF
      ENDIF      
      IF (CSPEC(1)) SCALE = SQRT(B1**2+B2**2)
      TOLMAX = MAX(TOL(1),TOL(2))
      TOL1 = MIN(TOL(1),TOLMIN)
      TOL2 = MIN(TOL(2),TOLMIN)
      ABSERR = TOL1
      RELERR = TOL2
      KLVL = 1
      DELTA = HALF
C
C     Begin the Main loop over LEVEL.
C
      DO 120 LEVEL = 0,MAXLVL
         IF (IPRINT .GE. 2) THEN
            WRITE (*,10) LEVEL,ITER
            WRITE (21,10) LEVEL,ITER
   10       FORMAT(' Level',I3,' of iteration',I3)
         ENDIF
         NRHO = 1
         KRHO = 0
         PSRHO = ZERO
         RHOSUM = ZERO
         ROLD = ZERO
         IASYMP = 0
         EVSAVE = -BIG
         NEV = 0
         JUMP = .FALSE.
C
C        Compute integrals needed in asymptotic formulas.
C
         QINT = ZERO
         RPINT = ZERO
         DO 20 I = 2,NXINIT
            XLEFT = STORE(I-1)
            H = (STORE(I)-XLEFT)/KLVL
            HALFH = HALF*H
            DO 15 J = 1,KLVL
               Z = XLEFT+HALFH
               XLEFT = XLEFT+H
               CALL PQRINT(Z,SQRTRP,QLNF)
               IF (FLAG .LT. 0) RETURN
               QINT = QINT+H*QLNF
               RPINT = RPINT+H*SQRTRP
   15       CONTINUE
   20    CONTINUE
         IF (QINT .GT. ONE/U) QINT = ZERO
         ZABS = MAX(MIN(ABSERR/100.0,RELERR)/10.0,TENU)
         ZREL = RELERR/10.0
         DELTA = MAX(DELTA/SIX,TOL1+TOL2)
C
C        Begin the secondary loop over NEV.
C
   25       IF (IASYMP .GE. 2) THEN
               AEV = ASYMEV(NEV,QINT,RPINT,ALPHA1,ALPHA2,BETA1,BETA2)
               EVHAT = AEV
               IF (IASYMP .LE. 3) GOTO 35
               RSUBN = ASYMR(NEV,RPINT,RPATA,RPATB,SCALE)
               GOTO 45
            ENDIF
            IF (HMIN/KLVL .LE. TENU) FLAG = -8
            IF (FLAG .LT. 0) RETURN
            IF ((LEVEL .GT. 0) .AND. (NEV .LT. MIN(MAXNEV,NSAVE))) THEN
               EV = MAX(HALF*TOL1,DELTA,HALF*TOL2*ABS(EVOLD(NEV+1)))
               EVLOW = EVOLD(NEV+1)-EV
               EVHIGH = EVOLD(NEV+1)+EV
            ELSE
               IF (LEVEL .EQ. 0) THEN
                  EV = ZERO
               ELSE
                  EV = ASYMEV(NEV,QINT,RPINT,ALPHA1,ALPHA2,BETA1,BETA2)
               ENDIF
               ETA = MAX(HALF*TOL1,DELTA,HALF*TOL2*ABS(EV))
               EVLOW = EV-ETA
               EVHIGH = EV+ETA
            ENDIF
            CALL BRCKET(NEV,EVLOW,EVHIGH,FLOW,FHIGH,ZABS,ZREL,STORE)
            IF ((LEVEL .EQ. 0) .AND. (NEV .EQ. 0)) DELTA = EVHIGH-EVLOW
            IF (FLAG .LT. 0) RETURN
            IF (ABS(EVHIGH-EVLOW) .GT. MAX(ZABS,ZREL*ABS(EVHIGH)))
     &          CALL ZZERO(EVLOW,EVHIGH,FLOW,FHIGH,ZABS,ZREL,J,STORE)
            EVHAT = MIN(EVLOW,EVHIGH)
            CALL GETEF(EVHAT,EFNORM,LPRINT,STORE,EFIN)
            IF (FLAG .LT. 0) RETURN
            IF (CSPEC(2)) THEN
               IF (REG(1) .OR. (PNU(1) .EQ. ZERO) .OR.
     &             (PNU(1) .EQ. ONE-EP(1))) THEN
                  UX = ABS(STORE(NXINIT+1))
                  PDU = ABS(STORE(2*NXINIT+1))
                  IF (A2-A2P*EVHAT .NE. ZERO) THEN
                     DENOM = SCALE*UX/ABS(A2-A2P*EVHAT)
                  ELSE
                     DENOM = SCALE*PDU/ABS(A1-A1P*EVHAT)
                  ENDIF
               ELSE
                  H = STORE(2)-STORE(1)
                  UX = ABS(STORE(NXINIT+2))
                  PDU = ABS(STORE(2*NXINIT+2))
                  IF (UX .GE. PDU) THEN
                     DENOM = UX/H**PNU(1)
                  ELSE
                     DENOM = PDU/(CP(1)*ABS(PNU(1))*
     &                           H**(EP(1)+PNU(1)-ONE))
                  ENDIF
               ENDIF
            ELSE
               IF (REG(2) .OR. (PNU(2) .EQ. ZERO) .OR.
     &             (PNU(2) .EQ. ONE-EP(2))) THEN
                  UX = ABS(STORE(2*NXINIT))
                  PDU = ABS(STORE(3*NXINIT))
                  IF (B2 .NE. ZERO) THEN
                     DENOM = SCALE*UX/ABS(B2)
                  ELSE
                     DENOM = SCALE*PDU/ABS(B1)
                  ENDIF
               ELSE
                  H = STORE(NXINIT)-STORE(NXINIT-1)
                  UX = ABS(STORE(2*NXINIT-1))
                  PDU = ABS(STORE(3*NXINIT-1))
                  IF (UX .GE. PDU) THEN
                     DENOM = UX/H**PNU(2)
                  ELSE
                     DENOM = PDU/(CP(2)*ABS(PNU(2))*
     &                            H**(EP(2)+PNU(2)-ONE))
                  ENDIF
               ENDIF
            ENDIF
            IF (BIG*DENOM .GE. ONE) THEN
               EFNORM = ONE/DENOM**2
            ELSE
               EFNORM = ONE/U**2
            ENDIF
C
C           Test for asymptotic EV.
C
            AEV = ASYMEV(NEV,QINT,RPINT,ALPHA1,ALPHA2,BETA1,BETA2) 
            IF (TEN*ABS(AEV-EVHAT) .GT. MAX(ABSERR,RELERR*ABS(EVHAT)))
     &      THEN
               IASYMP = 0
            ELSE
               IF (IASYMP .LT. 1) THEN
                  IASYMP = 1
               ELSE
                  IF (IPRINT .GE. 2) THEN
                     WRITE (*,30) NEV
                     WRITE (21,30) NEV
   30                FORMAT(' Switchover to asymptotic eigenvalues at',
     &                      ' Nev =',I8)
                  ENDIF
                  IASYMP = 2
               ENDIF
            ENDIF
            IF (EFNORM .LT. BIG) THEN
               RSUBN = ONE/(ALPHA+EFNORM)
            ELSE
               RSUBN = ZERO
            ENDIF
C 
C           Test for asymptotic RsubN.
C  
   35       ARN = ASYMR(NEV,RPINT,RPATA,RPATB,SCALE)
            IF (IASYMP .GE. 2) THEN
C 
C           Eigenvalues from asymptotic formulas; produce current RsubN.
C 
               CALL GETRN(AEV,ALPHA,CSPEC,SCALE,RSUBN,STORE)
               IF (FLAG .LT. 0) RETURN
               IF (ABS(ARN-RSUBN) .GT. MAX(ABSERR/100.0,RELERR*ARN))
     &         THEN
                  IASYMP = 2
               ELSE
                  IF (IASYMP .LT. 3) THEN
                     IASYMP = 3
                  ELSE
                     IF (IPRINT .GE. 2) THEN
                        WRITE (*,40) NEV
                        WRITE (21,40) NEV
   40                   FORMAT(' Switchover to asymptotic RsubN at ',
     &                         'Nev = ',I8)
                     ENDIF
                     IASYMP = 4
                  ENDIF
               ENDIF
            ENDIF
   45       IF (NEV .LT. NSAVE) EVOLD(NEV+1) = EVHAT
            IF (NEV .GT. 0) ETA = HALF*(EVHAT+EVSAVE)
            IF (EVHAT .LT. CUTOFF) PSRHO = PSRHO+RSUBN
   50       IF (T(NRHO) .LE. CUTOFF) THEN
C
C              Use step functions for Rho(t).
C
               IF (T(NRHO) .LE. EVHAT) THEN
                  RHO(NRHO) = RHOSUM
                  NRHO = NRHO+1
                  IF (NRHO .LE. NUMT) THEN
                     GOTO 50
                  ELSE
                     GOTO 85
                  ENDIF
               ENDIF
            ELSE
               IF (NEV .EQ. 0) GOTO 78
C
C              Use linear interpolation for Rho(t).
C
   55          IF (T(NRHO) .LE. ETA) THEN
                  IF (JUMP) THEN
   60                IF (T(NRHO) .LE. EVSAVE) THEN
                        RHO(NRHO) = RHOSUM-ROLD
                        NRHO = NRHO+1
                        IF (NRHO .LE. NUMT) THEN
                           GOTO 60
                        ELSE
                           JUMP = .FALSE.
                           GOTO 85
                        ENDIF
                     ELSE
                        JUMP = .FALSE.
   65                   IF (T(NRHO) .LE. ETA) THEN
                           RHO(NRHO) = RHOSUM
                           NRHO = NRHO+1
                           IF (NRHO .GT. NUMT) GOTO 85
                           GOTO 65
                        ELSE
                           GOTO 70
                        ENDIF
                     ENDIF
                  ENDIF
                  IF ((NEV .LE. 1) .OR. (CUTOFF .GT. EVSAVE)) THEN
                     UX = CUTOFF
                  ELSE
                     UX = ETAOLD
                  ENDIF
                  RHO(NRHO) = RHOSUM-ROLD*(ETA-T(NRHO))/(ETA-UX)
                  NRHO = NRHO+1
                  IF (NRHO .LE. NUMT) THEN
                     GOTO 55
                  ELSE
                     GOTO 85
                  ENDIF
               ENDIF
   70          IF ((RSUBN .GT. MAX(SIX*ROLD,TEN*TOLMAX,FOUR*AVGRHO))
     &             .AND. (EVHAT .GT. CUTOFF) .AND. (KRHO .GT. 4)) THEN
C
C                 Possible eigenvalue in the continuous spectrum.
C
                  LFLAG(3) = .TRUE.
                  IF (IPRINT .GE. 1) THEN
                     WRITE (*,75) EVHAT
                     WRITE (21,75) EVHAT
   75                FORMAT(' Large jump in the step spectral density',
     &                      ' function at',D17.10)
                     WRITE (*,76) ITER,LEVEL,RSUBN
                     WRITE (21,76) ITER,LEVEL,RSUBN
   76                FORMAT(18X,'Iteration =',I2,', level = ',I2,
     &                      ', jump = ',D17.10)
                  ENDIF
                  JUMP = .TRUE.
               ENDIF
            ENDIF
            IF  (NEV .GT. 0) ETAOLD = ETA
   78       RHOSUM = RHOSUM+RSUBN
            ROLD = RSUBN
            EVSAVE = EVHAT
C
C           Output requested information.
C
            IF (IPRINT .GE. 3) THEN
               IF ((NEV .LE. 25) .OR. ((IASYMP .LE. 1) .AND.
     &            (IPRINT .GE. 5))) THEN
                  WRITE (21,80) NEV,EVHAT,RSUBN
                  WRITE (*,80) NEV,EVHAT,RSUBN
   80             FORMAT(' Nev =',I7,', EvHat =',D15.6,', RHat =',D15.6)
               ELSE
                  IF ((NEV .LT. 100) .AND. (10*(NEV/10) .EQ. NEV)) THEN
                     WRITE (21,80) NEV,EVHAT,RSUBN
                     WRITE (*,80) NEV,EVHAT,RSUBN
                  ELSE
                     IF ((NEV .LT. 1000) .AND. (100*(NEV/100) .EQ. NEV))
     &                  THEN
                        WRITE (21,80) NEV,EVHAT,RSUBN
                        WRITE (*,80) NEV,EVHAT,RSUBN
                     ELSE
                        IF (1000*(NEV/1000) .EQ. NEV) THEN
                           WRITE (21,80) NEV,EVHAT,RSUBN
                           WRITE (*,80) NEV,EVHAT,RSUBN
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            NEV = NEV+1
            IF (EVHAT .GE. CUTOFF) THEN
               IF (KRHO .EQ. 5) THEN
                  RSAVE(1) = RSAVE(2)
                  RSAVE(2) = RSAVE(3)
                  RSAVE(3) = RSAVE(4)
                  RSAVE(4) = RSAVE(5)
               ELSE
                  KRHO = KRHO+1
               ENDIF
               RSAVE(KRHO) = RSUBN
               AVGRHO = ZERO
               DO 84 I = 1,KRHO
                  AVGRHO = AVGRHO+RSAVE(I)
   84          CONTINUE
               IF (KRHO .GT. 0) AVGRHO = AVGRHO/KRHO
            ENDIF
         GOTO 25
C
C        End of Nev loop ------------------
C
   85    MAXNEV = NEV
         NUMEV = NUMEV+MAXNEV
         IF (IPRINT .GE. 3) THEN
            WRITE (*,90) MAXNEV
            WRITE (21,90) MAXNEV
   90       FORMAT(' MaxNev = ',I8)
         ENDIF
         IF (IPRINT .GE. 4) THEN
            WRITE (*,95) 
            WRITE (21,95) 
   95       FORMAT(9X,'t',18X,'RhoHat(t)')
            DO 105 J = 1,NUMT
               WRITE (*,100) T(J),RHO(J)
               WRITE (21,100) T(J),RHO(J)
  100          FORMAT(F12.4,E31.15)
  105       CONTINUE
         ENDIF
C
C        Extrapolate interpolated approximations.
C
         DONE = .TRUE.
         DO 110 J = 1,NUMT
            XTOL = MAX(TOL1,ABS(RHO(J))*TOL2)
            CALL EXTRAP(RHO(J),XTOL,LEVEL+1,NEXTRP,.TRUE.,.FALSE.,1,
     &                  STORE((MAXLVL+1)*J+JS),IPRINT,ERROR,RDONE)
            IF (RHO(J) .LT. ZERO) RHO(J) = ZERO
            IF (J .GT. 1) RHO(J) = MAX(RHO(J),RHO(J-1)) 
            IF (ERROR .LE. HALF*XTOL) RDONE = .TRUE.
            DONE = RDONE .AND. DONE
  110    CONTINUE
         IF (DONE) RETURN
         ABSERR = MAX(HALF*ABSERR,TENU)
         RELERR = MAX(HALF*RELERR,TENU)
  120 CONTINUE
      FLAG = -3
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DSCRIP(LC,LPLC,TYPE,REG,CSPEC,CEV,CUTOFF,LASTEV,
     &                  A1,A1P,A2,A2P,B1,B2)
      INTEGER LASTEV
      LOGICAL LC(*),LPLC(*),TYPE(4,*),REG(*),CSPEC(*)
      DOUBLE PRECISION CEV(*),CUTOFF,A1,A1P,A2,A2P,B1,B2
C
C     Output (if requested) a description of the spectrum.
C
      WRITE (21,*)
      IF (TYPE(3,1) .AND. TYPE(3,2)) THEN
C
C     Category 1
C
         WRITE (21,123) 1
         WRITE (21,100)
         WRITE (21,121)
         WRITE (21,102)
         WRITE (21,110)
         IF (REG(1)) THEN
            WRITE (21,112)
         ELSE
            WRITE (21,113)
            WRITE (21,114)
            IF (LC(1)) THEN
               IF (LPLC(1)) WRITE (21,117)
            ELSE
               IF (LPLC(1)) WRITE (21,118)
            ENDIF
            WRITE (21,119) A1,A1P,A2,A2P
         ENDIF
         WRITE (21,111)
         IF (REG(2)) THEN
            WRITE (21,112)
         ELSE
            WRITE (21,113)
            WRITE (21,114)
            IF (LC(2)) THEN
               IF (LPLC(2)) WRITE (21,117)
            ELSE
               IF (LPLC(2)) WRITE (21,118)
            ENDIF
            WRITE (21,119) B1,B2
         ENDIF
      ENDIF
C
      IF ((TYPE(3,1) .AND. CSPEC(2)) .OR.
     &    (TYPE(3,2) .AND. CSPEC(1))) THEN
C
C        Category 2
C
         WRITE (21,123) 2
         WRITE (21,100)
         WRITE (21,103) CUTOFF
         IF (LASTEV .EQ. -5) THEN
            WRITE (21,105)
            WRITE (21,120)
         ELSE
            IF (LASTEV .EQ. 0) THEN
               WRITE (21,107)
            ELSE
               IF (LASTEV .EQ. 1) THEN
                  WRITE (21,108)
               ELSE
                  WRITE (21,109) LASTEV
               ENDIF
            ENDIF
         ENDIF
         WRITE (21,110) 
         IF (REG(1)) THEN
            WRITE (21,112)
         ELSE
            WRITE (21,113)
            IF (CSPEC(1)) THEN
               WRITE (21,116) CEV(1)
            ELSE
               WRITE (21,114)
               WRITE (21,119) A1,A1P,A2,A2P
            ENDIF
            IF (LC(1)) THEN
               IF (LPLC(1)) WRITE (21,117)
            ELSE
               IF (LPLC(1)) WRITE (21,118)
            ENDIF
         ENDIF
         WRITE (21,111)
         IF (REG(2)) THEN
            WRITE (21,112)
         ELSE
            WRITE (21,113)
            IF (CSPEC(2)) THEN
               WRITE (21,116) CEV(2)
            ELSE
               WRITE (21,114)
               WRITE (21,119) B1,B2
            ENDIF
            IF (LC(2)) THEN
               IF (LPLC(2)) WRITE (21,117)
            ELSE
               IF (LPLC(2)) WRITE (21,118)
            ENDIF
         ENDIF
      ENDIF
C
      IF ((TYPE(3,1) .AND. (TYPE(4,2) .AND. LC(2) .AND.
     &                           (.NOT. CSPEC(2)))) .OR.
     &    (TYPE(3,2) .AND. (TYPE(4,1) .AND. LC(1) .AND.
     &                           (.NOT. CSPEC(1))))) THEN
C
C        Category 3
C
         WRITE (21,123) 3
         WRITE (21,100)
         WRITE (21,106)
         WRITE (21,102)
         WRITE (21,110)
         IF (REG(1)) THEN
            WRITE (21,112)
         ELSE 
            WRITE (21,113)
            IF (TYPE(4,1)) THEN
               WRITE (21,115)
            ELSE
               WRITE (21,114)
               WRITE (21,119) A1,A1P,A2,A2P
            ENDIF
            IF (LC(1)) THEN
               IF (LPLC(1)) WRITE (21,117)
            ELSE
               IF (LPLC(1)) WRITE (21,118)
            ENDIF
         ENDIF
         WRITE (21,111)
         IF (REG(2)) THEN
            WRITE (21,112)
         ELSE 
            WRITE (21,113)
            IF (TYPE(4,2)) THEN
               WRITE (21,115)
            ELSE
               WRITE (21,114)
               WRITE (21,119) B1,B2
            ENDIF
            IF (LC(2)) THEN
               IF (LPLC(2)) WRITE (21,117)
            ELSE
               IF (LPLC(2)) WRITE (21,118)
            ENDIF
         ENDIF
      ENDIF
C
      IF ((TYPE(3,1) .AND. ((.NOT. LC(2)) .AND. TYPE(4,2) .AND.
     &     (.NOT. CSPEC(2)))) .OR.
     &    (TYPE(3,2) .AND. ((.NOT. LC(1)) .AND. TYPE(4,1) .AND.
     &     (.NOT. CSPEC(1))))) THEN
C
C        Category 4
C
         WRITE (21,123) 4
         WRITE (21,100)
         WRITE (21,104)
         WRITE (21,110)
         IF (REG(1)) THEN
            WRITE (21,112)
         ELSE
            WRITE (21,113)
            IF (TYPE(4,1)) THEN
               WRITE (21,115)
            ELSE
               WRITE (21,114)
               WRITE (21,119) A1,A1P,A2,A2P
            ENDIF
            IF (LC(1)) THEN
               IF (LPLC(1)) WRITE (21,117)
            ELSE
               IF (LPLC(1)) WRITE (21,118)
            ENDIF
         ENDIF
         WRITE (21,111)
         IF (REG(2)) THEN
            WRITE (21,112)
         ELSE
            WRITE (21,113)
            IF (TYPE(4,2)) THEN
               WRITE (21,115)
            ELSE
               WRITE (21,114)
               WRITE (21,119) B1,B2
            ENDIF
            IF (LC(2)) THEN
               IF (LPLC(2)) WRITE (21,117)
            ELSE
               IF (LPLC(2)) WRITE (21,118)
            ENDIF
         ENDIF
      ENDIF
C
      IF ((LC(1) .AND. TYPE(4,1)) .AND. (LC(2) .AND. TYPE(4,2))) THEN
C
C        Category 5
C
         WRITE (21,123) 5
         WRITE (21,100)
         WRITE (21,106)
         WRITE (21,102)
         WRITE (21,110)
         WRITE (21,113)
         WRITE (21,115)
         WRITE (21,117)
         WRITE (21,111)
         WRITE (21,113)
         WRITE (21,115)
         WRITE (21,117)
      ENDIF
C
      IF ((TYPE(4,1) .AND. (.NOT. LC(1)) .AND. (.NOT. CSPEC(1))) .AND.
     &    (TYPE(4,2) .AND. (.NOT. LC(2)) .AND. (.NOT. CSPEC(2)))) THEN
C
C        Category 6
C
         WRITE (21,123) 6
         WRITE (21,122)
         WRITE (21,110)
         WRITE (21,113)
         WRITE (21,115)
         WRITE (21,118)
         WRITE (21,111)
         WRITE (21,113)
         WRITE (21,115)
         WRITE (21,118)
      ENDIF
C
      IF ((TYPE(4,1) .AND. TYPE(4,2) .AND. .NOT. (CSPEC(1) .OR. 
     &     CSPEC(2))) .AND. ((LC(1) .AND. (.NOT. LC(2))) .OR.
     &    ((.NOT. LC(1)) .AND. LC(2)))) THEN
C
C        Category 7
C
         WRITE (21,123) 7
         WRITE (21,104)
         WRITE (21,110)
         WRITE (21,113)
         WRITE (21,115)
         IF (LC(1)) THEN         
            IF (LPLC(1)) WRITE (21,117)
         ELSE
            IF (LPLC(1)) WRITE (21,118)
         ENDIF
         WRITE (21,111)
         WRITE (21,113)
         WRITE (21,115)
         IF (LC(2)) THEN
            IF (LPLC(2)) WRITE (21,117)
         ELSE
            IF (LPLC(2)) WRITE (21,118)
         ENDIF
      ENDIF
C
      IF ((LC(1) .AND. TYPE(4,1) .AND. CSPEC(2)) .OR.
     &    (LC(2) .AND. TYPE(4,2) .AND. CSPEC(1))) THEN
C
C        Category 8
C
         WRITE (21,123) 8
         WRITE (21,100)
         WRITE (21,103) CUTOFF
         WRITE (21,110)
         WRITE (21,113)
         IF (CSPEC(1)) THEN
            WRITE (21,116) CEV(1)
         ELSE
            WRITE (21,115)
         ENDIF
         IF (LC(1)) THEN
            IF (LPLC(1)) WRITE (21,117)
         ELSE
            IF (LPLC(1)) WRITE (21,118)
         ENDIF
         WRITE (21,111)
         WRITE (21,113)
         IF (CSPEC(2)) THEN
            WRITE (21,116) CEV(2)
         ELSE
            WRITE (21,115)
         ENDIF
         IF (LC(2)) THEN
            IF (LPLC(2)) WRITE (21,117)
         ELSE
            IF (LPLC(2)) WRITE (21,118)
         ENDIF
      ENDIF
C
      IF ((((.NOT. LC(1)) .AND. TYPE(4,1) .AND. (.NOT. CSPEC(1)))
     &    .AND. CSPEC(2)) .OR. 
     &    (((.NOT. LC(2)) .AND. TYPE(4,2) .AND. (.NOT. CSPEC(2)))
     &    .AND. CSPEC(1))) THEN
C
C        Category 9
C
         WRITE (21,123) 9
         WRITE (21,101)
         WRITE (21,110)
         WRITE (21,113)
         IF (CSPEC(1)) THEN
            WRITE (21,116) CEV(1)
         ELSE
            WRITE (21,115)
         ENDIF
         IF (LC(1)) THEN
            IF (LPLC(1)) WRITE (21,117)
         ELSE
            IF (LPLC(1)) WRITE (21,118)
         ENDIF
         WRITE (21,111)
         WRITE (21,113)
         IF (CSPEC(2)) THEN
            WRITE (21,116) CEV(2)
         ELSE
            WRITE (21,115)
         ENDIF
         IF (LC(2)) THEN
            IF (LPLC(2)) WRITE (21,117)
         ELSE
            IF (LPLC(2)) WRITE (21,118)
         ENDIF
      ENDIF
C
      IF (CSPEC(1) .AND. CSPEC(2)) THEN
C
C        Category 10
C
         WRITE (21,123) 10
         WRITE (21,101)
         WRITE (21,103) CUTOFF
         IF (LASTEV .EQ. -5) THEN
            WRITE (21,120)
         ELSE
            IF (LASTEV .EQ. 0) THEN
               WRITE (21,107)
            ELSE
               IF (LASTEV .EQ. 1) THEN
                  WRITE (21,108)
               ELSE
                  WRITE (21,109) LASTEV
               ENDIF
            ENDIF
         ENDIF
         WRITE (21,110)
         WRITE (21,113)
         WRITE (21,116) CEV(1)
         IF (LC(1)) THEN
            IF (LPLC(1)) WRITE (21,117)
         ELSE
            IF (LPLC(1)) WRITE (21,118)
         ENDIF
         WRITE (21,111)
         WRITE (21,113)
         WRITE (21,116) CEV(2)
         IF (LC(2)) THEN
            IF (LPLC(2)) WRITE (21,117)
         ELSE
            IF (LPLC(2)) WRITE (21,118)
         ENDIF
      ENDIF
      WRITE (21,*)
      RETURN
  100 FORMAT(' This problem has simple spectrum.')
  101 FORMAT(' This problem may have non-simple spectrum.')
  102 FORMAT(' There is no continuous spectrum.')
  103 FORMAT(' There is continuous spectrum in [Ev, infinity) where'
     &      ,' Ev =',G15.6)
  104 FORMAT(' There is continuous spectrum consisting of the entire ',
     &       'real line.')
  105 FORMAT(' The set of eigenvalues is bounded below.')
  106 FORMAT(' There are infinitely many negative and infinitely many',
     &       ' positive',/,'    eigenvalues (unbounded in either',
     &       ' direction).')
  107 FORMAT(' There appear to be no eigenvalues below the start of',
     &       ' the continuous spectrum.')
  108 FORMAT(' There appears to be 1 eigenvalue below the start of the',
     &       ' continuous spectrum.')
  109 FORMAT(' There appear to be ',I12,' eigenvalues below the start'
     &       ,/,'    of the continuous spectrum.')
  110 FORMAT(' At endpoint A')
  111 FORMAT(' At endpoint B')
  112 FORMAT('    the problem is regular;')
  113 FORMAT('    the problem is singular;')
  114 FORMAT('    it is nonoscillatory for all Ev.')
  115 FORMAT('    it is oscillatory for all Ev.')
  116 FORMAT('    it is nonoscillatory for Ev <',G15.6,' and oscillatory
     & otherwise.')
  117 FORMAT('    It is limit circle.')
  118 FORMAT('    It is limit point.')
  119 FORMAT('    The constants for the Friedrichs boundary conditions',
     &       ' are',/,4E18.8)
  120 FORMAT(' There appear to be infinitely many eigenvalues below',
     &       ' the start',/,'    of the continuous spectrum.')
  121 FORMAT(' There are infinitely many eigenvalues, bounded below.') 
  122 FORMAT(' The nature of the spectrum is unknown; there is likely',
     &       ' to be ',/,' continuous spectrum.')
  123 FORMAT(' The spectral category is',I3,'.')
      END
C-----------------------------------------------------------------------
       SUBROUTINE EXTRAP(V,TOL,IROW,MAXCOL,FULL,TIGHT,MODE,VSAVE,
     &                   IPRINT,ERROR,DONE)
      INTEGER IROW,MAXCOL,MODE,IPRINT
      LOGICAL FULL,TIGHT,DONE
      DOUBLE PRECISION V,TOL,VSAVE(*),ERROR
C
C     Use Richardson's h**2 extrapolation (based on doubling) when 
C     suitable, otherwise use Wynn's acceleration scheme.
C
C     Input:
C        V       = real value at current level.
C        TOL     = real tolerance.
C        IROW    = integer giving current row index (1 .le. IROW).
C        MAXCOL  = integer giving maximum number of columns in table.
C        FULL    = logical, True iff entire table is to be computed.
C        TIGHT   = logical, True for conservative convergence tests.
C        MODE    = integer, value is 
C                  0   both Richardson and Wynn algorithms can be used;
C                  1   only Richardson is used;
C                  2   only Wynn is used.
C        IPRINT  = integer controlling amount of printing.
C     Output:
C        V       = real, best output estimate.
C        VSAVE   = real vector, holds previous level values.
C        DONE    = logical, True iff Error is sufficiently small.
C
C     If FULL is True, then the entire acceleration array is produced
C     (through row IROW); if False, then only the next row is appended.
C     Hence, the choice of FULL = True requires more work, but it may
C     save some global storage.  The vector VSAVE contains the V values
C     for levels 0 through max(IROW-1,MAXCOL).
C
      INTEGER I,IMIN,J,MAXJ,NCOL
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &                 DIFF,EPS,ETEMP,R(12,8),RATIO(12),RHIGH,RLOW,
     &                 RTOL,T,TOL1,TOL2,VTEMP,W(40,11),
     &                 ZERO,TENTH,HALF,ONE,TWO
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      SAVE R,W
C
C     The local arrays RATIO(*), R(*,*), and W(*,*) must be declared to 
C     have at least as many rows as the value of MAXLVL initialized in
C     routine START.
C
      PARAMETER (ZERO = 0.0, TENTH = 0.1D0, HALF = 0.5D0, ONE = 1.0,
     &           TWO = 2.0)
C
      EPS = 0.2
      DONE = .FALSE.
      ETEMP = ONE/U
      VTEMP = V
      VSAVE(IROW) = V
      TOL1 = TOL
      TOL2 = TOL/3.0
      IF (MODE .EQ. 2) THEN
         MAXJ = MAXCOL
         GOTO 40
      ELSE
         MAXJ = 11
      ENDIF
      RLOW = MAX(3.0,4.2-0.2*(IROW-1))
      RHIGH = MIN(5.0,3.8+0.2*(IROW-1))
      RTOL = U
C
C     Analyze the rate of convergence to determine NCOL and tolerances.
C
      NCOL = 2
      DO 10 I = 1,IROW
         R(I,1) = VSAVE(I)
         RATIO(I) = ONE/U
         IF (I .GE. 3) THEN
            T = R(I,1)-R(I-1,1)
            IF (T .NE. ZERO) THEN
               RATIO(I) = (R(I-1,1)-R(I-2,1))/T
            ELSE
               V = R(I,1)
               DONE = .TRUE.
               RETURN
            ENDIF
            IF (((RATIO(I) .GE. RLOW) .AND. (RATIO(I) .LE. RHIGH)) 
     &         .OR. (.NOT. TIGHT)) THEN
               RTOL = TOL1
               NCOL = MIN(MAXCOL,NCOL+1)
            ELSE
               RTOL = TOL2
               IF (RATIO(I) .LT. ZERO) RTOL = TENTH*TOL2
               IF (RATIO(I) .LT. TWO) NCOL = 2
            ENDIF
         ENDIF
   10 CONTINUE
      IF (FULL) THEN
         IMIN = 2
      ELSE
         IMIN = IROW
      ENDIF
C
C     Use Richardson's h^2 extrapolation.  The number of columns used
C     is a function of the amount of data (IROW), the requested order
C     (MAXCOL), the observed rate of convergence (NCOL), and the amount
C     of storage allocated to R(*,*).
C
      DO 30 I = IMIN,IROW
         DO 20 J = 2,MIN(I,NCOL,8)
            DIFF = (R(I,J-1)-R(I-1,J-1))/(4**(J-1)-1)
            R(I,J) = R(I,J-1)+DIFF
            IF ((.NOT. FULL) . OR. (I .EQ. IROW)) THEN
              T = ABS(DIFF)
               IF (T .LE. ETEMP) THEN
                  ETEMP = T
                  VTEMP = R(I,J)
                  IF (T .LE. RTOL) THEN
                     DONE = .TRUE.
                     V = VTEMP
                     ERROR = ETEMP
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
   20    CONTINUE
   30 CONTINUE
      V = VTEMP
      ERROR = ETEMP
      IF (IROW .LT. 4) RETURN
C
C     Test for rate of convergence other than second order.
C
      IF (ABS(RATIO(IROW)-RATIO(IROW-1))+ABS(RATIO(IROW-1)-
     &        RATIO(IROW-2)) .GT. EPS) RETURN
      IF ((3.5 .LT. RATIO(IROW)) .AND. (RATIO(IROW) .LT. 4.5)) RETURN
      IF (RATIO(IROW) .LT. ONE) RETURN
      IF (MODE .NE. 0) RETURN
C
C     Use Wynn's algorithm.
C
  40  IF (IPRINT .GE. 4) THEN
         DIFF = LOG(ABS(RATIO(IROW)))/LOG(TWO)
         WRITE (21,50) DIFF
  50     FORMAT(' In EXTRAP: using Wynn`s acceleration; rate = ',F8.5)
      ENDIF
      W(1,1) = VSAVE(1)
      ETEMP = ONE/U
      DO 60 I = 2,IROW
         W(I,1) = VSAVE(I)
         DIFF = W(I,1)-W(I-1,1)
         IF ((I .EQ. IROW) .OR. (MODE .EQ. 2)) THEN
            T = ABS(DIFF)
            IF (T .LE. ETEMP) THEN
               ETEMP = T
               VTEMP = W(I,1)
               IF (T .LE. TOL2) THEN
                  V = VTEMP
                  ERROR = ETEMP
                  DONE = .TRUE.
                  RETURN
               ENDIF
            ENDIF
         ENDIF    
         IF (DIFF .NE. ZERO) THEN
            W(I,2) = ONE/DIFF
         ELSE
            V = W(I,1)
           ERROR = ZERO
            DONE = .TRUE.
            RETURN
         ENDIF
   60 CONTINUE
      DO 80 J = 3,MIN(IROW,MAXJ)
         DO 70 I = J,IROW
            DIFF = W(I,J-1)-W(I-1,J-1)
            IF (DIFF .NE. ZERO) THEN
               DIFF = ONE/DIFF
            ELSE
               IF (MOD(J,2) .EQ. 0) THEN
                  V = W(I,J-1)
                  ERROR = ZERO
                  DONE = .TRUE.
                  RETURN
               ELSE
                  DIFF = ONE/U**2
               ENDIF
            ENDIF
            W(I,J) = W(I-1,J-2)+DIFF
            IF ((MOD(J,2) .EQ. 1) .AND. ((I .EQ. IROW) .OR.
     &          (MODE .EQ. 2))) THEN
               T = ABS(DIFF)
               IF (T .LE. ETEMP) THEN
                  ETEMP = T
                  VTEMP = W(I,J)
                  IF (T .LE. TOL2) THEN
                     V = VTEMP
                     ERROR = ETEMP
                     DONE = .TRUE.
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
   70    CONTINUE
   80 CONTINUE
      V = VTEMP
      ERROR = ETEMP
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GETEF(EV,EFNORM,IPRINT,X,EFIN)
      INTEGER IPRINT
      DOUBLE PRECISION EV,EFNORM,X(*)
      LOGICAL EFIN(2,2)
C
C     Compute an eigenfunction for one fixed mesh.
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &        KCLASS(2),I,J,JDU,JLAST,JS,JU,JX,KLVL,MIDDLE(2),MODE
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2),
     &        ALLOK,SYMM
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &        CP(2),CR(2),CUTOFF,D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     &        ETA(2,2),PNU(2),
     &        CHI,DPSI,DV,DW,FNORM,FSCALE,FSUM,H,HALFH,HOM,OM,OSCALE,
     &        PDV,PDW,PN,PROD,PSI,RATIO,RN,RNORM,RSCALE,RSUM,SCALE,
     &        T,TAU,TAUHH,TAUMAX,V,VNEW,W,WNEW,XLEFT,XRIGHT,Z,
     &        ZERO,C10M4,HALF,ONE,TWO,THREE,FIVE,C15,C21
      COMMON /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, C10M4 = 1.D-4, HALF = 0.5D0, ONE = 1.0,
     &           TWO = 2.0, THREE = 3.0, FIVE = 5.0, C15 = 15.0,
     &           C21 = 21.0)
C
      KLVL = 2**LEVEL
C
C     For this EV and mesh, calculate FSCALE, RSCALE, and MIDDLE(*).
C         FSCALE = sum of logs of scale factors 1 through m
C         RSCALE = sum of logs of scale factors m+1 through N
C         MIDDLE(1), MIDDLE(2) describe the coordinates of the matching
C            point M for the shooting; in particular,
C            M = X(MIDDLE(1)-1) + MIDDLE(2)*H(MIDDLE(1)-1)/2**LEVEL
C     The matching point is chosen to be roughly (a+b)/2 if either
C     a = -b and Tau(x) > 0 near 0, or if Tau(x) > 0 for all x;
C     otherwise, it is chosen to roughly maximize Tau(x).
C
      FSCALE = ZERO
      RSCALE = ZERO
      TAUMAX = -ONE/U
      ALLOK = .TRUE.
      SYMM = .FALSE.
      IF (A .EQ. -B) SYMM = .TRUE.
      EFIN(1,1) = .TRUE.
      EFIN(2,1) = .TRUE.
      EFIN(1,2) = .TRUE.
      EFIN(2,2) = .TRUE.
      JX = (NXINIT+1)/2
      DO 20 I = 2,NXINIT
         MODE = 0
         XLEFT = X(I-1)
         H = X(I)-XLEFT
         IF (I .EQ. 2) THEN
            IF (KCLASS(1) .GE. 9) MODE = 1
            IF (.NOT. AFIN) THEN
               MODE = 3
               XLEFT = ZERO
               H = -ONE/X(2)
            ENDIF
         ENDIF
         IF (I .EQ. NXINIT) THEN
            IF (KCLASS(2) .GE. 9) MODE = 2
            IF (.NOT. BFIN) THEN
               MODE = 4
               H = ONE/X(I-1)
               XLEFT = -H
            ENDIF
         ENDIF
         H = H/KLVL
         HALFH = HALF*H
         DO 10 J = 1,KLVL
            Z = XLEFT+HALFH
            XLEFT = XLEFT+H
            CALL STEP(Z,H,EV,PN,RN,TAU,OM,HOM,PSI,DPSI,SCALE,MODE)
            IF (TAU .LT. ZERO) ALLOK = .FALSE.
            IF (TAU .GT. TAUMAX) THEN
               TAUMAX = TAU
               MIDDLE(1) = I
               MIDDLE(2) = J
               FSCALE = FSCALE+RSCALE+SCALE
               OSCALE = SCALE
               RSCALE = ZERO
            ELSE
               RSCALE = RSCALE+SCALE
            ENDIF
   10    CONTINUE
         IF ((I .EQ. JX) .AND. (TAU .LT. ZERO)) SYMM = .FALSE.
   20 CONTINUE
      IF (ALLOK .OR. SYMM) THEN
         MIDDLE(1) = JX
         MIDDLE(2) = MAX(KLVL-1,1)
      ENDIF
      IF ((.NOT. AFIN) .OR. (.NOT. BFIN)) THEN
C
C        Don't split near infinity!
C
         IF ((.NOT. AFIN) .AND. (MIDDLE(1) .EQ. 2)) THEN
            IF (REG(2)) THEN
               MIDDLE(1) = NXINIT
            ELSE
               MIDDLE(1) = JX
            ENDIF
            MIDDLE(2) = MAX(KLVL-1,1)
         ENDIF
         IF ((.NOT. BFIN) .AND. (MIDDLE(1) .EQ. NXINIT)) THEN
            IF (REG(1)) THEN
               MIDDLE(1) = 2
            ELSE
               MIDDLE(1) = JX
            ENDIF
            MIDDLE(2) = 1
         ENDIF
      ENDIF
      IF ((LEVEL .GT. 1) .AND. (MIDDLE(2) .EQ. KLVL)) THEN
         MIDDLE(2) = KLVL-1
         FSCALE = FSCALE-OSCALE
         RSCALE = RSCALE+OSCALE
      ENDIF
      IF (IPRINT .GE. 4) WRITE (21,21) MIDDLE(1),MIDDLE(2)
   21 FORMAT(' Coordinates of matching point =',2I6)
      JU = NXINIT
      JDU = 2*NXINIT
      JS = 3*NXINIT
C
C     Shoot from x=A to the middle.
C
      V = A2-A2P*EV
      PDV = A1-A1P*EV
      FNORM = ZERO
      IF (KCLASS(1) .LT. 9) THEN
         SCALE = MAX(ABS(V),ABS(PDV))
         V = V/SCALE
         PDV = PDV/SCALE
         MODE = 0
      ELSE
         V = ONE
         PDV = D(4,1)
         MODE = 1
      ENDIF
      IF (.NOT. AFIN) MODE = 3
      FSUM = -FSCALE
      IF ((IPRINT .GE. 5) .AND. (MODE .EQ. 0)) WRITE (21,35)
     &                                         X(1),V,PDV,FSUM
      DO 40 I = 2,MIDDLE(1)
         X(JU+I-1) = V
         X(JDU+I-1) = PDV
         X(JS+I-1) = FSUM
         IF (MODE .EQ. 0) THEN
            XLEFT = X(I-1)
            H = X(I)-XLEFT
         ELSE
            XLEFT = ZERO
            IF (MODE .EQ. 1) THEN
               H = ((X(2)-X(1))/D(1,1))**ETA(1,1)
            ELSE
               H = -ONE/X(2)
            ENDIF
         ENDIF
         H = H/KLVL
         HALFH = HALF*H
         JLAST = KLVL
         IF (I .EQ. MIDDLE(1)) JLAST = MIDDLE(2)
         DO 30 J = 1,JLAST
            Z = XLEFT+HALFH
            XLEFT = XLEFT+H
            CALL STEP(Z,H,EV,PN,RN,TAU,OM,HOM,PSI,DPSI,SCALE,MODE)
            DV = PDV/PN
            IF (ABS(TAU)*H*H .GE. C10M4) THEN
               IF (TWO*(FSUM+SCALE) .GT. -UNDER) THEN
                  FNORM = FNORM+RN*(PSI*DPSI*(V*V-DV*DV/TAU)/TWO+
     &                              V*PSI*DV*PSI)*EXP(TWO*(FSUM+SCALE))
                  IF (TWO*FSUM .GT. -UNDER) FNORM = FNORM+RN*
     &                              EXP(TWO*FSUM)*H*(V*V+DV*DV/TAU)/TWO
               ENDIF
            ELSE
               IF (TWO*FSUM .GT. -UNDER) THEN
                  TAUHH = TAU*H*H
                  FNORM = FNORM+RN*H*EXP(TWO*FSUM)*(
     &                     V*V*(ONE+TAUHH*(TAUHH/FIVE-ONE)/THREE)
     &                    +H*V*DV*(ONE+TAUHH*(TWO*TAUHH/C15-ONE)/THREE)
     &                    +(H*DV)**2*(ONE+TAUHH*(TWO*TAUHH/C21-ONE)/
     &                                FIVE)/THREE )
               ENDIF
            ENDIF
            FSUM = FSUM+SCALE
            VNEW = DPSI*V+PSI*DV
            PDV = -PN*TAU*PSI*V+DPSI*PDV
            V = VNEW
   30    CONTINUE
         IF (MODE .EQ. 1) THEN
C
C           Convert from V(t) to u(x).
C
            IF (ETA(2,1) .LT. ZERO) THEN
               X(JU+1) = ONE/U
               EFIN(1,1) = .FALSE.
            ELSE
               IF (ETA(2,1) .EQ. ZERO) THEN
                  X(JU+1) = D(2,1)
               ELSE
                  X(JU+1) = ZERO
               ENDIF
            ENDIF
            Z = ETA(1,1)*ETA(2,1)+EP(1)-ONE
            IF (Z .LT. ZERO) THEN
               X(JDU+1) = SIGN(ONE/U,ETA(2,1))
               EFIN(2,1) = .FALSE.
            ELSE
               IF (Z .GT. ZERO) THEN
                  X(JDU+1) = ZERO
               ELSE
                  X(JDU+1) = D(2,1)*CP(1)*ETA(1,1)*ETA(2,1)/
     &                              D(1,1)**(ETA(1,1)*ETA(2,1))
               ENDIF
            ENDIF
            H = X(2)-A
            PN = CP(1)*H**(EP(1)-ONE)
            T = (H/D(1,1))**ETA(1,1)
            CHI = D(2,1)*T**ETA(2,1)
            V = V*CHI
            PDV = PN*ETA(1,1)*(ETA(2,1)*V+CHI*PDV*
     &                         T**(ONE-TWO*EMU(1)))
            IF (IPRINT .GE. 5) WRITE (21,35) X(1),X(JU+1),X(JDU+1)
         ENDIF
         MODE = 0
         IF (IPRINT .GE. 5) WRITE (21,35) XLEFT,V,PDV,FSUM
   35    FORMAT(G16.6,3D15.6)
   40 CONTINUE
C
C     Shoot from x=B to the middle.
C
      RNORM = ZERO
      MODE =  0
      IF (KCLASS(2) .LT. 9) THEN
         SCALE = MAX(ABS(B1),ABS(B2))
         W = -B2/SCALE
         PDW = B1/SCALE
         MODE = 0
      ELSE
         W = -ONE
         PDW = -D(4,2) 
         MODE = 2
      ENDIF
      IF (.NOT. BFIN) MODE = 4
      RSUM = -RSCALE
      IF ((IPRINT .GE. 5) .AND. (MODE .EQ. 0)) WRITE (21,35)
     &                                         X(NXINIT),W,PDW,RSUM
      DO 60 I = NXINIT,MIDDLE(1),-1
         X(JU+I) = W
         X(JDU+I) = PDW
         X(JS+I) = RSUM
         IF (MODE .EQ. 0) THEN
            XRIGHT = X(I)
            H = XRIGHT-X(I-1)
         ELSE
            XRIGHT = ZERO
            IF (MODE .EQ. 2) THEN
               H = ((X(NXINIT)-X(NXINIT-1))/D(1,2))**ETA(1,2)
            ELSE
               H = ONE/X(I-1)
            ENDIF
         ENDIF
         H = H/KLVL
         HALFH = HALF*H
         JLAST = KLVL
         IF (I .EQ. MIDDLE(1)) JLAST = JLAST-MIDDLE(2)
         DO 50 J = 1,JLAST
            Z = XRIGHT-HALFH
            XRIGHT = XRIGHT-H
            CALL STEP(Z,H,EV,PN,RN,TAU,OM,HOM,PSI,DPSI,SCALE,MODE)
            DW = PDW/PN
            IF (ABS(TAU)*H*H .GE. C10M4) THEN
               IF (TWO*(RSUM+SCALE) .GT. -UNDER) THEN
                  RNORM = RNORM+RN*(PSI*DPSI*(W*W-DW*DW/TAU)/TWO-
     &                              W*PSI*DW*PSI)*EXP(TWO*(RSUM+SCALE))
                  IF (TWO*RSUM .GT. -UNDER) RNORM = RNORM+RN*
     &                              EXP(TWO*RSUM)*H*(W*W+DW*DW/TAU)/TWO
               ENDIF
            ELSE
               IF (TWO*RSUM .GT. -UNDER) THEN
                  TAUHH = TAU*H*H
                  RNORM = RNORM+RN*H*EXP(TWO*RSUM)*(
     &                     W*W*(ONE+TAUHH*(TAUHH/FIVE-ONE)/THREE)
     &                    -H*W*DW*(ONE+TAUHH*(TWO*TAUHH/C15-ONE)/THREE)
     &                    +(H*DW)**2*(ONE+TAUHH*(TWO*TAUHH/C21-ONE)/
     &                     FIVE)/THREE )
               ENDIF
            ENDIF
            RSUM = RSUM+SCALE
            WNEW = DPSI*W-PSI*DW
            PDW = PN*TAU*PSI*W+DPSI*PDW
            W = WNEW
   50    CONTINUE
         IF (MODE .EQ. 2) THEN
C
C           Convert from V(t) to u(x).
C
            IF (ETA(2,2) .LT. ZERO) THEN
               X(JU+NXINIT) = -ONE/U
               EFIN(1,2) = .FALSE.
            ELSE
               IF (ETA(2,2) .EQ. ZERO) THEN
                  X(JU+NXINIT) = -D(2,2)
               ELSE
                  X(JU+NXINIT) = ZERO
               ENDIF
            ENDIF
            Z = ETA(1,2)*ETA(2,2)+EP(2)-ONE
            IF (Z .LT. ZERO) THEN
               X(JDU+NXINIT) = SIGN(ONE/U,ETA(2,2))
               EFIN(2,2) = .FALSE.
            ELSE
               IF (Z .GT. ZERO) THEN
                  X(JDU+NXINIT) = ZERO
               ELSE
                  X(JDU+NXINIT) = D(2,2)*CP(2)*ETA(1,2)*ETA(2,2)/
     &                                   D(1,2)**(ETA(1,2)*ETA(2,2))
               ENDIF
            ENDIF
            IF (IPRINT .GE. 5) WRITE(21,35) X(NXINIT),X(JU+NXINIT),
     &                                       X(JDU+NXINIT)
            H = X(NXINIT)-X(NXINIT-1)
            PN = CP(2)*H**(EP(2)-ONE)
            T = (H/D(1,2))**ETA(1,2)
            CHI = D(2,2)*T**ETA(2,2)
            W = CHI*W
            PDW = PN*ETA(1,2)*(CHI*PDW*T**(ONE-TWO*EMU(2))-ETA(2,2)*W)
         ENDIF
         IF (MODE .EQ. 4) THEN
            X(JU+NXINIT) = ZERO
            X(JDU+NXINIT) = ONE
            IF (IPRINT .GE. 5) WRITE(21,35) X(NXINIT),X(JU+NXINIT),
     &                                      X(JDU+NXINIT)
         ENDIF
         MODE = 0
         IF ((JLAST .NE. 0) .AND. (IPRINT .GE. 5)) 
     &                            WRITE (21,35) XRIGHT,W,PDW,RSUM
   60 CONTINUE
      IF (ABS(W) .GE. ABS(PDW)) THEN
         RATIO = V/W
         IF (IPRINT .GE. 5) WRITE (21,61) RATIO*PDW-PDV,RATIO
   61    FORMAT('  DuHat jump, ratio =',2D24.15)
      ELSE
         RATIO = PDV/PDW
         IF (V*W*RATIO .LT. ZERO) RATIO = -RATIO
         IF (IPRINT .GE. 5) WRITE (21,62) RATIO*W-V,RATIO
   62    FORMAT('  UHat jump, ratio =',2D24.15)
      ENDIF
C
C     Calculate weighted 2-norm and scale approximate eigenfunction.
C
      FSCALE = EXP(-FSUM)
      RSCALE = EXP(-RSUM)
      EFNORM = SQRT(FNORM*FSCALE**2+RNORM*(RATIO*RSCALE)**2)
      SCALE = LOG(EFNORM)
      IF (IPRINT .GE. 5) WRITE (21,65) EFNORM
   65 FORMAT('  EFnorm =',D24.15)
      DO 70 I = 1,NXINIT
         TAU = X(JS+I)-SCALE
         IF (TAU .LE. -UNDER) THEN
            X(JU+I) = ZERO
            X(JDU+I) = ZERO
         ELSE
            PROD = EXP(TAU)
            IF (I .GE. MIDDLE(1)) PROD = PROD*RATIO
            X(JU+I) = X(JU+I)*PROD
            X(JDU+I) = X(JDU+I)*PROD
         ENDIF
   70 CONTINUE
      IF (IPRINT .GE. 4) THEN
         WRITE (21,75)
   75    FORMAT(10X,'x',15X,'Uhat(x)',13X,'PUhat`(x)')
         DO 85 I = 1,NXINIT
            WRITE (21,80) X(I),X(JU+I),X(JDU+I)
   80       FORMAT(G16.6,2D20.8)
   85    CONTINUE
      ENDIF
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GETRN(EV,ALPHA,CSPEC,DENOM,RSUBN,X)
      DOUBLE PRECISION EV,ALPHA,DENOM,RSUBN,X(*)
      LOGICAL CSPEC(*)
C
C     Compute the RsubN value from the weighted eigenfunction 2-norm
C     when standard shooting is stable and an accurate eigenvalue is
C     available (from the asymptotic formulas).
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &        KCLASS(2),   I,J,KLVL
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2)
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,URN,UNDER,
     &        CP(2),CR(2),CUTOFF,D(4,2),EMU(2),EP(2),EQLNF(2),
     &        ER(2),ETA(2,2),PNU(2),
     &        DPHI,DPSI,DU,DUSAVE,FNORM,H,HALFH,HOM,HSAVE,OM,PDU,PHI,
     &        PN,PSI,RN,SCALE,TAU,TAUHH,U,UNEW,USAVE,XLEFT,Z,
     &        ZERO,C10M4,HALF,ONE,TWO,THREE,FIVE,C15,C21
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,URN,UNDER
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      PARAMETER (ZERO = 0.0, C10M4 = 1.D-4, HALF = 0.5D0, ONE = 1.D0,
     &           TWO = 2.0, THREE = 3.0, FIVE = 5.0, C15 = 15.0,
     &           C21 = 21.0)
C
      U = A2-A2P*EV
      PDU = A1-A1P*EV
      FNORM = ZERO
      KLVL = 2**LEVEL
C
C     Shoot from x=A to x=B.
C
      DO 20 I = 2,NXINIT
         XLEFT = X(I-1)
         H = (X(I)-XLEFT)/KLVL
         HALFH = HALF*H
         DO 10 J = 1,KLVL
            Z = XLEFT+HALFH
            XLEFT = XLEFT+H
            CALL STEP(Z,H,EV,PN,RN,TAU,OM,HOM,PSI,DPSI,SCALE,0)
            SCALE = EXP(SCALE)
            PHI = PSI*SCALE
            DPHI = DPSI*SCALE
            DU = PDU/PN
            IF (ABS(TAU)*H*H .GE. C10M4) THEN
               FNORM = FNORM+RN*(PHI*DPHI*(U*U-DU*DU/TAU)/TWO
     &                 +U*PHI*DU*PHI+H*(U*U+DU*DU/TAU)/TWO)
            ELSE
               TAUHH = TAU*H*H
               FNORM = FNORM+RN*H*(U*U*(ONE+TAUHH*(TAUHH/FIVE-ONE)/
     &                                  THREE)
     &                 +H*U*DU*(ONE+TAUHH*(TWO*TAUHH/C15-ONE)/THREE)
     &                 +(H*DU)**2*(ONE+TAUHH*(TWO*TAUHH/C21-ONE)/
     &                             FIVE)/THREE )
            ENDIF
            IF ((I .EQ. NXINIT) .AND. (J .EQ. KLVL) .AND. CSPEC(1)) THEN
               HSAVE = H
               USAVE = ABS(U)
               DUSAVE = ABS(PDU)
            ENDIF
            UNEW = DPHI*U+PHI*DU
            PDU = -PN*TAU*PHI*U+DPHI*PDU
            U = UNEW
            IF ((I .EQ. 2) .AND. (J .EQ. 1) .AND. CSPEC(2)) THEN
               HSAVE = H
               USAVE = ABS(U)
               DUSAVE = ABS(PDU)
            ENDIF
   10    CONTINUE
   20 CONTINUE
      IF (CSPEC(2)) THEN
         IF (REG(1) .OR. (PNU(1) .EQ. ZERO) .OR.
     &       (PNU(1) .EQ. ONE-EP(1))) THEN
             PHI = DENOM
         ELSE
            IF (USAVE .GE. DUSAVE) THEN
              PHI = USAVE/HSAVE**PNU(1)
            ELSE
              PHI = DUSAVE/(CP(1)*ABS(PNU(1))*HSAVE**(EP(1)+PNU(1)-ONE))
            ENDIF
         ENDIF
      ELSE
         IF (REG(2) .OR. (PNU(2) .EQ. ZERO) .OR.
     &       (PNU(2) .EQ. ONE-EP(2))) THEN
             PHI = DENOM
         ELSE
            IF (USAVE .GE. DUSAVE) THEN
              PHI = USAVE/HSAVE**PNU(2)
            ELSE
              PHI = DUSAVE/(CP(2)*ABS(PNU(2))*HSAVE**(EP(2)+PNU(2)-ONE))
            ENDIF
         ENDIF
      ENDIF
      RSUBN = ONE/(ALPHA+FNORM/PHI**2)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MESH(JOB,NEV,X,G,H,QLNF,Z,TOL,HMIN)
      LOGICAL JOB
      DOUBLE PRECISION X(*),G(*),H(*),QLNF(*),Z(*),TOL,HMIN
C
C     If JOB = True then calculate the initial mesh; redistribute so 
C     that H(*) is approximately equidistributed.  If JOB = False
C     then use the mesh input by the user.
C     
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &        KCLASS(2),I,ITS,J,JTOL,K,MAXITS,N,NADD
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &        CP(2),CR(2),CUTOFF,D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     &        ETA(2,2),PNU(2),
     &        DX,ENDA,ENDB,EPS,EQMAX,EQMIN,EV,GAMMA,P1,P2,P3,QMAX,QMIN,
     &        Q1,Q2,Q3,R1,R2,R3,WEIGHT,Y,Y1,Y2,Y3,
     &        ZERO,TENTH,HALF,ONE,TWO,FOUR
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2),
     &        DONE
      COMMON /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, TENTH = 0.1D0, HALF = 0.5D0, ONE = 1.0,
     &           TWO = 2.0, FOUR = 4.0)
      SAVE ENDA,ENDB
      DATA EPS/0.0001/
C
      IF (.NOT. JOB) THEN
         HMIN = B-A
         DO 5 I = 2,NXINIT
            HMIN = MIN(HMIN,X(I)-X(I-1))
    5    CONTINUE
         RETURN
      ENDIF
      N = NXINIT-1
      EV = NEV
C
C     Find an appropriate initial mesh.
C
      IF (AFIN) THEN
         X(1) = A
         IF (BFIN) THEN
            X(NXINIT) = B
            DX = (B-A)/N
            DO 10 I = 2,N
               X(I) = X(1)+(I-1)*DX
   10       CONTINUE
         ELSE
            IF (NEV .LT. 0) THEN
               ENDB = X(NXINIT)
            ELSE
               X(NXINIT) = (1+4*NEV)*ENDB
            ENDIF
            Y1 = X(1)/(ONE+ABS(X(1)))
            DX = (X(NXINIT)/(ONE+X(NXINIT))-Y1)/N
            DO 11 I = 2,N
               Y = Y1+(I-1)*DX
               X(I) = Y/(ONE-ABS(Y))
   11       CONTINUE
         ENDIF
      ELSE
         IF (NEV .LT. 0) THEN
            ENDA = X(1)
         ELSE
            X(1) = (1+4*NEV)*ENDA
         ENDIF
         IF (BFIN) THEN
            X(NXINIT) = B
            Y1 = X(NXINIT)/(ONE+ABS(X(NXINIT)))
            DX = (Y1-X(1)/(ONE-X(1)))/N
            DO 12 I = 2,N
               Y = Y1-(I-1)*DX
               X(NXINIT+1-I) = Y/(ONE-ABS(Y))
   12       CONTINUE
         ELSE
            Y1 = X(1)/(ONE-X(1))
            IF (NEV .LT. 0) THEN
               ENDB = X(NXINIT)
            ELSE
               X(NXINIT) = (1+4*NEV)*ENDB
            ENDIF
            Y2 = X(NXINIT)/(ONE+X(NXINIT))
            DX = (Y2-Y1)/N
            DO 13 I = 2,N
               Y = Y1+(I-1)*DX
               X(I) = Y/(ONE-ABS(Y))
   13       CONTINUE
            IF (ABS(X((NXINIT+1)/2)) .LT. TOL) X((NXINIT+1)/2) = ZERO
         ENDIF
      ENDIF
      JTOL = -LOG10(TOL)+HALF
      IF (REG(1) .AND. REG(2)) THEN
         MAXITS = 6
      ELSE
         MAXITS = 3
      ENDIF
C
C     Calculate H(*) and G(*).
C
      ITS = 1
      IF (.NOT. AFIN) THEN
         IF (X(2) .GE. ZERO) X(2) = MAX(HALF*X(1),-ONE)
      ENDIF
      IF (.NOT. BFIN) THEN
         IF (X(N) .LE. ZERO) X(N) = MIN(HALF*X(NXINIT),ONE)
      ENDIF
      QMIN = 1.E31
      QMAX = -QMIN
   20 GAMMA = ZERO
C
C    Equidistribute { [Qmax - Q]^2 * max[abs(p') , abs(q') , abs(r')] }.
C
      EQMAX = ZERO
      EQMIN = 1.E31
      DO 25 J = 1,N
         DX = X(J+1)-X(J)
         Y2 = X(J)+HALF*DX
         CALL COEFF(Y2,P2,Q2,R2)
         IF ((P2 .EQ. ZERO) .OR. (R2 .EQ. ZERO) .OR. (P2*R2 .LT. ZERO))
     &   THEN
            FLAG = -15
            RETURN
         ENDIF
         Y1 = MAX(Y2-EPS,X(1)+TWO*U*ABS(X(1)))
         CALL COEFF(Y1,P1,Q1,R1)
         Y3 = MIN(Y2+EPS,X(NXINIT)-TWO*U*ABS(X(NXINIT)))
         CALL COEFF(Y3,P3,Q3,R3)
         IF (LNF) THEN
            H(J) = ABS(Q3-Q1)/(Y3-Y1)
            QLNF(J) = Q2/R2
         ELSE
            H(J) = MAX(ABS(P3-P1),ABS(Q3-Q1),ABS(R3-R1))/(Y3-Y1)
            Y1 = SQRT(SQRT(R1*P1))
            Y2 = SQRT(SQRT(R2*P2))
            Y3 = SQRT(SQRT(R3*P3))
            Y = SQRT(P2/R2)
            QLNF(J) = Q2/R2 + Y*((Y3-Y1)*(SQRT(P3/R3)-SQRT(P1/R1))/FOUR
     &                          +(Y3-TWO*Y2+Y1)*Y)/(Y2*EPS**2)
            IF (ABS(QLNF(J)) .LE. EPS) QLNF(J) = ZERO
         ENDIF
         QMAX = MAX(QMAX,QLNF(J))
         QMIN = MIN(QMIN,QLNF(J))
   25 CONTINUE
      Y = MAX(QMAX-QMIN,ONE)
      EV = 100.0+MAX(ZERO,QMIN)
      DO 30 J = 1,N
         DX = X(J+1)-X(J)
         IF (QLNF(J) .LE. EV) THEN
            WEIGHT = 3.0*((QMAX-QLNF(J))/Y)**2+ONE
            H(J) = MAX(H(J)*WEIGHT,U)
         ELSE
            Y2 = TWO*DX*SQRT(QLNF(J)-EV)
            IF (Y2 .LE. UNDER) THEN
               WEIGHT = EXP(-Y2)
               H(J) = MAX(WEIGHT*H(J),U)
            ELSE
               H(J) = U
            ENDIF
         ENDIF
         IF ((.NOT. AFIN) .AND. (X(J+1) .LT. ZERO)) THEN
            H(J) = MAX(H(J)*EXP(X(J+1)),U)
            IF ((J .EQ. 1) .AND. (H(1) .EQ. U)) X(1) = HALF*(X(1)+X(2))
         ENDIF
         IF ((.NOT. BFIN) .AND. (X(J) .GT. ZERO)) THEN
            H(J) = MAX(H(J)*EXP(-X(J)),U)
            IF ((J .EQ. N) .AND. (H(N) .EQ. U)) X(NXINIT) =
     &                                          HALF*(X(N)+X(NXINIT))
         ENDIF
         EQMIN = MIN(EQMIN,H(J)*DX)
         EQMAX = MAX(EQMAX,H(J)*DX)
   30 CONTINUE
      NCOEFF = NCOEFF+3*N
      IF (EQMAX-EQMIN .LE. MAX(TENTH*EQMAX,U/TENTH)) GOTO 75
C
C     Use a roughly locally quasi-uniform mesh.
C
      GAMMA = ZERO
      DO 35 I = 1,N
         GAMMA = GAMMA+H(I)
   35 CONTINUE
      GAMMA = ONE/GAMMA
      DO 45 I = 1,N
         Y = ZERO
         DO 40 J = 1,N
            Y = MAX( Y,H(J)/(ONE+GAMMA*ABS(X(I)+X(I+1)-X(J)-X(J+1))
     &                *H(J)) )
   40    CONTINUE
         Z(I) = Y
   45 CONTINUE
      DO 50 I = 1,N
         H(I) = Z(I)
   50 CONTINUE
      G(1) = ZERO
      DO 55 J = 1,N
         G(J+1) = G(J)+H(J)*(X(J+1)-X(J))
   55 CONTINUE
      GAMMA = G(N+1)/N
C
C     Redistribution algorithm:
C
      Y = GAMMA
      I = 1
      DO 65 J = 1,N
   60    IF (Y .LE. G(J+1)) THEN
            I = I+1
            Z(I) = X(J)+(Y-G(J))/H(J)
            Y = Y+GAMMA
            GOTO 60
         ENDIF
   65 CONTINUE
      Z(1) = X(1)
      Z(NXINIT) = X(NXINIT)
      DONE = .TRUE.
      DO 70 J = 2,N
         IF (ABS(Z(J)-X(J)) .GT. TENTH*(Z(J+1)-Z(J-1))) DONE = .FALSE.
         X(J) = Z(J)
   70 CONTINUE
      IF (.NOT. DONE) THEN
         IF (ITS .LT. MAXITS) THEN
            ITS = ITS+1
            GOTO 20
         ENDIF
      ENDIF
   75 IF (.NOT. AFIN) THEN
         IF (X(2) .GE. ZERO) X(2) = -ONE
      ENDIF
      IF (.NOT. BFIN) THEN
         IF (X(N) .LE. ZERO) X(N) = ONE
      ENDIF
      DO 95 K = 1,2
         NADD = 0
         IF (KCLASS(K) .GT. 0) THEN
C
C           Add Nadd extra points near endpoint K.
C
            IF (KCLASS(K) .EQ. 1) THEN
               NADD = MIN(MAX((JTOL+2)/3,1),4)
               Y = TENTH
            ENDIF
            IF (KCLASS(K) .EQ. 2) THEN
               NADD = MIN(MAX(JTOL,3),4)
               Y = MIN(MAX(TENTH**(JTOL/3),1.D-3),1.D-2)
            ENDIF
            IF (KCLASS(K) .EQ. 3) NADD = MIN(MAX(JTOL/2,2),4)
            IF (KCLASS(K) .EQ. 4) THEN
               NADD = MIN(MAX((5+JTOL)/3,2),5)
               Y = TENTH**(NADD-1)
            ENDIF
            IF (KCLASS(K) .EQ. 5) THEN
               NADD = (2**JTOL)**0.4
               NADD = MIN(MAX(NADD,3),6)
               Y = MIN(MAX((0.1**JTOL)**(0.333),0.001),0.1)
            ENDIF
            IF (KCLASS(K) .EQ. 6) THEN
               NADD = MIN(MAX(JTOL,2),8)
               Y = 0.005
            ENDIF
            IF ((KCLASS(K) .EQ. 7) .OR. (KCLASS(K) .EQ. 10)) THEN
               NADD = MIN(MAX((2*JTOL+6)/3,3),8)
               Y = MIN(MAX(SQRT(TENTH*TOL),1.D-6),1.D-2)
            ENDIF
            IF (KCLASS(K) .EQ. 8) THEN
               NADD = (2**JTOL)**0.4
               NADD = MIN(MAX(NADD,2),6)
               Y = MIN(MAX((0.1**JTOL)**(0.333),0.001),0.05)
            ENDIF
            IF (KCLASS(K) .EQ. 9) THEN
               IF(LFLAG(5)) THEN
                  NADD = MIN(MAX((JTOL+4)/3,2),5)
                  Y = TENTH**(NADD-1)
               ELSE
                  NADD = MIN(MAX(2+JTOL*(JTOL-3)/40,2),4)
                  Y = 0.25
               ENDIF
            ENDIF
            IF (K .EQ. 1) THEN
               DO 80 I = NXINIT,2,-1
                  X(I+NADD) = X(I)
   80          CONTINUE
               NXINIT = NXINIT+NADD
               IF (AFIN) DX = X(2)-A
               DO 85 I = 1,NADD
                  IF (AFIN) THEN
                     X(I+1) = A+DX*Y**((NADD-I+ONE)/NADD)
                  ELSE
                     IF (KCLASS(1) .NE. 1) THEN
                        X(NADD+2-I) = X(NADD+3-I)-
     &                               (X(NADD+4-I)-X(NADD+3-I))*2.4
                     ELSE
                        X(NADD+2-I) = X(NADD+3-I)-
     &                               (X(NADD+4-I)-X(NADD+3-I))
                     ENDIF
                  ENDIF
   85          CONTINUE
            ELSE
               IF (BFIN) DX = B-X(NXINIT-1)
               N = NXINIT-1
               NXINIT = NXINIT+NADD
               X(NXINIT) = B
               DO 90 I = 1,NADD
                  IF (BFIN) THEN
                     X(NXINIT-I) = B-DX*Y**((NADD-I+ONE)/NADD)
                  ELSE
                     IF (KCLASS(2) .NE. 1) THEN
                        X(N+I) = X(N+I-1)+(X(N+I-1)-X(N+I-2))*2.4
                     ELSE
                        X(N+I) = X(N+I-1)+(X(N+I-1)-X(N+I-2))
                     ENDIF
                  ENDIF
   90          CONTINUE
            ENDIF
         ENDIF
   95 CONTINUE
      IF (.NOT. AFIN) X(1) = -ONE/U
      IF (.NOT. BFIN) X(NXINIT) = ONE/U
      HMIN = X(NXINIT)-X(1)
      DO 100 I = 2,NXINIT
         HMIN = MIN(HMIN,X(I)-X(I-1))
  100 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------------
      SUBROUTINE POWER(X,F,N,TOL,IPRINT,EF,CF,OSC,EXACT,Y,IFLAG)
      DOUBLE PRECISION X(*),F(*),TOL,EF,CF,Y(*)
      INTEGER N,IPRINT,IFLAG
      LOGICAL OSC,EXACT
C
C     Find the power function which "dominates" the tabled
C     coefficient function.  The output is Cf and Ef such that
C           f(x)  is asymptotic to  Cf*x^Ef .
C     The vectors X(*) and F(*) hold the N input points:
C           F(I) = f(X(I)) I = 1,...,N.  
C     Set IFLAG = 0 for normal return; 1 for uncertainty in Ef;
C     2 if uncertain about Cf (oscillatory).
C
      DOUBLE PRECISION ERROR,TOLAIT,TOLMIN,  ZERO,HALF,ONE
      INTEGER K,NY
      PARAMETER (ZERO = 0.0, HALF = 0.5D0, ONE = 1.0)
      DATA TOLMIN/1.D-6/
C
C     Estimate the exponent.
C
      OSC = .FALSE.
      NY = N-1
      ERROR = 1.E30
      TOLAIT = MIN(TOLMIN,TOL)
      DO 10 K = 1,NY
         IF ((F(K) .NE. ZERO) .AND. (F(K+1) .NE. ZERO)) THEN
            Y(K) = LOG(ABS(F(K+1)/F(K)))/LOG(ABS(X(K+1)/X(K)))
         ELSE
            Y(K) = ZERO
         ENDIF
  10  CONTINUE
      EF = Y(NY)
      IF (IPRINT .GE. 5) THEN
         WRITE (21,*) ' From POWER; E_k and c_k sequences:'
         WRITE (21,15) (Y(K),K=1,NY)
   15    FORMAT(4D19.10)
         WRITE (21,*)
      ENDIF
      CALL AITKEN(EF,TOLAIT,NY,Y,ERROR)
      K = EF+SIGN(HALF,EF)
      IF (ABS(K-EF) .LE. SQRT(TOL)) EF = K
      IF (ABS(ERROR) .GT. TOL*MAX(ONE,ABS(EF))) THEN
C
C        There is uncertainty in the exponent.
C
         IFLAG = 1
      ENDIF
      IF (ABS(EF) .LE. TOL) EF = ZERO
C
C     Estimate the coefficient.
C
      DO 20 K = 1,N-1
         Y(K) = F(K)/ABS(X(K))**EF
   20 CONTINUE
      CF = Y(N-1)
      IF (IPRINT .GE. 5) WRITE (21,15) (Y(K),K=1,N-1)
      CALL AITKEN(CF,TOLAIT,N-1,Y,ERROR)
      IF ((EF .GT. 20.) .AND. (ABS(CF) .LE. TOL)) THEN
C
C        Coefficient probably has exponential behavior.
C
         CF = SIGN(ONE,Y(N-1))
      ELSE
         IF ((ABS(ERROR) .GT. TOL*MAX(ONE,ABS(CF))) .OR.
     &       ((ABS(F(N)-CF*X(N)**EF) .GT. 20.0*TOL*ABS(F(N))) .AND.
     &        (EF .NE. ZERO ))) THEN
C
C           There is uncertainty in the coefficient; call such
C           cases oscillatory.
C
            IFLAG = 2
            OSC = .TRUE.
         ENDIF
      ENDIF
      IF (ABS(CF) .GT. 1.D7) THEN
         EXACT = .FALSE.
         RETURN
      ENDIF
      K = CF+HALF
      IF ((ABS(K-CF) .LE. SQRT(TOL))  .AND. (K .NE. 0)) CF = K
      EXACT = .TRUE.
      DO 30 K = 1,N
         IF (ABS(F(K)-CF*X(K)**EF) .GT. TOL*ABS(F(K))) EXACT = .FALSE.
   30 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE PQRINT(X,SQRTRP,QLNF)
      DOUBLE PRECISION X,SQRTRP,QLNF
C
C     Evaluate the integrands needed for the asymptotic formulas.
C       (1) The Liouville normal form potential Qlnf:
C             Qlnf(t) =  q/r + f"(t)/f   with   f = (pr)**.25  .
C       (2) The term in the change of independent variable:
C                        sqrt(r/p) .
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &                 EPS,FL,FM,FR,XDOTL,XDOTR,PX,QX,RX,Z,
     &                 ZERO,TWO,FOUR
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2)
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      PARAMETER (ZERO = 0.0, TWO = 2.0, FOUR = 4.0)
      IF (LNF) THEN
         CALL COEFF(X,PX,QX,RX)
         IF ((PX .EQ. ZERO) .OR. (RX .EQ. ZERO) .OR. (PX*RX .LT. ZERO))
     &   THEN
            FLAG = -15
            RETURN
         ENDIF
         NCOEFF = NCOEFF+1
         QLNF = QX/RX
         SQRTRP = SQRT(RX/PX)
      ELSE
         EPS = MIN(1.D-4,MIN(ABS(B-X),ABS(X-A))/TWO)
         Z = X-EPS
         CALL COEFF(Z,PX,QX,RX)
         IF ((PX .EQ. ZERO) .OR. (RX .EQ. ZERO) .OR. (PX*RX .LT. ZERO))
     &   THEN
            FLAG = -15
            RETURN
         ENDIF
         XDOTL = SQRT(PX/RX)
         FL = SQRT(SQRT(PX*RX))
         Z = X+EPS
         CALL COEFF(Z,PX,QX,RX)
         XDOTR = SQRT(PX/RX)
         FR = SQRT(SQRT(PX*RX))
         CALL COEFF(X,PX,QX,RX)
         SQRTRP = SQRT(RX/PX)
         FM = SQRT(SQRT(PX*RX))
         NCOEFF = NCOEFF+3
         QLNF = QX/RX + ((FR-FM)*(XDOTR-XDOTL)/FOUR+
     &                   (FR-TWO*FM+FL)/SQRTRP)/(EPS*EPS*FM*SQRTRP)
         IF (ABS(QLNF) .LE. EPS) QLNF = ZERO
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE REGULR(JOB,JOBMSH,TOL,NEV,EV,IPRINT,NEXTRP,XEF,EF,PDEF,
     &                  HMIN,STORE)
      INTEGER NEV,IPRINT,NEXTRP
      LOGICAL JOB,JOBMSH
      DOUBLE PRECISION TOL(*),EV,XEF(*),EF(*),PDEF(*),HMIN,STORE(*)
***********************************************************************
*                                                                     *
*     REGULR calculates Sturm-Liouville eigenvalue and (optionally)   *
*     eigenfunction estimates for the problem described initially.    *
*                                                                     *
***********************************************************************
C
C    Input parameters:
C      JOB       = logical variable describing tasks to be carried out.
C                  JOB = .True. iff an eigenfunction is to be calculated.
C      JOBMSH    = logical variable, JOBMSH = .True. iff initial mesh
C                  is a function of the eigenvalue index.
C      CONS(*)   = real vector of 8 input constants: A1, A1', A2, A2',
C                  B1, B2, A, B.
C      TOL(*)    = real vector of 6 tolerances. 
C                  TOL(1) is the absolute error tolerance for e-values,
C                  TOL(2) is the relative error tolerance for e-values,
C                  TOL(3) is the abs. error tolerance for e-functions,
C                  TOL(4) is the rel. error tolerance for e-functions,
C                  TOL(5) is the abs. error tolerance for e-function
C                         derivatives,
C                  TOL(6) is the rel. error tolerance for e-function
C                         derivatives.
C                  Eigenfunction tolerances need not be set if JOB is
C                  False.  All absolute error tolerances must be
C                  positive; all relative must be at least 100 times
C                  the unit roundoff.
C      NEV       = integer index for the eigenvalue sought; NEV .GE. 0 .
C      EV        = real initial guess for eigenvalue NEV; accuracy is
C                  not at all critical, but if a good estimate is
C                  available some time may be saved.
C      IPRINT    = integer controlling amount of internal printing done.
C
C    Output parameters:
C      EV        = real computed approximation to NEVth eigenvalue.
C      XEF(*)    = real vector of points for eigenfunction output.
C      EF(*)     = real vector of eigenfunction values: EF(i) is the
C                  estimate of u(XEF(i)).  If JOB is False then this
C                  vector is not referenced.
C      PDEF(*)   = real vector of eigenfunction derivative values:
C                  PDEF(i) is the estimate of (pu')(XEF(i)).  If JOB is
C                  False then this vector is not referenced.
C
C    Auxiliary storage:
C      STORE(*) = real vector of auxiliary storage, must be dimensioned
C                 at least max[100,26N]. (N the number of mesh points)
C
C     Storage allocation in auxiliary vector (currently Maxlvl = 10):
C         STORE(*)
C       1   ->      N           vector of mesh points X(*),
C     N+1   ->     2N           best current eigenfunction values,
C    2N+1   ->     3N           best current derivative values,
C    3N+1   ->     4N           scale factors in GETEF,
C    4N+1   ->  (6+2*Maxlvl)N   intermediate eigenfunction values.
C-----------------------------------------------------------------------
C     Local variables:
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &      I,II,J,JDU,JU,KDU,KK,KU
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &      ABSERR,ALPHA1,ALPHA2,BETA1,BETA2,DELTA,EFNORM,ERROR,
     &      EVEXT(20),EVHAT,EVHIGH,EVLOW,FHIGH,FLOW,H,PDUMAX,QINT,QLNF,
     &      RELERR,RPINT,SQRTRP,TOLEXT,TOLMIN,TOLPDU,TOLSUM,TOL1,TOL2,
     &      TOL3,TOL4,TOL5,TOL6,UMAX,Z,     ASYMEV,
     &      ZERO,HALF,ONE,THREE,FIVE,TEN
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2),
     &        DONE,EFDONE,EFIN(2,2),EVDONE,EXFULL
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, HALF = 0.5D0 ,ONE = 1.0, THREE = 3.0,
     &           FIVE = 5.0, TEN = 10.0, TOLMIN = 1.D-3)
C
      ALPHA1 = ZERO
      ALPHA2 = ONE
      BETA1 = ZERO
      BETA2 = ONE
      IF (NEV .LT. 0) THEN
         FLAG = -37
         RETURN
      ENDIF
      EVDONE = .FALSE.
      TOL1 = MIN(TOL(1),TOLMIN)
      TOL2 = MIN(TOL(2),TOLMIN)
      IF (.NOT. (REG(1) .AND. REG(2))) THEN
         TOL1 = TOL1/THREE
         TOL2 = TOL2/THREE
      ENDIF
      IF (JOB) THEN
         TOL3 = MIN(TOL(3),TOLMIN)
         TOL4 = MIN(TOL(4),TOLMIN)
         TOL5 = MIN(TOL(5),TOLMIN)
         TOL6 = MIN(TOL(6),TOLMIN)
         ABSERR = TOL1/FIVE
         RELERR = TOL2/FIVE
         EXFULL = .TRUE.
      ELSE
         EFDONE = .TRUE.
         ABSERR = TOL1/TEN
         RELERR = TOL2/TEN
         EXFULL = .FALSE.
      ENDIF
      TOLSUM = TOL1+TOL2      
      IF (JOBMSH) THEN
         IF (NEV .GE. 0) THEN
            II = NXINIT+16
            CALL MESH(.TRUE.,NEV,STORE,STORE(II),STORE(2*II+1),
     &                STORE(3*II+1),STORE(4*II+1),TOLSUM,HMIN)
         ENDIF
         IF (IPRINT .GE. 1) THEN
            WRITE (*,15) (STORE(I),I=1,NXINIT)
            WRITE (21,15) (STORE(I),I=1,NXINIT)
  15        FORMAT(' Level 0 mesh:',/,(5G15.6))
         ENDIF
      ENDIF
C
C     Compute estimates for integrals in asymptotic formulas (accuracy
C     is not all that critical).
C
      QINT = ZERO
      RPINT = ZERO
      DO 20 I = 2,NXINIT
         H = STORE(I)-STORE(I-1)
         Z = STORE(I-1)+HALF*H
         IF ((.NOT. AFIN) .AND. (I .EQ. 2)) THEN
            H = STORE(3)-STORE(2)
            Z = STORE(2)
         ENDIF
         IF ((.NOT. BFIN) .AND. (I .EQ. NXINIT)) THEN
            H = STORE(NXINIT-1)-STORE(NXINIT-2)
            Z = STORE(NXINIT-1)
         ENDIF
         CALL PQRINT(Z,SQRTRP,QLNF)
         IF (FLAG .LT. 0) RETURN
         QINT = QINT+H*QLNF
         RPINT = RPINT+H*SQRTRP
   20 CONTINUE
      IF (QINT .GT. ONE/U) QINT = ZERO
      IF (RPINT .GT. ONE/U) RPINT = ZERO
C
C     Loop over the levels.
C
      DO 60 LEVEL = 0,MAXLVL
         IF (HMIN/2**LEVEL .LE. TEN*U) THEN
            FLAG = -8
            GOTO 70
         ENDIF
C
C        Find a bracket for the Nevth eigenvalue.
C
         IF (LEVEL .EQ. 0) THEN
            EV = ASYMEV(NEV,QINT,RPINT,ALPHA1,ALPHA2,BETA1,BETA2)
            EV = MAX(EV,ZERO)
            DELTA = HALF
         ELSE
            DELTA = MAX(TOLSUM*ABS(EVHAT),HALF*HALF*DELTA)
            IF (LEVEL .GT. 1) THEN
               ERROR = (EVEXT(LEVEL)-EVEXT(LEVEL-1))/THREE
               IF (ABS(ERROR) .LE. 100.) THEN
                  EV = EVHAT+ERROR
               ELSE
                  DELTA = ONE
                  EV = EVHAT
               ENDIF
            ELSE
               EV = EVHAT
            ENDIF
         ENDIF
         EVLOW = EV-DELTA
         EVHIGH = EV+DELTA
         IF (IPRINT .GE. 4) WRITE (21,25) EVLOW,EVHIGH
   25    FORMAT('      In  bracket:',2D24.15)
         CALL BRCKET(NEV,EVLOW,EVHIGH,FLOW,FHIGH,ABSERR,RELERR,STORE)
         IF (IPRINT .GE. 4) WRITE (21,30) EVLOW,EVHIGH
   30    FORMAT('      Out bracket:',2D24.15)
         DELTA = HALF*(EVHIGH-EVLOW)
         IF (FLAG .LT. 0) RETURN
         IF (ABS(EVHIGH-EVLOW) .GT. MAX(ABSERR,RELERR*ABS(EVHIGH)))
     &      THEN
            CALL ZZERO(EVLOW,EVHIGH,FLOW,FHIGH,ABSERR,RELERR,J,STORE)
            IF (J .NE. 0) THEN
               FLAG = -7
               RETURN
            ENDIF
         ENDIF
         EVHAT = MIN(EVLOW,EVHIGH)
         IF (IPRINT .GE. 1) THEN
            WRITE (*,40) LEVEL,EVHAT
            WRITE (21,40) LEVEL,EVHAT
   40       FORMAT(' Level ',I3,' ;    EvHat = ',D24.15)
         ENDIF
         EV = EVHAT
         TOLEXT = MAX(TOL1,ABS(EV)*TOL2)
         CALL EXTRAP(EV,TOLEXT,LEVEL+1,NEXTRP,EXFULL,.TRUE.,0,EVEXT,
     &               IPRINT,ERROR,EVDONE)
         IF (JOB) THEN
            CALL GETEF(EVHAT,EFNORM,IPRINT,STORE,efin)
            IF (LEVEL .EQ. 0) THEN
               UMAX = ONE
               PDUMAX = ONE
C
C              Set pointers to STORE(*).
C
               JU = NXINIT
               JDU = 2*NXINIT
               KU = 4*NXINIT
               KDU = (MAXLVL+5)*NXINIT
               KK = MAXLVL+1
            ENDIF
C
C           Extrapolate eigenfunction values.
C
            TOLEXT = MAX(TOL3,UMAX*TOL4)
            TOLPDU = MAX(TOL5,PDUMAX*TOL6)
            EFDONE = .TRUE.
            IF (AFIN) THEN
               UMAX = ABS(STORE(JU+1))
               PDUMAX = ABS(STORE(JDU+1))
            ELSE
               UMAX = ZERO
               PDUMAX = ZERO
            ENDIF
            IF (EFIN(1,1)) THEN
               CALL EXTRAP(STORE(JU+1),TOLEXT,LEVEL+1,NEXTRP,EXFULL,
     &                     .TRUE.,0,STORE(KU+1),0,ERROR,DONE)
               EFDONE = EFDONE .AND. DONE
            ENDIF
            IF (EFIN(2,1)) THEN
               CALL EXTRAP(STORE(JDU+1),TOLPDU,LEVEL+1,NEXTRP,EXFULL,
     &                     .TRUE.,0,STORE(KDU+1),0,ERROR,DONE)
               EFDONE = EFDONE .AND. DONE
            ENDIF
            DO 50 I = 2,NXINIT-1
               II = KK*(I-1)+1
               CALL EXTRAP(STORE(JU+I),TOLEXT,LEVEL+1,NEXTRP,EXFULL,
     &                     .TRUE.,0,STORE(KU+II),0,ERROR,DONE)
               EFDONE = EFDONE .AND. DONE 
               CALL EXTRAP(STORE(JDU+I),TOLPDU,LEVEL+1,NEXTRP,EXFULL,
     &                     .TRUE.,0,STORE(KDU+II),0,ERROR,DONE)
               EFDONE = EFDONE .AND. DONE
               UMAX = MAX(UMAX,ABS(STORE(JU+I)))
               PDUMAX = MAX(PDUMAX,ABS(STORE(JDU+I)))
   50       CONTINUE 
            II = KK*(NXINIT-1)+1
            IF (EFIN(1,2)) THEN
               CALL EXTRAP(STORE(JU+NXINIT),TOLEXT,LEVEL+1,NEXTRP,
     &                     EXFULL,.TRUE.,0,STORE(KU+II),0,ERROR,DONE)
                EFDONE = EFDONE .AND. DONE
            ENDIF
            IF (EFIN(2,2)) THEN
               CALL EXTRAP(STORE(JDU+NXINIT),TOLPDU,LEVEL+1,NEXTRP,
     &                     EXFULL,.TRUE.,0,STORE(KDU+II),0,ERROR,DONE)
               EFDONE = EFDONE .AND. DONE
            ENDIF
            IF (BFIN) THEN
               UMAX = MAX(UMAX,ABS(STORE(JU+NXINIT)))
               PDUMAX = MAX(PDUMAX,ABS(STORE(JDU+NXINIT)))
            ENDIF
            ABSERR = MAX(HALF*ABSERR,TEN*U)
            RELERR = MAX(HALF*RELERR,TEN*U)
         ENDIF
         IF (EVDONE .AND. (LEVEL .GE. 2) .AND. EFDONE) GOTO 70
   60 CONTINUE
      IF (.NOT. EVDONE) THEN
         FLAG = -1
         RETURN
      ENDIF
C
C     Unload eigenfunction values.
C
   70 IF (JOB) THEN
         DO 80 I = 1,NXINIT
            XEF(I) = STORE(I)
            EF(I) = STORE(JU+I)
            PDEF(I) = STORE(JDU+I)
   80    CONTINUE
         IF ((FLAG .GE. 0) .AND. (.NOT. EFDONE)) FLAG = -2
      ENDIF
      RETURN
      END
C---------------------------------------------------------------------
      SUBROUTINE SHOOT(EV,X,MU,FEV)
      INTEGER MU
      DOUBLE PRECISION EV,X(*),FEV
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &        KCLASS(2),I,IMAX,J,K1,K2,KLVL,MODE,NSAVE,NZERO
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2)
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &       CP(2),CR(2),CUTOFF,D(4,2),EMU(2),EP(2),EQLNF(2),
     &       ER(2),ETA(2,2),PNU(2),
     &       CHI,DPSI,DV,H,HALFH,HOMEGA,OMEGA,PDV,PHASE,PN,PSI,RN,
     &       SA1,SA2,SB1,SB2,SCALE,SGN,T,TAU,V,VNEW,XLEFT,X2,Z,
     &       ZERO,HALF,ONE,TWO,PI
      COMMON /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, HALF = 0.5D0, ONE = 1.0, TWO = 2.0,
     &           PI = 3.141592653589793D0)
      DATA IMAX/1000000/
C
C     Make one shot across (A,B) for the current mesh using the scaled
C     variable v(x). Count zeros if COUNTZ is True.
C
C     Shoot from x=A to x=B.
C
      V = A2-A2P*EV
      PDV = A1-A1P*EV
      NZERO = 0
      NSAVE = NSGNF
      SA1 = A1
      SA2 = A2
      SB1 = B1
      SB2 = B2
      KLVL = 2**LEVEL
      MODE = 0
C
C     Modify base sign count if necessary.
C
      IF (((KCLASS(1) .GE. 9) .OR. (KCLASS(2) .GE. 9)) .AND.
     &    (EV .LT. CUTOFF)) THEN
         IF (KCLASS(1) .GE. 9) THEN
            SA1 = D(4,1)
            SA2 = ONE
            V = SA2
            PDV = SA1
            MODE = 1
         ENDIF
         IF (KCLASS(2) .GE. 9) THEN
            SB1 = D(4,2)
            SB2 = ONE
         ENDIF
         SGN = A2*B2
         IF (SGN .EQ. ZERO) THEN
            SGN = SA1*SB2+SA2*SB1
            IF (SGN .EQ. ZERO) SGN = A1*B1
         ENDIF
         NSGNF = SIGN(ONE,SGN)               
      ENDIF
      IF (.NOT. AFIN) MODE = 3
      SCALE = MAX(ABS(V),ABS(PDV))
      V = V/SCALE
      PDV = PDV/SCALE
      DO 20 I = 2,NXINIT
         XLEFT = X(I-1)
         H = X(I)-XLEFT
         IF (MODE .EQ. 1) THEN
            H = (H/D(1,1))**ETA(1,1)
            XLEFT = ZERO
         ENDIF
         IF (MODE .EQ. 3) THEN
            XLEFT = ZERO
            H = -ONE/X(2)
         ENDIF
         IF ((KCLASS(2) .GE. 9) .AND. (I .EQ. NXINIT) .AND.
     &       (EV .LT. CUTOFF)) THEN
C
C           Convert from u(x) to V(t) near x=b.
C
            MODE = 2
            T = (H/D(1,2))**ETA(1,2)
            CHI = D(2,2)*T**ETA(2,2)
            V = V/CHI
            PN = CP(2)*H**(EP(2)-ONE)
            PDV = (PDV/(PN*CHI*ETA(1,2))+ETA(2,2)*V)*T**(TWO*EMU(2)-ONE)
            H = T
            XLEFT = -H
         ENDIF
         IF ((.NOT. BFIN) .AND. (I .EQ. NXINIT)) THEN
            MODE = 4
            H = ONE/X(I-1)
            XLEFT = -H
         ENDIF
         H = H/KLVL
         HALFH = HALF*H
         DO 10 J = 1,KLVL
            Z = XLEFT+HALFH
            CALL STEP(Z,H,EV,PN,RN,TAU,OMEGA,HOMEGA,PSI,DPSI,SCALE,MODE)
            IF (FLAG .LT. 0) THEN
               FEV = ZERO
               RETURN
            ENDIF
            DV = PDV/PN
            VNEW = DPSI*V+PSI*DV
            XLEFT = XLEFT+H
            IF (COUNTZ) THEN
C
C              Count zeros of v(x).
C
               IF (TAU .LE. ZERO) THEN
                  IF (VNEW*V .LT. ZERO) NZERO = NZERO+1
               ELSE
                  IF (DV .EQ. ZERO) THEN
                     NZERO = NZERO+INT(HALF+HOMEGA/PI)
                  ELSE
                     PHASE = ATAN(V*OMEGA/DV)
                     K1 = PHASE/PI
                     X2 = (PHASE+HOMEGA)/PI
                     IF (X2 .LT. IMAX) THEN
                        K2 = X2
                        NZERO = NZERO+K2-K1
                     ELSE
                        NZERO = IMAX
                     ENDIF
                     IF (PHASE*(PHASE+HOMEGA) .LT. ZERO) NZERO = NZERO+1
                  ENDIF
                  NZERO = MIN(IMAX,NZERO)
               ENDIF
            ENDIF
            PDV = -PN*TAU*PSI*V+DPSI*PDV
            V = VNEW
   10    CONTINUE
         IF (MODE .EQ. 1) THEN
C
C           Convert from V(t) back to u(x) near x=a.
C
            PN = CP(1)*(X(2)-A)**(EP(1)-ONE)
            T = ((X(2)-A)/D(1,1))**ETA(1,1)
            CHI = D(2,1)*T**ETA(2,1)
            PDV = PN*ETA(1,1)*CHI*(ETA(2,1)*V+PDV*T**(ONE-TWO*EMU(1)))
            V = CHI*V
         ENDIF
         MODE = 0
   20 CONTINUE
      FEV = SB1*V+SB2*PDV
      IF (COUNTZ) THEN
C
C        Adjust zero count.
C
         MU = NZERO
         IF (A2P .NE. ZERO) THEN
            IF (EV .GE. SA2/A2P) MU = NZERO+1
         ENDIF
         IF (SB2 .NE. ZERO) THEN
            SGN = SIGN(ONE,NSGNF*FEV)
            IF (MOD(MU,2) .EQ. 1) SGN = -SGN
            IF (SGN .LT. ZERO) MU = MU+1
         ENDIF
      ENDIF
      NSGNF = NSAVE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE START(JOB,CONS,TOL,NEV,INDXEV,N,XEF,NUMT,T,NEXTRP,X)
      INTEGER NEV,INDXEV(*),N,NUMT,NEXTRP
      LOGICAL JOB(*)
      DOUBLE PRECISION CONS(*),TOL(*),XEF(*),T(*),X(*)
C
C     This routine tests the input data, initializes the labeled
C     common blocks, and generates the first mesh.  Check
C        eigenfunction tolerances iff JOB(1) is True ,
C        XEF(*) iff JOB(2) is True,
C        NUMT, T(*) iff JOB(3) is True ,
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &        I,K
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2)
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &                 UFLOW,URN,
     &                 ZERO,HALF,ONE,PIOVR2,TWO,FIVE,HUNDRD
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
      PARAMETER (ZERO = 0.0, HALF = 0.5D0, ONE = 1.0,
     &           PIOVR2 = 1.5707963267948966D0, TWO = 2.0, FIVE = 5.0,
     &           HUNDRD = 100.0)
************************************************************************
*   In the following DATA statement, initialize URN to an estimate for *
*   the unit roundoff and UFLOW to a value somewhat less than          *
*   -ln(underflow).  E.g.,                                             *
*      for IEEE double precision use URN = 2.D-16, UFLOW = 650;        *
*      for VAXen double precision use URN = 1.D-17, UFLOW = 85;        *
*      for Crays (single precision) use URN = 7.D-15, UFLOW = 5000.    *
*   Exact values are not at all critical.                              *
*   Here, we assume IEEE double precision:                             *
*                                                                      *
      DATA URN/2.D-16/,UFLOW/650.0/
C***********************************************************************
      U = URN
      UNDER = UFLOW
C
C     Initialize.
C
      A1 = CONS(1)
      A1P = CONS(2)
      A2 = CONS(3)
      A2P = CONS(4)
      A = CONS(7)
      B1 = CONS(5)
      B2 = CONS(6)
      B = CONS(8)
      NCOEFF = 0
      FLAG = 0
C
C     Test input.
C
      IF (AFIN .AND. BFIN .AND. (A .GE. B)) FLAG = -34
      IF (TOL(1) .LE. ZERO) FLAG = -35
      IF (TOL(2) .LT. HUNDRD*U) FLAG = -36
      IF (JOB(1)) THEN
         IF (TOL(3) .LE. ZERO) FLAG = -35
         IF (TOL(4) .LT. HUNDRD*U) FLAG = -36
         IF (TOL(5) .LE. ZERO) FLAG = -35
         IF (TOL(6) .LT. HUNDRD*U) FLAG = -36
      ENDIF
      IF (JOB(2)) THEN
         IF (N .EQ. 1) THEN
            FLAG = -30
         ELSE
            IF (AFIN .AND. (XEF(2) .LE. A)) FLAG = -39
            DO 10 I = 3,N-1
              IF (XEF(I-1) .GE. XEF(I)) FLAG = -39
   10       CONTINUE
            IF (BFIN .AND. (XEF(N-1) .GT. B)) FLAG = -39
         ENDIF
         IF ((.NOT. AFIN) .AND. (XEF(2) .GE. ZERO)) FLAG = -39
         IF ((.NOT. BFIN) .AND. (XEF(N-1) .LE. ZERO)) FLAG = -39
      ENDIF 
      IF ((JOB(2) .OR. JOB(3)) .AND. (NEV .GT. 0)) THEN
         DO 20 I = 1,NEV
            IF (INDXEV(I) .LT. 0) FLAG = -37
   20    CONTINUE
      ENDIF
      IF (JOB(3)) THEN
         IF (NUMT .LE. 0) FLAG = -38
         DO 30 I = 2,NUMT
            IF (T(I) .LE. T(I-1)) FLAG = -39
   30    CONTINUE
      ENDIF
      IF (FLAG .LT. 0) RETURN
C
C     Set MAXEXT, the maximum number of extrapolations allowed, and
C         MAXINT, the maximum number of intervals in X(*) allowed when
C         mesh is chosen by START.
C
C     IMPORTANT NOTE: the size of various fixed arrays in this package
C     depends on the value of MAXEXT in this FORTRAN77 implementation.
C     If MAXEXT is increased, then more storage may have to be allocated
C     to the columns of R(*,*) in EXTRAP.
C
      MAXEXT = 6
      MAXINT = 31
C
C     Calculate maximum number of columns in extrapolation table
C     and the maximum number of levels allowed.
C
      K = -LOG10(TOL(2))
      I = MAX(K+3,0)
      NEXTRP = MIN(MAX(3,I/2),MAXEXT)
C
C     Calculate the initial mesh.
C
      IF (N .GT. 0) THEN
         NXINIT = N
      ELSE
         NXINIT = MIN(2*NEXTRP+3,MAXINT)
         IF (JOB(1)) N = NXINIT
      ENDIF
      IF (JOB(2)) THEN
         A = XEF(1)
         B = XEF(NXINIT)
         DO 40 I = 1,NXINIT
            X(I) = XEF(I)
   40    CONTINUE
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE STEP(X,H,EV,PX,RX,TAU,OMEGA,HOMEGA,PSI,DPSI,SCLOG,MODE)
      INTEGER MODE
      DOUBLE PRECISION X,H,EV,PX,RX,TAU,OMEGA,HOMEGA,PSI,DPSI,SCLOG
C
C    Evaluate the coefficient functions, the scaled basis function PSI,
C    its derivative DPSI, and the log of the scale factor SCLOG.
C
      INTEGER FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT,
     &        KCLASS(2)
      LOGICAL AFIN,BFIN,COUNTZ,LFLAG(6),LNF,LC(2),OSC(2),REG(2)
      DOUBLE PRECISION A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER,
     &        CP(2),CR(2),CUTOFF,D(4,2),EMU(2),EP(2),EQLNF(2),
     &        ER(2),ETA(2,2),PNU(2)
      COMMON /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      COMMON /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      COMMON /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      COMMON /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C
      DOUBLE PRECISION DX,FP,FR,OVER,QX,T,TMU,Z,
     &            ZERO,HNDRTH,QUART,HALF,ONE,TWO,SIX,TWELVE,TWENTY
      PARAMETER (ZERO = 0.0, HNDRTH = .01, QUART = 0.25D0,
     &           HALF = 0.5D0, ONE = 1.0, TWO = 2.0, SIX = 6.0,
     &           TWELVE = 12.0, TWENTY = 20.0)
      DATA OVER/1.D8/
C
C    Evaluate the coefficient functions at X and calculate TAU.  The
C    error flag FLAG is zero for a successful calculation; if p(x)
C    or r(x) are zero, then FLAG is set to -15.  If the argument for
C    a trig function exceeds OVER then FLAG is set to -10.
C         Proceed normally when Mode = 0; when Mode = 1 or 2 use the
C    change of variable for "hard" problems; when Mode = 3 or 4 use
C    the change of variable t = -1/x near infinity.
C
      IF (MODE .EQ. 0) THEN
         CALL COEFF(X,PX,QX,RX)
      ELSE
         IF (MODE .GE. 3) THEN
            T = X
            X = -ONE/T
            CALL COEFF(X,PX,QX,RX)
            T = T*T
            PX = T*PX
            QX = QX/T
            RX = RX/T
         ELSE
            T = X
            IF (ETA(1,MODE) .EQ. ONE) THEN
               DX = D(1,MODE)*ABS(T)
            ELSE
               DX = D(1,MODE)*(ABS(T))**(ONE/ETA(1,MODE))
            ENDIF
            IF (MODE .EQ. 1) THEN
               Z = A+DX
            ELSE
               Z = B-DX
            ENDIF
            CALL COEFF(Z,PX,QX,RX)
         ENDIF
      ENDIF
      NCOEFF = NCOEFF+1
      IF ((PX .EQ. ZERO) .OR. (RX .EQ. ZERO)) THEN
         FLAG = -15
         RETURN
      ENDIF
      IF ((MODE .EQ. 1) .OR. (MODE .EQ. 2)) THEN
         IF (EMU(MODE) .EQ. HALF) THEN
            TMU = ABS(T)
         ELSE
            TMU = (T*T)**EMU(MODE)
         ENDIF
         FP = CP(MODE)
         IF (EP(MODE) .NE. ZERO) FP = FP*DX**EP(MODE)
         FR = CR(MODE)
         IF (ER(MODE) .NE. ZERO) FR = FR*DX**ER(MODE)
         PX = PX*TMU/FP
         QX = (QX/FR - D(3,MODE)/T**2)*TMU
         RX = RX*TMU/FR
      ENDIF
      TAU = (EV*RX-QX)/PX
      OMEGA = SQRT(ABS(TAU))
      HOMEGA = H*OMEGA
      SCLOG = ZERO
C
C     Evaluate the scaled basis functions.
C
      IF (HOMEGA .GT. HNDRTH) THEN
         IF (TAU .GT. ZERO) THEN
            IF (HOMEGA . GT. OVER) THEN
               FLAG = -10
               RETURN
            ENDIF
            DPSI = COS(HOMEGA)
            PSI = SIN(HOMEGA)/OMEGA
         ELSE
            SCLOG = HOMEGA
            IF (HOMEGA .LT. UNDER) THEN
               T = TANH(HOMEGA)
               DPSI = ONE/(ONE+T)
               PSI = T*DPSI/OMEGA
            ELSE
               SCLOG = MIN(SCLOG,TWO*UNDER)
               DPSI = HALF
               PSI = DPSI/OMEGA
            ENDIF
         ENDIF
      ELSE
         T = TAU*H*H
         DPSI = ONE+T*(T/TWELVE-ONE)/TWO
         PSI = H*(ONE+T*(T/TWENTY-ONE)/SIX)
         IF (T .LT. ZERO) THEN
            T = MAX(ABS(PSI),ABS(DPSI))
            SCLOG = LOG(T)
            PSI = PSI/T
            DPSI = DPSI/T
         ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE ZZERO(B,C,FB,FC,ABSERR,RELERR,IFLAG,X)
C
      DOUBLE PRECISION B,C,FB,FC,ABSERR,RELERR,X(*)
      INTEGER IFLAG
C
C  ZZERO computes a root of F.  The method used is a combination of
C  bisection and the secant rule.  This code is adapted from one in
C  the text "Foundations of Numerical Computing" written by Allen,
C  Pruess, and Shampine.
C
C  Input parameters:
C     B,C   = values of X such that F(B)*F(C) .LE. 0.
C     FB,FC = values of F at input B and C, resp.
C     ABSERR,RELERR = absolute and relative error tolerances.  The
C             stopping criterion is:
C               ABS(B-C) .LE. 2.0*MAX(ABSERR,ABS(B)*RELERR).
C  Output parameters:
C     B,C   = see IFLAG returns.
C     FB    = value of final residual F(B).
C     IFLAG = 0 for normal return; F(B)*F(C) .LT. 0 and the
C               stopping criterion is met (or F(B)=0).  B always
C               satisfies ABS(F(B)) .LE. ABS(F(C)).
C           = 1 if too many function evaluations were made; in this version
C               200 are allowed.
C           =-2 if F(B)*F(C) is positive on input.
C
C  Local variables:
      DOUBLE PRECISION A,ACMB,CMB,FA,P,Q,TOL,WIDTH
      INTEGER KOUNT,MAXF,MU,NF
C
C  Internal constants
C
      DOUBLE PRECISION ZERO,ONE,TWO,EIGHT
      PARAMETER (ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0,
     &           MAXF = 200)
C
C  Initialization.
C
      KOUNT = 0
      WIDTH = ABS(B-C)
      A = C
      FA = FC
      IF (SIGN(ONE,FA) .EQ. SIGN(ONE,FB)) THEN
         IFLAG = -2
         RETURN
      ENDIF
      FC = FA
      NF = 2
   20 IF (ABS(FC) .LT. ABS(FB)) THEN
C
C  Interchange B and C so that ABS(F(B)) .LE. ABS(F(C)).
C
         A = B
         FA = FB
         B = C
         FB = FC
         C = A
         FC = FA
      ENDIF
      CMB = (C-B)/TWO
      ACMB = ABS(CMB)
      TOL = MAX(ABSERR,ABS(B)*RELERR)
C
C  Test stopping criterion and function count.
C
      IF (ACMB .LE. TOL) THEN
         IFLAG = 0
         RETURN
      ENDIF
      IF (NF .GE. MAXF) THEN
         IFLAG = 1
         RETURN
      ENDIF
C
C  Calculate new iterate implicitly as B+P/Q where we arrange
C     P .GE. 0.  The implicit form is used to prevent overflow.
C
      P = (B-A)*FB
      Q = FA-FB
      IF (P .LT. ZERO) THEN
         P = -P
         Q = -Q
      ENDIF
C
C  Update A; check if reduction in the size of bracketing interval is
C     satisfactory.  If not, bisect until it is.
C
      A = B
      FA = FB
      KOUNT = KOUNT+1
      IF (KOUNT .GE. 4) THEN
         IF (EIGHT*ACMB .GE. WIDTH) THEN
            B = B+CMB
            GOTO 30
         ENDIF
         KOUNT = 0
         WIDTH = ACMB
      ENDIF
C
C  Test for too small a change.
C
      IF (P .LE. ABS(Q)*TOL) THEN
C
C  Increment by tolerance.
C
         B = B+SIGN(TOL,CMB)
      ELSE
C
C  Root ought to be between B and (C+B)/2.
C
         IF (P .LT. CMB*Q) THEN
C
C  Use secant rule.
C
            B = B+P/Q
         ELSE
C
C  Use bisection.
C
            B = B+CMB
         ENDIF
      ENDIF
C
C  Have completed computation for new iterate B.
C
   30 CALL SHOOT(B,X,MU,FB)
      NF = NF+1
      IF (ABS(FB) .EQ. ZERO) THEN
         IFLAG = 0
         C = B
         FC = FB
         RETURN
      ENDIF
      IF (SIGN(ONE,FB) .EQ. SIGN(ONE,FC)) THEN
         C = A
         FC = FA
      ENDIF
      GOTO 20
      END
