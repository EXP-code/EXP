         SUBROUTINE DTORH1(Z,M,NMAX,PL,QL,NEWN)         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  INPUT :                                                         C         
C    Z        ARGUMENT OF THE FUNCTIONS                            C       
C    M        ORDER OF THE FUNCTIONS                               C                                                                            
C    NMAX     MAXIMUM DEGREE OF THE FUNCTIONS :                    C       
C             WE GET  FUNCTIONS OF ALL THE ORDERS BELOW            C       
C             MIN(NEWN,NMAX). NEWN IS DEFINED BELOW .              C       
C  OUTPUT :                                                        C       
C   *IF MODE IS EQUAL TO 0:                                        C                                                                               
C    PL(N)                                                         C       
C             THESE VALUES ARE KEPT IN AN ARRAY                    C       
C    QL(N)                                                         C       
C             THESE VALUES ARE KEPT IN AN ARRAY                    C       
C                                                                  C       
C    NEWN     MAXIMUM ORDER OF FUNCTIONS CALCULATED WHEN           C        
C             PL (NMAX+1) IS LARGER THAN 1/TINY                    C       
C             (OVERFLOW LIMIT = 1/TINY, TINY IS DEFINED BELOW)     C       
C    NOTE1: FOR A PRECISION OF 10**(-12), IF Z>5 AND (Z/M)>0.22    C
C           THE CODE USES A SERIES EXPANSION FOR PL(0).            C
C           WHEN Z<20 AND (Z/M)<0.22 A CONTINUED FRACTION IS       C
C           APPLIED.                                               C
C    NOTE2: FOR A PRECISION OF 10**(-8), IF Z>5 AND (Z/M)>0.12 THE C
C           CODE USES A SERIES EXPANSION FOR PL(0).                C
C           WHEN Z<20 AND (Z/M)<0.12 A CONTINUED FRACTION IS       C
C           APPLIED.                                               C
C   *IF MODE IS EQUAL TO 1:                                        C       
C      THE SET OF FUNCTIONS EVALUATED IS:                          C       
C           PL(N)/GAMMA(M+1/2),QL(N)/GAMMA(M+1/2),                 C       
C      WHICH ARE RESPECTIVELY STORED IN THE ARRAYS PL(N),QL(N)     C                                                                                   
C      NEWN REFERS TO THIS NEW SET OF FUNCTIONS                    C                                                                                                                                      
C      NOTE1 AND NOTE2 ALSO APPLY IN THIS CASE                     C                                                            
C   *IF MODE IS EQUAL TO 2:                                        C       
C       THE CODE PERFORMS AS FOR MODE 1, BUT THE RESTRICTION       C
C       Z<20 FOR THE EVALUATION OF THE CONTINUED FRACTION IS NOT   C
C       CONSIDERED                                                 C
C       WARNING: USE ONLY IF HIGH M'S FOR Z>20 ARE REQUIRED. THE   C       
C       EVALUATION OF THE CF MAY FAIL TO CONVERGE FOR TOO HIGH Z'S C       
C  PARAMETERS:                                                     C       
C   MODE: ENABLES THE ENLARGEMENT OF THE RANGE OF ORDERS AND       C
C         DEGREES THAT CAN BE EVALUATED.                           C                    
C   EPS:  CONTROLS THE ACCURACY OF THE CONTINUED FRACTIONS         C
C         AND SERIES.                                              C
C   IPRE: REQUIRED PRECISION IN THE EVALUATION OF TOROIDAL         C
C         HARMONICS.                                               C
C           *IF IPRE=1, PRECISION=10**(-12) (TAKING EPS<10**(-12)) C                                                                                
C           *IF IPRE=2, PRECISION=10**(-8)  (TAKING EPS<10**(-8))  C       
C   TINY: SMALL PARAMETER NEAR THE UNDERFLOW LIMIT.                C       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                                         
C   DECLARATION OF VARIABLES   C                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       INTEGER M,NMAX,NMAXP,MODE,IPRE,MP,NP,N,ICAL,NEWN,I
       DOUBLE PRECISION Z,PI,EPS,TINY,OVER,TINYSQ,QZ,PISQ,DPPI,
     *   FL,CC,AR,GAMMA,FC,QDC1,QARGU,ARGU1,DFACQS,FCP,DD,
     *   GAMMAH,ELLIP1,ELLIP2,D1,QM0,DFAC3,DFAC4,PL0,FACTCO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C THE DIMENSION OF QLMM (INTERNAL ARRAY) MUST BE GREATER THAN M  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
       DOUBLE PRECISION  QLMM(0:1001),PL(0:NMAX+1),QL(0:NMAX+1),PR(2)
       PARAMETER(PI=3.14159265358979323D0,EPS=1.D-14,TINY=1.D-290,
     *           MODE=0,IPRE=1)
       OVER=1.D0/TINY
       TINYSQ=DSQRT(TINY)
       IF ((IPRE.NE.1).AND.(IPRE.NE.2)) THEN
         WRITE(6,*)'IPRE MUST BE 1 OR 2'
         STOP
       END IF
       PR(1)=.22D0
       PR(2)=.12D0
       NMAXP=NMAX
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   EPS: REQUIRED ACCURACY FOR THE CONTINUED FRACTION     C
C        (MODIFIED LENTZ)                                 C
C   TINY: SMALL PARAMETER TO PREVENT OVERFLOWS IN THE CF  C
C          (CLOSE TO THE UNDERFLOW LIMIT)                 C                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (Z.LE.1.D0) THEN
          WRITE(6,*)'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
          STOP
       END IF
       QZ=Z
       PISQ=DSQRT(PI)
       DPPI=DSQRT(2.D0)/PISQ
       FL=M/2.D0
       CC=ABS(FLOAT(INT(FL))-FL)
         IF (CC.LT.0.4D0) THEN
           AR=1.
         ELSE
           AR=-1.
         END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  WE CHOOSE EXPANSION OR CF FOR PL(0) DEPENDING ON THE VALUES C
CC  OF Z,M AND MODE                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       ICAL=1
       IF (M.NE.0) THEN
        IF ((Z/M).GT.PR(IPRE)) ICAL=2
        IF (Z.LT.5.D0) THEN
          ICAL=1
        ELSE IF (Z.GT.20.D0) THEN
          IF ((MODE.NE.2).AND.(ICAL.EQ.1)) ICAL=0
        END IF
        IF (ICAL.EQ.0) THEN
          WRITE(*,*)'YOU MUST CHOOSE MODE=2'
          STOP
        END IF
       ELSE
        IF (Z.LT.5.D0) THEN
          ICAL=1
        ELSE
          ICAL=2
        END IF
       END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                                       C
C   WE USE THE CODE IF NMAX IS GREATER THAN OR EQUAL TO 2  C               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF(NMAXP.LT.2) NMAXP=2
       IF (MODE.EQ.0) THEN
          GAMMA=GAMMAH(M,OVER)*AR*PISQ  
          IF (ABS(GAMMA).LT.TINY) THEN
             WRITE(*,*)'M IS TOO LARGE FOR MODE=0'
             WRITE(*,*)'BETTER TRY MODE=1'
             STOP
          END IF              
        ELSE 
          GAMMA=AR
        END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      WE EVALUATE THE CONTINUED FRACTION  USING  C             
C      LENTZ-THOMPSON                             C             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        CALL FRAC(Z,M,0,EPS,TINYSQ,FC)
        QDC1=QZ*QZ-1.D0
        QARGU=QZ/DSQRT(QDC1)
        DFAC1=DPPI*GAMMA/PI
        DFAC2=GAMMA/DPPI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC   WE EVALUATE Q_{-1/2},Q^{1}_{-1/2}              C    
CC   USING SLATEC ROUTINES FOR ELLIPTIC FUNCTIONS   C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        ARGU1=DSQRT(2.D0/(Z+1.D0))
        QLMM(0)=ARGU1*ELLIP1(ARGU1)
        QLMM(1)=-1.D0/DSQRT(2.D0*(QZ-1.D0))
     *          *ELLIP2(ARGU1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    WE APPLY FORWARD RECURRENCE IN M FOR Q'S  C          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        MP=1
         IF (MODE.EQ.0) THEN
1          IF ((MP.LE.M).AND.(ABS(QLMM(MP)).LT.OVER)) THEN
            QLMM(MP+1)=-2.D0*MP*QARGU*QLMM(MP) 
     *      -(MP-0.5D0)*(MP-0.5D0)*QLMM(MP-1) 
             MP=MP+1
            GOTO 1 
           ENDIF       
           IF ((MP-1).NE.M) THEN
            WRITE(*,*)'M IS TOO LARGE FOR MODE=0'
            WRITE(*,*)'BETTER TRY MODE=1'
            STOP
           END IF              
        ELSE
            QLMM(0)=QLMM(0)/PISQ
            QLMM(1)=QLMM(1)*2.D0/PISQ
2           IF ((MP.LE.M).AND.(ABS(QLMM(MP)).LT.OVER)) THEN
             D1=MP+0.5D0
             QLMM(MP+1)=-2.D0*MP*QARGU*QLMM(MP)/D1
     *       -(MP-0.5D0)*QLMM(MP-1)/D1               
             MP=MP+1
             GOTO 2
            ENDIF
            IF ((MP-1).NE.M) THEN
             WRITE(*,*)'M IS TOO LARGE FOR MODE=1,2'
             STOP
            END IF              
        END IF
        NMMAX=M
        DFACQS=-(GAMMA/QLMM(M+1))*GAMMA/PI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  EVALUATION OF PL(0)                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IF (ICAL.EQ.1) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
CC   WE CALCULATE THE CF FOR P'S                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
           CALL FRACPS(QARGU,M+1,0,EPS,TINYSQ,FCP) 
           DD=M+0.5D0
           IF (MODE.NE.0) THEN
            FCP=FCP/DD
            DFACQS=DFACQS/DD
           END IF
           PL(0)=DFACQS/DSQRT(QDC1)/(1.D0-FCP*QLMM(M)/QLMM(M+1))
        ELSE
           CALL EXPAN(Z,MODE,IPRE,OVER,QARGU,M,PL0)
           PL(0)=PL0
        END IF
        QM0=QLMM(M)
        DFAC3=(GAMMA/QM0)*GAMMA/PI/(0.5D0-M)
        PL(1)=PL(0)*FC+DFAC3
        NP=1 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                                         C
C   WE USE THE RECURRENCE RELATIONS FOR P'S   C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
3       IF ((NP.LE.NMAXP).AND.(ABS((NP-M+0.5D0)*PL(NP)).LT.OVER)) THEN
         PL(NP+1)=(2.D0*NP*Z*PL(NP)-(NP+M-0.5D0)*PL(NP-1))/(NP-M+0.5D0)
         NP=NP+1
         GOTO 3
        ENDIF      
        NMAXP=NP-1
        NEWN=NMAXP
        DFAC4=(FACTCO(NMAXP,PL(NMAXP+1),M)*GAMMA)*GAMMA/PI/
     *        (NMAXP+M+0.5D0)     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  WE EVALUATE THE CONTINUED FRACTION USING LENTZ-THOMPSON  C                                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        CALL FRAC(Z,M,NMAXP,EPS,TINYSQ,FC)         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     EVALUATION OF PL(NMAX+1) AND PL(NMAX) USING   C                    
C     THE WRONSKIAN W{PL(NMAX),QL(NMAX)},           C                    
C     THE KNOWN VALUES OF PL(NMAX+1) AND PL(NMAX)   C                    
C     THE VALUE OF H = QL(NMAX+1)/QL(NMAX)          C                    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       QL(NMAXP)=DFAC4/(1.D0-FC*PL(NMAXP)/PL(NMAXP+1))
       QL(NMAXP+1)=QL(NMAXP)*FC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   WE USE THE BACKWARD RECURRENCE RELATION   C                                                                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DO 4 I=1,NMAXP
          NP=NMAXP+1-I
          N=NP-1
          QL(N)=((NP+NP)*Z*QL(NP)-(NP-M+0.5D0)*QL(NP+1))/(NP+M-0.5D0)
4       CONTINUE
        RETURN
        END
       
