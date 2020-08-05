C====>==================================================================
C====> @Fortran-file{
C====>    author        = "Werner Rath",
C====>    version       = "1.0",
C====>    date          = "14 July 1995",
C====>    time          = "10:50:44 MESZ",
C====>    filename      = "demo.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "09736 389 1643 12415",
C====>    email         = "rath@mathematik.tu-chemnitz.de",
C====>    codetable     = "ISO/ASCII",
C====>    keywords      = "",
C====>    supported     = "yes",
C====>    abstract      = "",
C====>    docstring     = "The checksum field above contains a CRC-16
C====>                     checksum as the first value, followed by the
C====>                     equivalent of the standard UNIX wc (word
C====>                     count) utility output of lines, words, and
C====>                     characters.  This is produced by Robert
C====>                     Solovay's checksum utility.",
C====> }
C====>==================================================================
C     DEMO: Demonstration program for DGELDA.
C
C     The problem is a 2x2 strangeness index 1 DAE wich is equivalent
C     to a purely algebraic equation.
C
C     REVISIONS
C
C     1995, July 14 (W. Rath).
C
C     *****************************************************************
C
C     .. Parameters ..
      INTEGER          NEQ, LRW, LIW
      PARAMETER        (NEQ =2, LRW = 5000, LIW = 100)
      INTEGER          KPRINT
      PARAMETER        (KPRINT = 3)
      INTEGER          LUN
      PARAMETER        (LUN = 6)
      INTEGER          NOUT
      PARAMETER        (NOUT = 10)
      DOUBLE PRECISION TSTART
      PARAMETER        (TSTART = 0.0D0)
C     .. Local Scalars ..
      DOUBLE PRECISION DSEC, ERO, HU, T, TOUT
      INTEGER          I, IDID, IOUT, IWARN, METHOD, NFE, NJE, NQU, NST
C     .. Local Arrays ..
      DOUBLE PRECISION ATOL(1), DTOUT(NOUT), ERROR(NEQ), RPAR(1),
     $                 RTOL(1), RWORK(LRW), X(NEQ), XPRIME(NEQ)
      INTEGER          CVAL(4), INFO(20), IPAR(1), IWORK(LIW)
C     .. External Functions ..
      DOUBLE PRECISION DNRM2, DSECND
      EXTERNAL         DNRM2, DSECND
      EXTERNAL         EDIF, ADIF, FDIF
C     .. Data Statements ..
      DATA             DTOUT / 0.1D0, 0.2D0, 0.3D0, 0.4D0, 0.5D0,
     $                         0.6D0, 0.7D0, 0.8D0, 0.9D0, 1.0D0 /
C     .. Executable Statements ..
C
C     First, we have to choose between BDF and RK method.
C
   1  PRINT *,'Choose method: (1 = BDF, 2 = RK)'
      READ *, METHOD
      IF (METHOD .LT. 1 .OR. METHOD .GT. 2) THEN
         PRINT *,'ERROR: Method out of range !!'
         PRINT *,'Try again'
         GOTO 1
      ENDIF
C
C     Set the INFO array to tell the code how to solve the problem.
C     We use the standard setup by setting all entries to zero.
C
      DO 10 I=1,20
         INFO(I) = 0
  10  CONTINUE
C
C     Set the starting time and the initial values.
C     Note, that we choose inconsistent initial values.
C
      T = TSTART
      X(1)=0.0D0
      X(2)=0.0D0
C
C     Now we have to set the tolerances for DGELDA to indicate how
C     accurate we want the solution to be computed.
C     In this case we use a combined error test with the same absolute
C     and relative tolerance.
C
      ATOL(1) = 1.0D-5
      RTOL(1) = 1.0D-5
C
C     Write some information about the problem to be solved.
C
      IF (KPRINT .GE. 2) THEN
         WRITE (LUN,100) NEQ,RTOL,ATOL,METHOD
      ENDIF
C
 100  FORMAT(/1X,' DEMONSTRATION PROGRAM FOR DGELDA',///
     $     1X,' PROBLEM 1 FROM P. KUNKEL AND V. MEHRMANN,',/
     $     1X,' A NEW CLASS OF DISCRETIZATION METHODS FOR THE',/
     $     1X,' SOLUTION OF LINEAR DIFFERENTIAL-ALGEBRAIC EQUATIONS:',//
     $     1X,' NEQ =',I3,   /1X,' RTOL =',E10.1,'   ATOL =',E10.1,/
     $     1X,' METHOD =',I3,' (1 = BDF, 2 = RK, -1 = DASSL)')
C
C     Before we solve the problem, we compute the characteristic values
C     of the DAE and consistent initial values by setting TOUT = T and
C     calling DGELDA.
C
      TOUT = T
      CALL DGELDA(EDIF, ADIF, FDIF,
     $     NEQ, T, TOUT, X, XPRIME, CVAL,
     $     IPAR, RPAR, IWORK, LIW, RWORK, LRW,
     $     RTOL, ATOL, METHOD, INFO, IWARN, IDID)
      IF (IDID.LT.0) THEN
         CALL DGEERM(EDIF, ADIF, FDIF,
     $        NEQ, T, TOUT, X, XPRIME, CVAL,
     $        IPAR, RPAR, IWORK, LIW, RWORK, LRW,
     $        RTOL, ATOL, METHOD, INFO, IWARN, IDID)
         STOP
      ENDIF
C
C     Now we can print some statistics about the characteristic values
C     of the DAE and the computed consistent values.
C
C
      IF (KPRINT .GT. 2) THEN
         WRITE (LUN,115) (CVAL(I),I=1,4)
         WRITE (LUN,116) (X(I),I=1,NEQ)
      ENDIF
C
 115  FORMAT(//' PROBLEM DESCRIPTION'/
     $     '   STRANGENESS INDEX   ',I2,4X,
     $     '   DIFFERENTIAL COMPONENTS   ',I2/29X,
     $     '   ALGEBRAIC    COMPONENTS   ',I2/29X,
     $     '   UNDETERMINED COMPONENTS   ',I2/)
 116  FORMAT(' CORRECTED INITIAL VALUES'/(' ',6E12.4/))
C
C     The next step is to solve the problem.
C     We will use DGELDA to compute NOUT intermediate solutions from
C     DTOUT(1) to DTOUT(10) (see above).
C
C     For each intermediate solution the weighted error is computed and
C     some statistics are displayed, where
C     T    is the actual time,
C     X(1) is the computed first solution component at time T,
C     ERO  is the actual weighted error,
C     ORD  is the order of the BDF method used in the last step (if METHOD=1),
C     H    is the stepsize used in the last step.
C
      IF (KPRINT .GT. 2) THEN
         WRITE (LUN,110)
      ENDIF
C
 110  FORMAT(///
     $     10X,'T',14X,'X(1)',12X,'ERO',8X,'ORD',8X,'H'/)
      DSEC = DSECND()
      DO 200 IOUT = 1,NOUT
         TOUT = DTOUT(IOUT)
 120     CALL DGELDA(EDIF, ADIF, FDIF,
     $        NEQ, T, TOUT, X, XPRIME, CVAL,
     $        IPAR, RPAR, IWORK, LIW, RWORK, LRW,
     $        RTOL, ATOL, METHOD, INFO, IWARN, IDID)
         IF (IDID.LT.0) THEN
            CALL DGEERM(EDIF, ADIF, FDIF,
     $        NEQ, T, TOUT, X, XPRIME, CVAL,
     $        IPAR, RPAR, IWORK, LIW, RWORK, LRW,
     $        RTOL, ATOL, METHOD, INFO, IWARN, IDID)
            STOP
         ENDIF
         ERROR(1) = ((1+t)*EXP(-t) - X(1))
     $              /(ABS((1+t)*EXP(-t))+1)
         ERROR(2) = (EXP(-t) - X(2))
     $              /(EXP(-t)+1)
         ERO = DMAX1(ERO,DNRM2(NEQ,ERROR,1)/DSQRT(DFLOAT(NEQ)))
         HU = RWORK(7)
         NQU = IWORK(8)
         IF (KPRINT .GT. 2) THEN
            WRITE (LUN,130)
     $           T,X(1),DNRM2(NEQ,ERROR,1)/DSQRT(DFLOAT(NEQ)),NQU,HU
         ENDIF
C
 130     FORMAT(1X,E15.5,E16.5,E16.5,I6,E14.3)
C
 200  CONTINUE
C
C     Finally, we display some final statistics
C
      DSEC = DSECND()-DSEC
      NST = IWORK(11)
      NFE = IWORK(12)
      NJE = IWORK(13)
      IF (KPRINT .GT. 2) THEN
         WRITE (LUN,210) NST,NFE,NJE,IWORK(14),ERO,DSEC
      ENDIF
C
 210  FORMAT(//1X,' FINAL STATISTICS FOR THIS RUN..',/
     $     1X,' NUMBER OF STEPS               =',I5/
     $     1X,' NUMBER OF EVALUATIONS         =',I5/
     $     1X,' NUMBER OF FACTORIZATIONS      =',I5/
     $     1X,' NUMBER OF ERROR TEST FAILURES =',I5/
     $     1X,' MAX ERROR                     =',E10.2/
     $     1X,' RUNTIME                       =',E10.2)
C
      STOP
      END
      SUBROUTINE EDIF(N,T,IDIF,E,LDE,IPAR,RPAR,IERR)
C
C     ARGUMENT LIST
C
C        ARGUMENTS IN
C
C             N - INTEGER.
C             T - DOUBLE PRECISION.
C          IDIF - INTEGER.
C           LDE - INTEGER.
C                 The leading dimension of array E as declared in the
C                 calling program.
C                 LDE .GE. N
C
C        ARGUMENTS OUT
C
C             E - DOUBLE PRECISION array of DIMENSION (LDE,N).
C                 The leading N by N part of this array contains the
C                 IDIF-th derivative of E(t) at time T.
C
C        ERROR INDICATOR
C
C          IERR - INTEGER.
C                 Unless the routine detects an error (see next section),
C                 IERR contains 0 on exit.
C
C        WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C        IERR = -1 : Failed to compute the IDIF-th derivative.
C        IERR = -2 : On entry, IDIF is larger than the highest
C                    derivative of E(t) the subroutine provides
C
C     REVISIONS
C
C     1995, July 14 (W. Rath).
C
C     *****************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION T, RPAR
      INTEGER          IDIF, IERR, IPAR, LDE, N
C     .. Array Arguments ..
      DOUBLE PRECISION E(LDE,*)
C     .. Local Scalars ..
      INTEGER          I, J
C     .. Executable Statements ..
      IERR = 0
      DO 110 I=1,N
        DO 100 J=1,N
           E(I,J)=0.0D0
 100    CONTINUE
 110  CONTINUE
      IF (IDIF .EQ. 0) THEN
         E(2,1)=1.0D0
         E(2,2)=-T
      ELSEIF (IDIF .EQ. 1) THEN
         E(2,2)=-1.0D0
      ENDIF
      RETURN
C *** Last line of EDIF ***
      END
      SUBROUTINE ADIF(N,T,IDIF,A,LDA,IPAR,RPAR,IERR)
C
C     ARGUMENT LIST
C
C        ARGUMENTS IN
C
C             N - INTEGER.
C             T - DOUBLE PRECISION.
C          IDIF - INTEGER.
C           LDA - INTEGER.
C                 The leading dimension of array A as declared in the
C                 calling program.
C                 LDA .GE. N
C
C        ARGUMENTS OUT
C
C             A - DOUBLE PRECISION array of DIMENSION (N,N).
C                 The leading N by N part of this array contains the
C                 IDIF-th derivative of A(t) at time T.
C
C        ERROR INDICATOR
C
C          IERR - INTEGER.
C                 Unless the routine detects an error (see next section),
C                 IERR contains 0 on exit.
C
C        WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C        IERR = -1 : Failed to compute the IDIF-th derivative.
C        IERR = -2 : On entry, IDIF is larger than the highest
C                    derivative of E(t) the subroutine provides
C
C     REVISIONS
C
C     1995, July 14 (W. Rath).
C
C     *****************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION RPAR, T
      INTEGER          IDIF, IERR, IPAR, LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)
C     .. Local Scalars ..
      INTEGER          I, J
C     .. Executable Statements ..
      IERR = 0
      DO 110 I=1,N
         DO 100 J=1,N
            A(I,J)=0.D0
 100     CONTINUE
 110  CONTINUE
      IF (IDIF .EQ. 0) THEN
         A(1,1)=-1.0D0
         A(1,2)=T
      ELSEIF (IDIF .EQ. 1) THEN
         A(1,2)=1.0D0
      ENDIF
      RETURN
C *** Last line of ADIF ***
      END
      SUBROUTINE FDIF(N,T,IDIF,F,IPAR,RPAR,IERR)
C
C     ARGUMENT LIST
C
C        ARGUMENTS IN
C
C             N - INTEGER.
C             T - DOUBLE PRECISION.
C          IDIF - INTEGER.
C
C        ARGUMENTS OUT
C
C             F - DOUBLE PRECISION array of DIMENSION N.
C                 The first N components of this array contain the
C                 IDIF-th derivative of f(t) at time T.
C
C        ERROR INDICATOR
C
C          IERR - INTEGER.
C                 Unless the routine detects an error (see next section),
C                 IERR contains 0 on exit.
C
C        WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C        IERR = -1 : Failed to compute the IDIF-th derivative.
C        IERR = -2 : On entry, IDIF is larger than the highest
C                    derivative of f(t) the subroutine provides
C
C     REVISIONS
C
C     1995, July 14 (W. Rath).
C
C     *****************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION T, RPAR
      INTEGER          IDIF, IERR, IPAR, N
C     .. Local Scalars ..
      INTEGER          I
C     .. Local Arrays ..
      DOUBLE PRECISION F(*)
C     .. Executable Statements ..
      IERR = 0
      DO 100 I=1,N
         F(I)=0.0D0
 100  CONTINUE
      IF (IDIF .EQ. 0 .OR. IDIF .EQ. 2) THEN
         F(1)=DEXP(-T)
      ELSEIF (IDIF .EQ. 1 .OR. IDIF .EQ. 3) THEN
         F(1)=-DEXP(-T)
      ELSE
         IERR = 1
      ENDIF
      RETURN
C *** Last line of SUBROUTINE FDIF ***
      END
