C   Copyright (C) 2020  Peter Kunkel, Volker Mehrmann, Werner Rath, Jörg Weickert
C
C   This program is free software: you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation, either version 3 of the License, or
C   (at your option) any later version.
C
C   This program is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program.  If not, see <https://www.gnu.org/licenses/>.
C
C   source: https://dx.doi.org/10.5281/zenodo.3972144
C
C====>==================================================================
C====> @Fortran-file{
C====>    author        = "Werner Rath",
C====>    version       = "1.1.01",
C====>    date          = "21 July 1995",
C====>    time          = "16:56:52 MESZ",
C====>    filename      = "dbdstp.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "45179 599 2300 18832",
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
      SUBROUTINE DBDSTP (EDIF, ADIF, FDIF, NEQ, T, H, HOLD, HMIN,
     $                   CJ, CJOLD, K, KOLD, M, ID, IA, IU, IREQ, LMAX,
     $                   NS, IPHASE, ROWEQU, COLEQU, JSTART, X, XPRIME,
     $                   E, LDE, A, LDA, Z1, LDZ1, Z2Q, LDZ2Q, WT, W,
     $                   LDW, PHI, LDPHI, ALPHA, BETA, GAMMA, PSI,
     $                   SIGMA, DELTA, ERRV, RS, CS, F, EQ, LDEQ, AQ,
     $                   LDAQ, FQ, AH, LDAH, IPAR, RPAR, IWM, WORK,
     $                   LWORK, MCONST, IWARN, IDID)
C
C     PURPOSE
C
C     Performs one step of the BDF integration of GELDA [2].
C
C     METHOD
C
C     DBDSTP solves a system of differential-algebraic equations of
C     the form E(T) XPRIME(T) = A(T) X(T) + F(T), for one step
C     (normally from T to T+H).
C
C     The methods used are modified divided difference, fixed leading
C     coefficient forms of backward differentiation formulas. The code
C     adjusts the stepsize and order to control the local error per
C     step (see [1] for details). Most of the parameters contain
C     information which is needed internally by DBDSTP to continue from
C     step to step, so DBDSTP cannot be used alone.
C
C     DBDSTP is a modiefied version of the DDASTP subroutine of DASSL
C     [3].
C
C     REFERENCES
C
C     [1] K. E. Brenan and S. L. Campbell and L. R. Petzold.
C         Numerical Solution of Initial-Value Problems in Differential
C         Algebraic Equations.
C         Elsevier, North Holland, New York, 1989.
C
C     [2] P. Kunkel, V. Mehrmann, W. Rath and J. Weickert.
C         GELDA: A software package for the solution of general linear
C         differential algebraic equations.
C         Preprint SPC 95_8, TU Chemnitz-Zwickau, February 1995.
C
C     [3] L. R. Petzold.
C         A Description of DASSL: A Differential/Algebraic System
C         Solver.
C         Scientific Computing, R. S. Stepleman et al. (Eds.),
C         North-Holland, pp. 65-68, 1983.
C
C     CONTRIBUTORS
C
C     W. Rath (TU Chemnitz, Germany).
C
C     REVISIONS
C
C     1995, July 21 [Version 1.1.1]
C       There was a nonstandard use of the INT function, i.e. converting
C       a LOGICAL into an INTEGER. Now this is done by hand. (W. Rath)
C
C     1995, July 16 [Version 1.1]
C       Changed order of IPAR, RPAR, and IWORK, RWORK to meet
C       SLICOT interface standard.
C       Added equilibration of linear systems for the BDF solver,
C       therefore we need additional space in RWORK (2*NEQ). (W. Rath)
C
C     1995, July 10 [Version 1.00]
C       First release. (W. Rath)
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER        (ONE=1.0D0)
      INTEGER          OCTF, OETF, OIPVT, OMXORD, ONEV, ONFA, ONST
      PARAMETER        (OMXORD=3, ONST=11, ONEV=12, ONFA=13, OETF=14,
     $                  OCTF=15, OIPVT=21)
C     .. Subroutine Arguments ..
      EXTERNAL         ADIF, EDIF, FDIF
C     .. Scalar Arguments ..
      DOUBLE PRECISION CJ, CJOLD, H, HOLD, HMIN, T
      INTEGER          COLEQU, IA, ID, IDID, IPHASE, IU, IWARN, JSTART,
     $                 K, KOLD, LDA, LDAH, LDAQ, LDE, LDEQ, LDPHI, LDW,
     $                 LDZ1, LDZ2Q, LMAX, LWORK,  NEQ, NS, M, ROWEQU
      LOGICAL          COL, MCONST, ROW
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), AH(LDAH,*), ALPHA(*), AQ(LDAQ,*),
     $                 BETA(*), CS(*), DELTA(*), E(LDE,*), EQ(LDEQ,*),
     $                 ERRV(*), F(*), FQ(*), GAMMA(*), PHI(LDPHI,*),
     $                 PSI(*), RS(*), RPAR(*), SIGMA(*), X(*),
     $                 XPRIME(*), W(LDW,*), WORK(*),WT(*),
     $                 Z1(LDZ1,*), Z2Q(LDZ2Q,*)
      INTEGER          IPAR(*), IWM(*)
C     .. Local Scalars ..
      DOUBLE PRECISION AMAX, ALPHA0, ALPHAS, CK, COLCND, ENORM, ERK,
     $                 ERKM1,ERKM2,ERKP1, ERR, EST, HNEW, R, ROWCND,
     $                 TEMP1, TEMP2 ,TERK,TERKM1,TERKM2, TERKP1, TOLD
      INTEGER          I, IAQ, IDQ, IER, IRANK, IRED, IREQ, IUQ, J,
     $                 KDIFF, KM1, KNEW, KP1, KP2, MQ, NCF,
     $                 NEF, NSF, NSP1
      LOGICAL          FACTOR
      CHARACTER*1      EQUED
C     .. External Functions ..
      DOUBLE PRECISION DGENRM
      LOGICAL          LSAME
      EXTERNAL         DGENRM, LSAME
C     .. External Subroutines ..
      EXTERNAL         DAXPY, DBDTRP, DGEEQU, DGETRF, DGETRS, DSCAL
C     .. Executable Statements ..
C
C----------------------------------------------------------------------
C     Block 1.
C     Initialize. On the first call, set
C     the order to 1 and initialize
C     other variables.
C----------------------------------------------------------------------
C
C     Initializations for all calls.
      IDID=1
      TOLD=T
      NCF=0
      NSF=0
      NEF=0
      IRED=0
C
C     If this is the first step, perform
C     other initializations.
      IF(JSTART .NE. 1) THEN
        K=1
        KOLD=0
        HOLD=0.0D0
        PSI(1)=H
        CJOLD = 0.0D0
	FACTOR = .TRUE.
        IPHASE = 0
        NS=0
        JSTART=1
      ENDIF
C
C----------------------------------------------------------------------
C     Block 2
C     Compute coefficients of formulas for
C     this step.
C----------------------------------------------------------------------
C
 200  CONTINUE
      KP1=K+1
      KP2=K+2
      KM1=K-1
      TOLD=T
      IF (H.NE.HOLD.OR.K .NE. KOLD) NS = 0
      NS=MIN(NS+1,KOLD+2)
      NSP1=NS+1
C
      IF (KP1 .GE. NS) THEN
        BETA(1)=1.0D0
        ALPHA(1)=1.0D0
        TEMP1=H
        GAMMA(1)=0.0D0
        SIGMA(1)=1.0D0
        DO 210 I=2,KP1
           TEMP2=PSI(I-1)
           PSI(I-1)=TEMP1
           BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
           TEMP1=TEMP2+H
           ALPHA(I)=H/TEMP1
           SIGMA(I)=(I-1)*SIGMA(I-1)*ALPHA(I)
           GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
 210    CONTINUE
        PSI(KP1)=TEMP1
      ENDIF
C
C     Compute ALPHAS, ALPHA0.
      ALPHAS = 0.0D0
      DO 220 I=1,K
        ALPHAS = ALPHAS - 1.0D0/I
 220  CONTINUE
      ALPHA0 = 0.0D0
      DO 230 I = 1,K
        ALPHA0 = ALPHA0 - ALPHA(I)
 230  CONTINUE
C
C     Compute leading coefficient CJ
C     NOTE: If MCONST = TRUE and IU = 0 we can use the old factorization
C           if CJ = CJOLD.
      CJ = -ALPHAS/H
      IF (MCONST .AND. IU.EQ.0) THEN
        IF (CJ.EQ.CJOLD) THEN
          FACTOR = .FALSE.
        ELSE
          FACTOR = .TRUE.
        ENDIF
      ENDIF
C
C     Compute variable stepsize error coefficient CK.
      CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0)
      CK = MAX(CK,ALPHA(KP1))
C
C     Change PHI to PHI STAR.
      IF (KP1 .GE. NSP1) THEN
        DO 240 J=NSP1,KP1
	  CALL DSCAL(NEQ,BETA(J),PHI(1,J),1)
 240   CONTINUE
      ENDIF
C
C     Update time.
      T=T+H
C
C----------------------------------------------------------------------
C     Block 3
C     Predict the solution and derivative,
C     and solve the corrector equation.
C----------------------------------------------------------------------
C
C     First, predict the solution and derivative.
 300  CONTINUE
C
C     Set X(*)      = PHI(*,1)
C         XPRIME(*) = 0.0D0
      CALL DCOPY(NEQ,PHI(1,1),1,X,1)
      DO 310 I=1,NEQ
        XPRIME(I)=0.0D0
 310  CONTINUE
C
C     FOR J=2,KP1
C            X(*)     =X(*)+PHI(*,J)
C            XPRIME(*)=XPRIME(*)+GAMMA(J)*PHI(*,J)
      DO 320 J=2,KP1
        CALL DAXPY(NEQ,ONE,PHI(1,J),1,X,1)
        CALL DAXPY(NEQ,GAMMA(J),PHI(1,J),1,XPRIME,1)
 320  CONTINUE
C
C     Compute E(T), A(T), F(T) and set up linear system.
      IWM(ONEV) = IWM(ONEV) + 1
      CALL DNFIX1(EDIF, ADIF, FDIF, NEQ, T, LMAX, M, ID, IA, IU, IREQ,
     $            MQ, IDQ, IAQ, IUQ, E, LDE, A, LDA, F, EQ, LDEQ, AQ,
     $            LDAQ, FQ, Z1, LDZ1, Z2Q, LDZ2Q, AH, LDAH, IPAR, RPAR,
     $            WORK, LWORK, MCONST, IRED)
      IF (IRED .LT. 0) GOTO 380
      IF (FACTOR) THEN
        DO 340 I = 1,NEQ
          DO 330 J = 1,NEQ
	     W(I,J) = -CJ*E(I,J) + A(I,J)
 330      CONTINUE
 340    CONTINUE
      ENDIF
      CALL DCOPY(NEQ,XPRIME,1,ERRV,1)
      CALL DAXPY(NEQ,-CJ,X,1,ERRV,1)
      CALL DCOPY(NEQ,F,1,DELTA,1)
      CALL DGEMV('N',NEQ,NEQ,1.0D0,E,LDE,ERRV,1,-1.0D0,DELTA,1)
C
C     First case: IU = 0 and we have to solve a linear system.
      IF (IU.EQ.0) THEN
         IF (FACTOR) THEN
            CALL DGEEQU(NEQ,NEQ,W,NEQ,RS,CS,ROWCND,COLCND,AMAX,IER)
            CALL DLAQGE(NEQ,NEQ,W,NEQ,RS,CS,ROWCND,COLCND,AMAX,EQUED)
            ROW = (LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' ))
            IF (ROW) THEN
               ROWEQU = 1
            ELSE
               ROWEQU = 0
            ENDIF
            COL = (LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' ))
            IF (COL) THEN
               COLEQU = 1
            ELSE
               COLEQU = 0
            ENDIF
            CALL DGETRF(NEQ,NEQ,W,NEQ,IWM(OIPVT),IER)
            IF (IER .NE. 0) GOTO 380
            CJOLD = CJ
            IWM(ONFA) = IWM(ONFA) + 1
         ENDIF
         IF(ROWEQU.EQ.1) THEN
            DO 930 I = 1, NEQ
               DELTA(I) = RS(I)*DELTA(I)
 930        CONTINUE
         ENDIF
         CALL DGETRS('N',NEQ,1,W,NEQ,IWM(OIPVT),DELTA,NEQ,IER)
         IF(COLEQU.EQ.1) THEN
            DO 970 I = 1, NEQ
               DELTA(I) = CS(I)*DELTA(I)
 970        CONTINUE
         ENDIF
C
C           Update ERRV, X and XPRIME:
C           If solution of the system is stored in DELTA (not in ERRV)
C           we get.
C             ERRV   = DELTA - X
C             X      = DELTA
C             XPRIME = XPRIME * CJ * (DELTA - X)
         CALL DCOPY(NEQ,DELTA,1,ERRV,1)
         CALL DAXPY(NEQ,-1.0D0,X,1,ERRV,1)
         CALL DCOPY(NEQ,DELTA,1,X,1)
         CALL DAXPY(NEQ,CJ,ERRV,1,XPRIME,1)
C
C     Second case: IU > 0 and we have to solve a least sqaure problem.
C
      ELSE
         CALL DCOPY(NEQ,DELTA,1,ERRV,1)
C        In PHI(:,1) we have the old X.

         CALL DGEMV('N',NEQ,NEQ,-1.0D0,W,NEQ,PHI(1,1),1,
     $        1.0D0,ERRV,1)
C
C        Use SVD for the LS problem.
         CALL DGELSS(ID+IA,NEQ,1,W,NEQ,ERRV,NEQ,DELTA,-1.0D0,IRANK,
     $               WORK,LWORK,IER)
         IF (IER.NE.0) GOTO 380
C
C          Update ERRV,X and XPRIME:
C          If solution of the ls problem is stored in B (not in ERRV) we get
C             ERRV   = B
C             X      = X + B
C             XPRIME = XPRIME + CJ * B
C
         CALL DAXPY(NEQ,1.0D0,PHI(1,1),1,ERRV,1)
         CALL DAXPY(NEQ,-1.0D0,X,1,ERRV,1)
         CALL DAXPY(NEQ,1.0D0,ERRV,1,X,1)
         CALL DAXPY(NEQ,CJ,ERRV,1,XPRIME,1)
      ENDIF
 380  IF((IER .NE. 0) .OR. (IRED .LT. 0)) GOTO 600
C
C----------------------------------------------------------------------
C     Block 4
C     Estimate the errors at orders K,K-1,K-2
C     As if constant stepsize was used. Estimate
C     the local error at order K and test
C     whether the current step is successful.
C----------------------------------------------------------------------
C
C     Estimate errors at orders K,K-1,K-2.
      ENORM = DGENRM(NEQ,ERRV,WT)
      ERK = SIGMA(K+1)*ENORM
      TERK = (K+1)*ERK
      EST = ERK
      KNEW=K
      IF(K .NE. 1) THEN
C
C       Set DELTA(*) = PHI(*,KP1) + ERRV(*).
        CALL DCOPY(NEQ,PHI(1,KP1),1,DELTA,1)
        CALL DAXPY(NEQ,ONE,ERRV,1,DELTA,1)
        ERKM1=SIGMA(K)*DGENRM(NEQ,DELTA,WT)
        TERKM1 = K*ERKM1
        IF(K .LE. 2) THEN
          IF(TERKM1 .LE. 0.5D0*TERK) THEN
C
C           Lower the ordder.
            KNEW=K-1
            EST = ERKM1
	  ENDIF
        ELSE
C
C         Set DELTA(*) = PHI(*,K) + DELTA(*).
          CALL DAXPY(NEQ,ONE,PHI(1,K),1,DELTA,1)
          ERKM2=SIGMA(K-1)*DGENRM(NEQ,DELTA,WT)
          TERKM2 = (K-1)*ERKM2
          IF(MAX(TERKM1,TERKM2).LE.TERK) THEN
C
C           Lower the order.
            KNEW=K-1
            EST = ERKM1
	  ENDIF
	ENDIF
      ENDIF
C
C
C     Calculate the local error for the current step
C     to see if the step was successful.
      ERR = CK * ENORM
      IF(ERR .GT. 1.0D0)GO TO 600
C
C----------------------------------------------------------------------
C     BLOCK 5
C     The step is successful. Determine
C     the best order and stepsize for
C     the next step. Update the differences
C     for the next step.
C----------------------------------------------------------------------
      IDID=1
      IWM(ONST)=IWM(ONST)+1
      KDIFF=K-KOLD
      KOLD=K
      HOLD=H
C
C
C     Estimate the error at order K+1 unless:
C        Already decided to lower order, or
C        already using maximum order, or
C        stepsize not constant, or
C        order raised in previous step.
      IF (KNEW.EQ.KM1.OR.K.EQ.IWM(OMXORD)) IPHASE=1
      IF (IPHASE .NE. 0) THEN
        IF (KNEW.EQ.KM1) GOTO 540
        IF (K.EQ.IWM(OMXORD)) GOTO 550
        IF (KP1.GE.NS.OR.KDIFF.EQ.1) GOTO 550
C
C       Set DELTA(*)=ERRV(*)-PHI(*,KP2).
        CALL DCOPY(NEQ,ERRV,1,DELTA,1)
        CALL DAXPY(NEQ,-ONE,PHI(1,KP2),1,DELTA,1)
        ERKP1 = (1.0D0/(K+2))*DGENRM(NEQ,DELTA,WT)
        TERKP1 = (K+2)*ERKP1
        IF(K.LE.1) THEN
          IF(TERKP1.GE.0.5D0*TERK)GO TO 550
	ELSE
          IF(TERKM1.LE.MIN(TERK,TERKP1))GO TO 540
          IF(TERKP1.GE.TERK.OR.K.EQ.IWM(OMXORD))GO TO 550
	ENDIF
C
C       Raise order.
        K=KP1
        EST = ERKP1
        GOTO 550
C
C       Lower order.
 540    K=KM1
        EST = ERKM1
        GOTO 550
C
C
C       Determine the appropriate stepsize for
C       the next step.
 550    HNEW=H
        TEMP2=K+1
        R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
        IF (R .GE. 2.0D0) THEN
          HNEW = 2.0D0*H
        ELSEIF(R .LE. 1.0D0) THEN
          R = MAX(0.5D0,MIN(0.9D0,R))
          HNEW = H*R
        ENDIF
        H=HNEW
      ELSE
C
C       If IPHASE = 0, increase order by one and multiply stepsize by
C       factor two.
        K = KP1
        HNEW = H*2.0D0
        H = HNEW
      ENDIF
C
C
C     Update differences for next step.
      IF (KOLD.NE.IWM(OMXORD)) CALL DCOPY(NEQ,ERRV,1,PHI(1,KP2),1)
C
C     Set PHI(*,KP1)=PHI(*,KP1)+ERRV(*).
      CALL DAXPY(NEQ,ONE,ERRV,1,PHI(1,KP1),1)
C
C     FOR J=KP1-1,1,-1 DO
C       SET PHI(*,J)=PHI(*,J)+PHI(*,J+1).
      DO 560 J=KP1-1,1,-1
         CALL DAXPY(NEQ,ONE,PHI(1,J+1),1,PHI(1,J),1)
 560  CONTINUE
      RETURN
C
C----------------------------------------------------------------------
C     Block 6
C     The step is unsuccessful. Restore X, PSI, PHI
C     determine appropriate stepsize for
C     continuing the integration, or exit with
C     an error flag if there have been many
C     failures.
C----------------------------------------------------------------------
 600  IPHASE = 1
C
C     restore T, PHI, PSI.
      T=TOLD
      IF(KP1.GE.NSP1) THEN
C
C       FOR J=NSP1,KP1 DO
C         SET PHI(*,J)=PHI(*,J)/BETA(J).
        DO 610 J=NSP1,KP1
           TEMP1 = ONE/BETA(J)
           CALL DSCAL(NEQ,TEMP1,PHI(1,J),1)
 610    CONTINUE
      ENDIF
      DO 620 I=2,KP1
        PSI(I-1)=PSI(I)-H
 620  CONTINUE
C
C
C     Test whether failure is due to problems with the
C     linear or least sqaure system, or the error test.
      IF((IER .EQ. 0).AND.(IRED.EQ.0))GO TO 660
      IWM(OCTF)=IWM(OCTF)+1
      IF(IER.EQ.0) GOTO 650
C
C
C     The lineare system is (maybe) singular. We got
C     an error in the factorization or solution subroutine
C     Reduce the stepsize by a factor of 4. If
C     this happens three times in a row on
C     the same step, return with an error flag.
      NSF=NSF+1
      R = 0.25D0
      H=H*R
      IF (NSF .LT. 3 .AND. ABS(H) .GE. HMIN) GO TO 690
      IDID=-8
      GO TO 675
C
C
C     The computation of the solution of the system failed
C     for a reason other than a singular system. If IRED=-2,
C     then return. Otherwise, if IRED=-1 reduce the stepsize
C     and try again, unless too many failures have occured.
C     In all other cases we must return control to the calling
C     program.
 650  CONTINUE
      IF (IRED .EQ. -1) THEN
         NCF = NCF + 1
         R = 0.25D0
         H = H*R
         IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GOTO 690
         IDID = -10
      ELSEIF (IRED .EQ. -2) THEN
         IDID = -21
      ELSE
         IDID = IRED
          IF (IRED .EQ. -26) THEN
            M = MQ
            ID = IDQ
            IA = IAQ
            IU = IUQ
         ENDIF
      ENDIF
      GOTO 675
C
C
C     The cause of the failure was the error estimate
C     exceeding the tolerance.
 660  NEF=NEF+1
      IWM(OETF)=IWM(OETF)+1
C
C     On first error test failure, keep current order or lower
C     order by one.  Compute new stepsize based on differences
C     of the solution.
      IF (NEF .EQ. 1) THEN
        K = KNEW
        TEMP2 = K + 1
        R = 0.90D0*(2.0D0*EST+0.0001D0)**(-ONE/TEMP2)
        R = MAX(0.25D0,MIN(0.9D0,R))
        H = H*R
C
C     On second error test failure, use the current order or
C     decrease order by one.  Reduce the stepsize by a factor of
C     four.
      ELSEIF (NEF .EQ. 2) THEN
        K = KNEW
        H = 0.25D0*H
C
C     On third and subsequent error test failures, set the order to
C     one and reduce the stepsize by a factor of four.
      ELSE
        K = 1
        H = 0.25D0*H
      ENDIF
C
C     Check if stepsize is to small.
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
C
C
C
C
C     For all crashes, restore X to its last value,
C     interpolate to find XPRIME at last T, and return.
 675  CONTINUE
      CALL DBDTRP(T,T,NEQ,K,PHI,LDPHI,PSI,X,XPRIME)
      RETURN
C
C
C     go back and try this step again.
 690  GO TO 200
C *** Last line of DBDSTP ***
      END
