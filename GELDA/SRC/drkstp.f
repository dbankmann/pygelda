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
C====>    author        = "Joerg Weickert",
C====>    version       = "1.01",
C====>    date          = "19 July 1995",
C====>    time          = "12:55:00 MESZ",
C====>    filename      = "drkstp.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "46852 427 1602 14209",
C====>    email         = "weickert@mathematik.tu-chemnitz.de",
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
      SUBROUTINE DRKSTP(EDIF, ADIF, FDIF, N, N3, T, H, HOLD, HMIN,
     $                  ERRACC, UROUND, SAFE, FACL, FACR, QUOT1, QUOT2,
     $                  M, ID, IA, IU, IREQ, LMAX, NSING, JSTART, LAST,
     $                  X, XPRIME, E0, LDE0, A0, LDA0, E, LDE, A, LDA,
     $                  EQ, LDEQ, AQ, LDAQ, Z, LDZ, ZQ, LDZQ,
     $                  AH, LDAH, W, LDW, CONT, LDCONT, WT, Z1, Z2, Z3,
     $                  Y, B, F, FQ, IPAR, RPAR, IWORK, WORK, LWORK,
     $                  PRED, MCONST, IWARN, IDID)
C
C     PURPOSE
C
C     Performs one step  of the Runge-Kutta integration of GELDA [2].
C
C     METHOD
C
C     DRKSTP solves a system  of  differential-algebraic equations of
C     the form     E(T) XPRIME(T) = A(T) X(T) + F(T),    for one step
C     (normally from T to T+H).
C
C     The method used is a three-stage implicit Runge Kutta scheme of
C     fifth order accuracy.  The code adjusts the stepsize to control
C     the local error per step (see [1] for details). Most of the pa-
C     rameters contain information  which is needed internally by the
C     subroutine DRKSTP to continue from step to step, so DRKSTP can-
C     not be used alone.
C
C     DRKSTP is  a modified version  of the Runge-Kutta system solver
C     RADAU5 written by Hairer and Wanner [1].
C
C     REFERENCES
C
C     [1] E. Hairer and G. Wanner.
C         Solving Ordinary Differential Equations II.
C         Springer-Verlag, Berlin 1991.
C
C     [2] P. Kunkel, V. Mehrmann, W. Rath and J. Weickert.
C         GELDA: A software package for the solution of general linear
C         differential algebraic equations.
C         Preprint SPC 95\_8, TU Chemnitz-Zwickau, February 1995.
C
C     CONTRIBUTORS
C
C     J. Weickert (TU Chemnitz, Germany).
C
C     REVISIONS
C
C     1995, July 19 [Version 1.1]
C       Changed order of IPAR, RPAR, and IWORK, RWORK to meet
C       SLICOT interface standard.
C
C     1995, July 10 [Version 1.0]
C       First release.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      INTEGER OCTF, OETF, OIPVT, ONEV, ONFA, ONST
      PARAMETER (ONST   = 11)
      PARAMETER (ONEV   = 12)
      PARAMETER (ONFA   = 13)
      PARAMETER (OETF   = 14)
      PARAMETER (OCTF   = 15)
      PARAMETER (OIPVT  = 21)
C     .. Scalar Arguments ..
      DOUBLE PRECISION FACL, FACR, H, HMIN, HOLD,
     $                 QUOT1, QUOT2, SAFE, UROUND, T
      INTEGER          IA, ID, IDID, IREQ, IU, IWARN, JSTART,
     $                 LDA, LDAH, LDAQ, LDA0, LDCONT, LDE, LDEQ,
     $                 LDE0, LDW, LDZ, LDZQ,
     $                 LMAX, LWORK, M, N, N3
      LOGICAL          MCONST, PRED
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), AH(LDAH,*), AQ(LDAQ,*),
     $                 A0(LDA0,*), B(*), CONT(LDCONT,*), E(LDE,*),
     $                 EQ(LDEQ,*), E0(LDE0,*), F(*), FQ(*),
     $                 RPAR(*),
     $                 W(LDW,* ), WORK(*), WT(*),
     $                 X(*), XPRIME(*), Y(*), Z(LDZ,*),
     $                 ZQ(LDZQ,*), Z1(*), Z2(*), Z3(*)
      INTEGER          IPAR(*), IWORK(*)
C     .. User Supplied Functions ..
      EXTERNAL         ADIF, EDIF, FDIF
C     .. Local Scalars ..
      DOUBLE PRECISION A11, A12, A13, A21, A22, A23, A31, A32, A33,
     $                 C1, C1M1, C1MC2, C2, C2M1,
     $                 DD1, DD2, DD3, ERR, ERRACC, FAC, FACGUS, FAC1,
     $                 HA11, HA12, HA13, HA21, HA22, HA23, HA31,
     $                 HA32, HA33, HEE1, HEE2, HEE3, HHFAC, HNEW,
     $                 POSNEG, QUOT, QT, SQ6, U1, TPH
      INTEGER          I, IAQ, IDQ, IER, IRANK, IUQ, J, MQ,
     $                 NSING, N2
      LOGICAL          FACTOR, LAST, REJECT
C     .. External Subroutines ..
      EXTERNAL         DAXPY, DCOPY, DNFIX1, DGEMV, DGELSX,
     $                 DRKCXS, DRKRLS, DSCAL
C     .. Executable Statements ..
C
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** ***
      IF (JSTART .EQ. 0) THEN
         IWORK(OETF) = 0
         IWORK(OCTF) = 0
         NSING = 0
      END IF
C
C ---------- Constants ---------
      SQ6=DSQRT(6.D0)
      C1=(4.D0-SQ6)/10.D0
      C2=(4.D0+SQ6)/10.D0
      C1M1=C1-1.D0
      C2M1=C2-1.D0
      C1MC2=C1-C2
      DD1=-(13.D0+7.D0*SQ6)/3.D0
      DD2=(-13.D0+7.D0*SQ6)/3.D0
      DD3=-1.D0/3.D0
      U1=(6.D0+81.D0**(1.D0/3.D0)-9.D0**(1.D0/3.D0))/30.D0
      U1=ONE/U1
      A11=( 88.D0-  7.D0*SQ6)/360.D0
      A12=(296.D0-169.D0*SQ6)/ 18.D2
      A13=( -2.D0+  3.D0*SQ6)/225.D0
      A21=(296.D0+169.D0*SQ6)/ 18.D2
      A22=( 88.D0+  7.D0*SQ6)/360.D0
      A23=( -2.D0-  3.D0*SQ6)/225.D0
      A31=( 16.D0-       SQ6)/ 36.D0
      A32=( 16.D0+       SQ6)/ 36.D0
      A33= 1.D0/9.D0
C ***
      IDID=0
      REJECT=.FALSE.
      IF (JSTART.GT.0 .AND. MCONST .AND. IU.EQ.0 .AND. H.EQ.HOLD) THEN
         FACTOR=.FALSE.
      ELSE
         FACTOR=.TRUE.
      END IF
      N2=2*N
      HHFAC=H
      IF (.NOT. MCONST .OR. JSTART .EQ. 0) THEN
         DO 2 I=1,N
            DO 1 J=1,N
               E0(I,J)=E(I,J)
               A0(I,J)=A(I,J)
   1        CONTINUE
   2     CONTINUE
      END IF
C
C *** *** *** *** *** *** ***
C  BASIC INTEGRATION STEP
C *** *** *** *** *** *** ***
  10  CONTINUE
      IF (ABS(H).LE.HMIN) GOTO 56
      TPH=T+H
      HA11=H*A11
      HA12=H*A12
      HA13=H*A13
      HA21=H*A21
      HA22=H*A22
      HA23=H*A23
      HA31=H*A31
      HA32=H*A32
      HA33=H*A33
C
C --- Evaluate E,A and F with respect to the three RK stages
C ---   and compute solution matrix and right-hand side
C
C --- First Runge-Kutta stage
      CALL DNFIX1(EDIF, ADIF, FDIF, N, T+C1*H, LMAX, M, ID, IA, IU,
     $            IREQ, MQ, IDQ, IAQ, IUQ, E, LDE, A, LDA, F,
     $            EQ, LDEQ, AQ, LDAQ, FQ, Z, LDZ, ZQ, LDZQ,
     $            AH, LDAH, IPAR, RPAR, WORK, LWORK, MCONST, IER)
      IF (IER .NE. 0) GOTO 60
      IWORK(ONEV)=IWORK(ONEV)+1
      CALL DGEMV('N',N,N,ONE,A,LDA,X,1,ONE,F,1)
      CALL DCOPY(N,F,1,B(1),1)
      IF ((.NOT. MCONST) .AND. FACTOR) THEN
        DO 12 I=1,N
          DO 11 J=1,N
            W(I,J   )=E(I,J) -HA11*A(I,J)
            W(I,N+J )=       -HA12*A(I,J)
            W(I,N2+J)=       -HA13*A(I,J)
  11      CONTINUE
  12    CONTINUE
      END IF
C
C --- Second Runge-Kutta stage
      CALL DNFIX1(EDIF, ADIF, FDIF, N, T+C2*H, LMAX, M, ID, IA, IU,
     $            IREQ, MQ, IDQ, IAQ, IUQ, E, LDE, A, LDA, F,
     $            EQ, LDEQ, AQ, LDAQ, FQ, Z, LDZ, ZQ, LDZQ,
     $            AH, LDAH, IPAR, RPAR, WORK, LWORK, MCONST, IER)
      IF (IER .NE. 0) GOTO 60
      IWORK(ONEV)=IWORK(ONEV)+1
      CALL DGEMV('N',N,N,ONE,A,LDA,X,1,ONE,F,1)
      CALL DCOPY(N,F,1,B(N+1),1)
      IF ((.NOT. MCONST) .AND. FACTOR) THEN
        DO 14 I=1,N
          DO 13 J=1,N
            W(N+I,J   )=       -HA21*A(I,J)
            W(N+I,N+J )=E(I,J) -HA22*A(I,J)
            W(N+I,N2+J)=       -HA23*A(I,J)
  13      CONTINUE
  14    CONTINUE
      END IF
C
C --- Third Runge-Kutta stage
      CALL DNFIX1(EDIF, ADIF, FDIF, N, TPH, LMAX, M, ID, IA, IU,
     $            IREQ, MQ, IDQ, IAQ, IUQ, E, LDE, A, LDA, F,
     $            EQ, LDEQ, AQ, LDAQ, FQ, Z, LDZ, ZQ, LDZQ,
     $            AH, LDAH, IPAR, RPAR, WORK, LWORK, MCONST, IER)
      IF (IER .NE. 0) GOTO 60
      IWORK(ONEV)=IWORK(ONEV)+1
      CALL DGEMV('N',N,N,ONE,A,LDA,X,1,ONE,F,1)
      CALL DCOPY(N,F,1,B(N2+1),1)
      IF ((.NOT. MCONST) .AND. FACTOR) THEN
        DO 16 I=1,N
          DO 15 J=1,N
            W(N2+I,J   )=       -HA31*A(I,J)
            W(N2+I,N+J )=       -HA32*A(I,J)
            W(N2+I,N2+J)=E(I,J) -HA33*A(I,J)
  15      CONTINUE
  16    CONTINUE
      END IF
C
C --- Decomposition and solution of the Runge-Kutta system
      FAC1=U1/H
      IF (MCONST) THEN
         CALL DRKCXS(N,N2,N3,IU,FAC1,H,UROUND,AH,LDAH,W,N,
     $               E0,LDE0,A0,LDA0,B,Y,IWORK,WORK,FACTOR,IER)
      ELSE
         CALL DRKRLS(N,N3,IU,FAC1,UROUND,AH,LDAH,W,LDW,
     $               E0,LDE0,A0,LDA0,B,Y,IWORK,WORK,IER)
      END IF
      IF (IER.ne.0) GOTO 50
C
C --- Compute the vectors Z1, Z2, Z3
      CALL DSCAL(N,ZERO,Z1,1)
      CALL DAXPY(N,HA11,Y      ,1,Z1,1)
      CALL DAXPY(N,HA12,Y(N+1) ,1,Z1,1)
      CALL DAXPY(N,HA13,Y(N2+1),1,Z1,1)
      CALL DSCAL(N,ZERO,Z2,1)
      CALL DAXPY(N,HA21,Y      ,1,Z2,1)
      CALL DAXPY(N,HA22,Y(N+1) ,1,Z2,1)
      CALL DAXPY(N,HA23,Y(N2+1),1,Z2,1)
      CALL DSCAL(N,ZERO,Z3,1)
      CALL DAXPY(N,HA31,Y      ,1,Z3,1)
      CALL DAXPY(N,HA32,Y(N+1) ,1,Z3,1)
      CALL DAXPY(N,HA33,Y(N2+1),1,Z3,1)
C
C *** *** *** *** *** *** ***
C ERROR ESTIMATION
C *** *** *** *** *** *** ***
      HEE1=DD1/H
      HEE2=DD2/H
      HEE3=DD3/H
      CALL DCOPY(N,XPRIME,1,F,1)
      CALL DAXPY(N,HEE1,Z1,1,F,1)
      CALL DAXPY(N,HEE2,Z2,1,F,1)
      CALL DAXPY(N,HEE3,Z3,1,F,1)
C
C --- Consider enlarged system
      CALL DGEMV('N',N,N,ONE,E0,LDE0,F,1,ZERO,CONT(1,1),1)
      CALL DGEMV('N',N,N,ONE,A0,LDA0,F,1,ZERO,CONT(1,2),1)
      IF (IU .EQ. 0) THEN
         CALL DGETRS ('No transpose',N,2,AH,LDAH,IWORK(OIPVT+N3),
     $                CONT,LDCONT,IER)
         IF (IER .NE. 0) GOTO 56
      ELSE
         DO 22 I=1,N
            DO 21 J=1,N
               AH(I,J) = FAC1*E0(I,J) - A0(I,J)
  21        CONTINUE
            IWORK(OIPVT+N3+I-1)=0
  22     CONTINUE
         CALL DGELSX(N,N,2,AH,LDAH,CONT,LDCONT,IWORK(OIPVT+N3),UROUND,
     $               IRANK,WORK,IER)
         IF (IER .NE. 0) GOTO 56
      END IF
C
C --- The error norm (Index-2-components must be multiplied by H)
      ERR=ZERO
      DO 24 I=1,N
         ERR=ERR + (CONT(I,1)*WT(I))**2
     $           + (CONT(I,2)*HHFAC*WT(I))**2
  24  CONTINUE
      ERR=MAX(SQRT(ERR/N2),1.D-10)
      IF (ERR .GE. ONE .AND. (JSTART .EQ. 0 .OR. REJECT)) THEN
C
C ---    Use improved error estimator
         CALL DAXPY(N,ONE,CONT(1,2),1,F,1)
         CALL DGEMV('N',N,N,ONE,E0,LDE0,F,1,ZERO,CONT(1,1),1)
         CALL DGEMV('N',N,N,ONE,A0,LDA0,F,1,ZERO,CONT(1,2),1)
         IF (IU .EQ. 0) THEN
            CALL DGETRS ('No transpose',N,2,AH,LDAH,IWORK(OIPVT+N3),
     $                   CONT,LDCONT,IER)
            IF (IER .NE. 0) GOTO 56
         ELSE
            DO 26 I=1,N
               DO 25 J=1,N
                  AH(I,J) = FAC1*E0(I,J) - A0(I,J)
  25           CONTINUE
               IWORK(OIPVT+N3+I-1)=0
  26        CONTINUE
            CALL DGELSX(N,N,2,AH,LDAH,CONT,N,IWORK(OIPVT+N3),UROUND,
     $                  IRANK,WORK,IER)
            IF (IER .NE. 0) GOTO 56
         END IF
C
C ---    The error norm (Index-2-components must be multiplied by H)
         ERR=ZERO
         DO 28 I=1,N
            ERR=ERR + (CONT(I,1)*WT(I))**2
     $              + (CONT(I,2)*HHFAC*WT(I))**2
  28     CONTINUE
         ERR=MAX(SQRT(ERR/N2),1.D-10)
      END IF
C
C --- Computation of HNEW (we require 0.2 <= HNEW/H <= 8.0)
      FAC=SAFE
      QUOT=MAX(FACR,MIN(FACL,ERR**.25D0/FAC))
      HNEW=H/QUOT
C
C --- Is the error small enough ?
      IF (ERR.LT.ONE) THEN
C *** *** *** *** *** *** ***
C STEP IS ACCEPTED
C *** *** *** *** *** *** ***
         IF (JSTART .EQ. 0) JSTART=1
         IWORK(ONST)=IWORK(ONST)+1
         IF (PRED) THEN
C
C ---    Predictive controller of Gustafsson
            IF (IWORK(ONST).GT.1) THEN
               FACGUS=(HOLD/H)*(ERR**2/ERRACC)**0.25D0/SAFE
               FACGUS=MAX(FACR,MIN(FACL,FACGUS))
               QUOT=MAX(QUOT,FACGUS)
               HNEW=H/QUOT
            END IF
            ERRACC=MAX(1.0D-2,ERR)
         END IF
         HOLD=H
C
C ---    Update T, solution vector X and derivative vector XPRIME
         T=TPH
         CALL DAXPY(N,ONE,Z3,1,X,1)
         CALL DCOPY(N,Y(N2+1),1,XPRIME,1)
         IF (LAST) RETURN
         IF (REJECT) THEN
            POSNEG=SIGN(ONE,H)
            HNEW=POSNEG*MIN(ABS(HNEW),ABS(H))
         END IF
         QT=HNEW/H
         IF (.NOT.MCONST .OR. IU.NE.0 .OR. QT.LT.QUOT1 .OR.
     $          QT.GT.QUOT2) H=HNEW
         RETURN
      ELSE
C *** *** *** *** *** *** ***
C STEP IS REJECTED
C *** *** *** *** *** *** ***
         REJECT=.TRUE.
         LAST=.FALSE.
         FACTOR=.TRUE.
         IF (JSTART .EQ. 0) THEN
             H=H*0.1D0
             HHFAC=0.1D0
         ELSE
             HHFAC=HNEW/H
             H=HNEW
         END IF
         IWORK(OETF)=IWORK(OETF)+1
         GOTO 10
      END IF
C
C --- Unexpected step rejection
  50  CONTINUE
      IF (IER.NE.0) THEN
          NSING=NSING+1
          IF (NSING.GE.5) GOTO 58
      END IF
      H=H*0.5D0
      HHFAC=0.5D0
      REJECT=.TRUE.
      LAST=.FALSE.
      FACTOR=.TRUE.
      GOTO 10
C *** *** *** *** *** *** ***
C FAIL EXIT
C *** *** *** *** *** *** ***
  56  CONTINUE
      IDID=-6
      RETURN
  58  CONTINUE
      IDID=-8
      RETURN
  60  CONTINUE
      IDID=-10
      RETURN
C
C --- Last line of subroutine DRKLIN
      END
