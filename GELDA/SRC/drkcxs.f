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
C====>    version       = "1.00",
C====>    date          = "12 October 1995",
C====>    time          = "13:29:32 MEZ",
C====>    filename      = "drkcxs.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "59500 195 640 6867",
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
      SUBROUTINE DRKCXS(N, N2, N3, IU, FAC1, H, UROUND,
     $                  AH, LDAH, W, LDW, E0, LDE0, A0, LDA0, B, YY,
     $                  IWORK, WORK, FACTOR, IDID)
C
C     PURPOSE
C
C     Reduce the system to be solved in an (n*n) real and a (2n*2n)
C     conjugate complex one [1].
C     Solve these systems and compose solution of the received vectors.
C

C     REFERENCES
C
C     [1] P. Kunkel, V. Mehrmann, W. Rath and J. Weickert.
C         GELDA: A software package for the solution of general linear
C         differential algebraic equations.
C         Preprint SPC 95\_8, TU Chemnitz-Zwickau, February 1995.
C
      IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      INTEGER          OIPVT, ONFA
      PARAMETER (ONFA   = 13)
      PARAMETER (OIPVT  = 21)
C     .. Scalar Arguments ..
      DOUBLE PRECISION FAC1, H, UROUND
      INTEGER          IDID, IU, LDAH, LDA0, LDE0, LDW, N, N2, N3
      LOGICAL          FACTOR
C     .. Array Arguments ..
      DOUBLE PRECISION AH(LDAH,*), A0(LDA0,*), E0(LDE0,*),
     $                 B(*), WORK(*), YY(*)
      COMPLEX*16       W(LDW,*)
      INTEGER          IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION ALPH, ALPHN, BETA, BETAN, CNO,
     $                 T11, T12, T13, T21, T22, T23, T31, TI11, TI12,
     $                 TI13, TI21, TI22, TI23, TI31, TI32, TI33
      INTEGER          I, IER, IRANK, J
C     .. External Subroutines ..
      EXTERNAL         DAXPY, DCOPY, DGELSX, DGETRF, DGETRS, DSCAL,
     $                 ZGELSX, ZGETRF, ZGETRS
C     .. Executable Statements ..
C
      ALPH=(12.D0-81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))/60.D0
      BETA=(81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))*DSQRT(3.D0)/60.D0
      CNO=ALPH**2+BETA**2
      ALPH=ALPH/CNO
      BETA=BETA/CNO
      T11=9.1232394870892942792D-02
      T12=-0.14125529502095420843D0
      T13=-3.0029194105147424492D-02
      T21=0.24171793270710701896D0
      T22=0.20412935229379993199D0
      T23=0.38294211275726193779D0
      T31=0.96604818261509293619D0
      TI11=4.3255798900631553510D0
      TI12=0.33919925181580986954D0
      TI13=0.54177053993587487119D0
      TI21=-4.1787185915519047273D0
      TI22=-0.32768282076106238708D0
      TI23=0.47662355450055045196D0
      TI31=-0.50287263494578687595D0
      TI32=2.5719269498556054292D0
      TI33=-0.59603920482822492497D0
C
      ALPHN=ALPH/H
      BETAN=BETA/H
C
C --- Multiply right hand side by the inverse of T
      CALL DSCAL(N,ZERO,YY(1),1)
      CALL DAXPY(N,TI11,B(1),1,YY(1),1)
      CALL DAXPY(N,TI12,B(N+1),1,YY(1),1)
      CALL DAXPY(N,TI13,B(N2+1),1,YY(1),1)
      CALL DSCAL(N,ZERO,YY(N+1),1)
      CALL DAXPY(N,TI21,B(1),1,YY(N+1),1)
      CALL DAXPY(N,TI22,B(N+1),1,YY(N+1),1)
      CALL DAXPY(N,TI23,B(N2+1),1,YY(N+1),1)
      CALL DSCAL(N,ZERO,YY(N2+1),1)
      CALL DAXPY(N,TI31,B(1),1,YY(N2+1),1)
      CALL DAXPY(N,TI32,B(N+1),1,YY(N2+1),1)
      CALL DAXPY(N,TI33,B(N2+1),1,YY(N2+1),1)
C
C --- Multiply right hand side by the block diagonal matrix Lambda
      CALL DCOPY(N,YY(1),1,B(1),1)
      CALL DSCAL(N,FAC1,B(1),1)
      CALL DSCAL(N,ZERO,B(N+1),1)
      CALL DAXPY(N,ALPHN,YY(N+1),1,B(N+1),1)
      CALL DAXPY(N,-BETAN,YY(N2+1),1,B(N+1),1)
      CALL DSCAL(N,ZERO,YY(1),1)
      CALL DAXPY(N,BETAN,YY(N+1),1,YY(1),1)
      CALL DAXPY(N,ALPHN,YY(N2+1),1,YY(1),1)
C
      DO 11 I=N,1,-1
         B(N+2*I-1)=B(N+I)
         B(N+2*I)=YY(I)
  11  CONTINUE
C
C --- Compute the system matrices
      IF (FACTOR) THEN
         DO 13 I=1,N
            DO 12 J=1,N
               AH(I,J) = FAC1*E0(I,J) - A0(I,J)
               W (I,J) = DCMPLX(ALPHN,BETAN)*E0(I,J) - A0(I,J)
  12        CONTINUE
  13     CONTINUE
      END IF
C
C --- Decomposition and solution
      IF (IU .EQ. 0) THEN
         IF (FACTOR) THEN
            CALL DGETRF (N,N,AH,LDAH,IWORK(OIPVT+N3),IER)
            IF (IER.NE.0) GOTO 98
            CALL ZGETRF (N,N,W,LDW,IWORK(OIPVT),IER)
            IF (IER.NE.0) GOTO 98
            IWORK(ONFA)=IWORK(ONFA)+1
         END IF
         CALL DGETRS ('N',N,1,AH,LDAH,IWORK(OIPVT+N3),B(1),N,IER)
         IF (IER .NE. 0) GOTO 99
         CALL ZGETRS ('N',N,1,W,LDW,IWORK(OIPVT),B(N+1),N,IER)
         IF (IER .NE. 0) GOTO 99
      ELSE
         DO 17 I=1,N
            IWORK(OIPVT+I-1)=0
            IWORK(OIPVT+N3+I-1)=0
  17     CONTINUE
         CALL DGELSX(N,N,1,AH,LDAH,B(1),N,IWORK(OIPVT+N3),UROUND,IRANK,
     $               WORK,IER)
         IF (IER .NE. 0) GOTO 99
         CALL ZGELSX(N,N,1,W,LDW,B(N+1),N,IWORK(OIPVT),UROUND,IRANK,
     $               WORK(1),WORK(6*N+1),IER)
         IF (IER .NE. 0) GOTO 99
         IWORK(ONFA)=IWORK(ONFA)+1
      END IF
C
      DO 19 I=1,N
         YY(I) =B(N+2*I)
         B(N+I)=B(N+2*I-1)
  19  CONTINUE
      CALL DCOPY(N,YY(1),1,B(N2+1),1)
C
C --- Compute y_j (solutions of RK system) from v_j (stored in B), j=1,2,3
      CALL DSCAL(N,ZERO,YY(1),1)
      CALL DAXPY(N,T11,B(1),1,YY(1),1)
      CALL DAXPY(N,T12,B(N+1),1,YY(1),1)
      CALL DAXPY(N,T13,B(N2+1),1,YY(1),1)
      CALL DSCAL(N,ZERO,YY(N+1),1)
      CALL DAXPY(N,T21,B(1),1,YY(N+1),1)
      CALL DAXPY(N,T22,B(N+1),1,YY(N+1),1)
      CALL DAXPY(N,T23,B(N2+1),1,YY(N+1),1)
      CALL DSCAL(N,ZERO,YY(N2+1),1)
      CALL DAXPY(N,T31,B(1),1,YY(N2+1),1)
      CALL DAXPY(N,ONE,B(N+1),1,YY(N2+1),1)
C
      RETURN
C
C --- Exit in decomposition routine
  98  CONTINUE
      IDID=2
      RETURN
C
C --- Exit in solution routine
  99  CONTINUE
      IDID=1
      RETURN
C
C --- Last line of subroutine DRKCXS
      END
