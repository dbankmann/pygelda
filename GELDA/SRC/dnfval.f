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
C====>    time          = "13:30:53 MEZ",
C====>    filename      = "dnfval.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "62322 296 1162 9354",
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
      SUBROUTINE DNFVAL(N, NQ, EQ, LDEQ, AQ, LDAQ, R0, A0, S0, D0,
     $                  U0, ZQ, LDZQ, RWORK, LRWORK, IERR)
C
C     PURPOSE
C
C     DNFVAL is part of the code DGELDA [2] which solves linear DAEs
C      with variable coefficients of the form
C
C                               .
C                           E(t)x(t) = A(t)x(t) + f(t)
C                             x(t_0) = x_0
C
C     DNFVAL computes the characteristic values of a matrix pair (EQ,AQ)
C     and the right orthogonal factor ZQ of the SVD EQ = ZQ*S*VT.
C
C     ARGUMENT LIST
C
C       ARGUMENTS IN
C
C         N - INTEGER.
C             The number of equations of the original DAE.
C             N .GE. 1.
C
C        NQ - INTEGER.
C             The size of the square matricies EQ and AQ.
C             NQ .GE. N.
C
C        EQ - DOUBLE PRECISION array of DIMENSION (LDEQ,*).
C             The leading NQ by NQ part of this array must contain the
C             matrix EQ.
C
C      LDEQ - INTEGER.
C             The leading dimension of array EQ as declared in the
C             calling program.
C             LDEQ .GE. NQ.
C
C        AQ - DOUBLE PRECISION array of DIMENSION (LDAQ,*).
C             The leading NQ by NQ part of this array must contain the
C             matrix AQ.
C
C      LDAQ - INTEGER.
C             The leading dimension of array AQ as declared in the
C             calling program.
C             LDAQ .GE. NQ.
C
C      LDZQ - INTEGER.
C             The leading dimension of array ZQ as declared in the
C             calling program.
C             LDAQ .GE. NQ.
C
C       ARGUMENTS OUT
C
C        R0 - INTEGER.
C             The rank of EQ.
C
C        A0 - INTEGER.
C             The size of the algebraic part.
C
C        S0 - INTEGER.
C             The size of the strangeness part.
C
C        D0 - INTEGER.
C             The size of the differential part.
C
C        U0 - INTEGER.
C             The size of the undetermind part.
C
C        ZQ - DOUBLE PRECISION array of DIMENSION (LDZQ,*).
C             The leading NQ by NQ part of this array contains the right
C             orthogonal factor of the SVD EQ = ZQ*S*VT.
C
C     WORK SPACE
C
C      RWORK - DOUBLE PRECISION array of DIMENSION at least (LRWORK).
C
C     LRWORK - NTEGER.
C              The length  of RWORK.
C              LRWORK .GE. 3*NQ*NQ + 6*NQ.
C
C              NOTE that for good performance, LRWORK should generally
C              be larger.
C
C     ERROR INDICATOR
C
C       IERR - INTEGER.
C              Unless the routine detects an error (see next section),
C              IERR contains 0 on exit.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C     IERR = -3 : An argument of DGESVD had an illegal value.
C     IERR = -4 : DGESVD failed to converge.
C     IERR = -5 : LRWORK < 3NQ*NQ + 6*NQ.
C
C     METHOD
C
C     The characteristic values of a NQ by NQ matrix pair (EQ,AQ) are
C     defined by
C
C                R0 = rank EQ
C                A0 = rank (Z^* AQ T)
C                S0 = rank(V^* Z^* AQ T')
C                D0 = R0 - S0
C                U0 = NQ - R0 - A0 - S0
C
C     where
C
C                T basis of kernel EQ
C                Z basis of corange EQ
C                T' basis of cokernel EQ
C                V basis of corange(Z^* AQ T).
C
C     The characteristic values R0,A0,S0,D0,U0 are computed via three
C     rank decisions, which are obtained by means of singular value
C     decompositions (SVDs). The LAPACK subroutine DGESVD is used to
C     compute these SVDs.
C
C     REFERENCES
C
C     [1] Peter Kunkel, Volker Mehrmann.
C         A New Class of Discretization Methods for the Solution of
C         Linear Differential-Algebraic Equations with Variable
C         Coefficients.
C         Materialien LXII, FSP Schwerpunkt Mathematisierung,
C         Universitaet Bielefeld, FRG.
C         To appear in SIAM J. Numer. Anal.
C
C     [2] P. Kunkel, V. Mehrmann, W. Rath and J. Weickert.
C         GELDA: A software package for the solution of general linear
C         differential algebraic equations.
C         Preprint SPC 95_8, TU Chemnitz-Zwickau, February 1995.
C
C     CONTRIBUTORS
C
C     W. Rath, J. Weickert (TU Chemnitz, GERMANY)
C
C     REVISIONS
C
C     1995, July 21 [Version 1.1.1]
C       Nonstandard use of debug code left. Removed this code. (W. Rath)
C
C     1995, July 17 [Version 1.1]
C       Changed documentation to meet SLICOT standard.
C
C     1995, July 10 [Version 1.0]
C       First release. (W. Rath, J. Weickert)
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER        (ZERO = 0.0D0, ONE = 1.0D0)
      INTEGER          OWQ
      PARAMETER        (OWQ = 1)
C     .. Scalar Arguments ..
      INTEGER          A0, D0, IERR, LDEQ, LDAQ, LDZQ, LRWORK, N, NQ,
     $                 R0, S0, U0
C     .. Array Arguments ..
      DOUBLE PRECISION AQ(LDAQ,*), EQ(LDEQ,*), RWORK(*), ZQ(LDZQ,*)
C     .. Local Scalars ..
      DOUBLE PRECISION UMACH
      INTEGER          I, LWORK, M, OS, OU, OVT, OWORK
C     .. External Subroutines ..
      EXTERNAL         DGEMM, DGESVD, DLACPY
C     .. Exectutable Statements ..
C
C     Set up rwork storage.
      OS = OWQ + NQ*NQ
      OU = OS + NQ
      OVT = OU + NQ*NQ
      OWORK = OVT + NQ*NQ
      LWORK = LRWORK - 3*NQ*NQ - 6*NQ
      IF (LWORK.LT.0) THEN
        IERR = -5
	RETURN
      ENDIF
      UMACH=1.D-12
C
C     Determine R0.
      CALL DLACPY('N',NQ,NQ,EQ,LDEQ,RWORK(OWQ),NQ)
      CALL DGESVD('A','A',NQ,NQ,RWORK(OWQ),NQ,RWORK(OS),ZQ,LDZQ,
     $            RWORK(OVT),NQ,RWORK(OWORK),LWORK,IERR)
      IF (IERR.NE.0) THEN
         IF (IERR.LT.0) THEN
            IERR = -3
         ELSE
            IERR = -4
         ENDIF
         RETURN
      ENDIF
      M=NQ
      R0=0
      IF (RWORK(OS).GT.UMACH) THEN
         DO 10 I=0,M-1
            IF (RWORK(OS+I).GT.UMACH*RWORK(OS)) THEN
               R0=R0+1
            END IF
 10      CONTINUE
      END IF
C
C     Determine A0.
      A0=NQ-R0
      IF (A0.NE.0) THEN
C
C       WQ(1:A0,1:A0)=U(1:NQ,R0+1:R0+A0)^T AQ(1:NQ,1:N) VT(R0+1:R0+A0,1:N)^T
        CALL DGEMM('T','N',A0,N,NQ,ONE,ZQ(1,R0+1),LDZQ,AQ,LDAQ,
     $             ZERO,RWORK(OU),NQ)
	CALL DGEMM('N','T',A0,A0,N,ONE,RWORK(OU),NQ,RWORK(OVT+R0),NQ,
     $             ZERO,RWORK(OWQ),NQ)
        CALL DGESVD('A','N',A0,A0,RWORK(OWQ),NQ,RWORK(OS),RWORK(OU),NQ,
     $               RWORK(OVT),NQ,RWORK(OWORK),LWORK,IERR)
        IF (IERR.NE.0) THEN
           IF (IERR.LT.0) THEN
              IERR = -3
           ELSE
              IERR = -4
           ENDIF
           RETURN
      ENDIF
        M=A0
        A0=0
        IF (RWORK(OS).GT.UMACH) THEN
          DO 20 I=0,M-1
            IF (RWORK(OS+I).GT.UMACH*RWORK(OS)) THEN
              A0=A0+1
            END IF
 20       CONTINUE
        END IF
      END IF
C
C     Determine S0.
      S0=NQ-R0-A0
      IF (S0.NE.0 .AND. R0.NE.0) THEN
C
C        WQ(1:S0,1:R0)=UQ(1:NQ-R0,A0+1:A0+S0)^T U(1:NQ,R0+1:R0+NQ)^T
C                     *AQ(1:NQ,1:N) VT(1:R0,1:N)^T
         CALL DGEMM('T','T',S0,NQ,NQ-R0,ONE,RWORK(OU+A0*NQ),NQ,
     $              ZQ(1,R0+1),LDZQ,ZERO,RWORK(OWQ),NQ)
         CALL DGEMM('N','N',S0,N,NQ,ONE,RWORK(OWQ),NQ,AQ,LDAQ,
     $              ZERO,RWORK(OU),NQ)
	 CALL DGEMM('N','T',S0,R0,N,ONE,RWORK(OU),NQ,RWORK(OVT),NQ,
     $              ZERO,RWORK(OWQ),NQ)
         CALL DGESVD('N','N',S0,R0,RWORK(OWQ),NQ,RWORK(OS),RWORK(OU),NQ,
     $                RWORK(OVT),NQ,RWORK(OWORK),LWORK,IERR)
         IF (IERR.NE.0) THEN
            IF (IERR.LT.0) THEN
               IERR = -3
            ELSE
               IERR = -4
            ENDIF
            RETURN
         ENDIF
         M=S0
         S0=0
         IF (RWORK(OS).GT.UMACH) THEN
            DO 30 I=0,M-1
               IF (RWORK(OS+I).GT.UMACH*RWORK(OS)) THEN
                  S0=S0+1
               END IF
 30         CONTINUE
         END IF
      END IF
C
C     Finish computation and print results.
      D0=R0-S0
      U0=NQ-R0-A0-S0
      RETURN
C *** Last line of DNFVAL ***
      END
