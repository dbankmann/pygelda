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
C====>    version       = "1.01",
C====>    date          = "18 July 1995",
C====>    time          = "15:22:45 MESZ",
C====>    filename      = "dnfotr.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "49975 237 1113 8238",
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
      SUBROUTINE DNFOTR(N, NQ, ID, IA, IREQ, EQ, LDEQ, AQ, LDAQ, Z2Q,
     $                  LDZ2Q, Z1, LDZ1, AH, LDAH, RWORK, LRWORK, IERR)
C
C     PURPOSE
C
C     DNFOTR is part of the code DGELDA [2] which solves linear DAEs
C     with variable coefficients of the form
C
C                               .
C                           E(t)x(t) = A(t)x(t) + f(t)
C                             x(t_0) = x_0
C
C     DNFOTR computes the (local) orthogonal transoformations which
C     extract a strangeness free DAE from the extended derivative array
C     at a fixed time.
C
C     ARGUMENT LIST
C
C         N - INTEGER.
C             The number of equations in the DAE system.
C             N .GE. 1.
C
C        NQ - INTEGER.
C             The size of EQ and AQ, i.e. NQ=(MU+1)*N, where MU is the
C             strangeness index of the DAE.
C             NQ .GE. 1.
C
C        ID - INTEGER.
C             The number DMU of differential components.
C
C        IA - INTEGER.
C             The number AMU of algebraic components.
C
C        IU - INTEGER.
C             The number UAU of undetermined components.
C
C      IREQ - INTEGER.
C             The rank of the matrix EQ.
C
C        EQ - DOUBLE PRECISION array of DIMENSION (LDEQ,*).
C             The leading NQ by NQ part of this array contains the
C             matrix EQ.
C
C      LDEQ - INTEGER.
C             The leading dimension of array EQ as declared in the
C             calling program.
C             LDEQ .GE. NQ.
C
C        AQ - DOUBLE PRECISION array of DIMENSION (LDAQ,*).
C             The leading NQ by NQ part of this array contains the
C             matrix AQ.
C
C      LDAQ - INTEGER.
C             The leading dimension of array AQ as declared in the
C             calling program.
C             LDAQ .GE. NQ.
C
C       Z2Q - DOUBLE PRECISION array of DIMENSION (LDZ2Q,*)
C             The leading NQ by NQ part of this array must contain the
C             right orthogonal factor of the SVD EQ=ZQ*S*VT.
C             NOTE that this array is overwritten.
C
C     LDZ2Q - INTEGER.
C             The leading dimension of array Z2Q as declared in the
C             calling program.
C             LDZ2Q .GE. NQ.
C
C      LDZ1 - INTEGER.
C             The leading dimension of array Z1 as declared in the
C             calling program.
C             LDZ1 .GE. N.
C
C      LDAH - INTEGER.
C             The leading dimension of array AH as declared in the
C             calling program.
C             LDZ1 .GE. NQ.
C
C       ARGUMENTS OUT
C
C       Z2Q - DOUBLE PRECISION array of DIMENSION (LDZ2Q,*).
C             The leading NQ by IA part of this array contains the
C             transformation matrix Z2Q which extracts the algebraic
C             part.
C
C        Z1 - DOUBLE PRECISION array of DIMENSION (LDZ1,*).
C             The leading N by ID part of this array contains the
C             transformation matrix Z1 which extracts the differential
C             part.
C
C        AH - DOUBLE PRECISION array of DIMENSION (LDAH,*).
C             This array is used as workarray.
C
C     WORK SPACE
C
C      RWORK - DOUBLE PRECISION array of DIMENSION at least (LRWORK).
C
C     LRWORK - NTEGER.
C              The length of RWORK.
C              LRWORK .GE. 3*N + 2*NQ + 2*NQ*NQ.
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
C     IERR = -5 : LRWORK < 3*N + 2*NQ + 2*NQ*NQ.
C
C     METHOD
C
C     All information we need to transform the original DAE is hidden in
C     the extended derivative array EQ, AQ anf FQ. The two orthogonal
C     transformation matricies Z1 and Z2Q are computed via two singular
C     value decompositions (SVDs). The LAPACK subroutine DGESVD is used
C     to compute these SVDs.
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
C     1995, July 18 [Version 1.1]
C       Changed documentation to meet SLICOT standard.
C
C     1995, July 10 [Version 1.0]
C       First release. (W. Rath, J. Weickert)
C
C     ******************************************************************
C
      IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER        (ZERO = 0.0D0, ONE = 1.0D0)
      INTEGER          OSIGMA
      PARAMETER        (OSIGMA = 1)
C     .. Scalar Arguments ..
      INTEGER          IA, ID, IERR, IREQ, LDEQ, LDAH, LDAQ, LDZ1,
     $                 LDZ2Q, LRWORK, N, NQ
C     .. Array Arguments ..
      DOUBLE PRECISION EQ(LDEQ,*), AH(LDAH,*), AQ(LDAQ,*), RWORK(*),
     $                 Z1(LDZ1,*), Z2Q(LDZ2Q,*)
C     .. Local Scalars
      INTEGER          INFO, LWORK, OU, OVT, OWORK
C     .. External Subroutines ..
      EXTERNAL         DGEMM, DGESVD, DLACPY, DLASET
C     .. Executable Statements ..
C
C     SET UP RWORK STORAGE.
      OU = OSIGMA + NQ
      OVT = OU + NQ*NQ
      OWORK = OVT + NQ*NQ
      LWORK = LRWORK - 3*N - 2*NQ -2*NQ*NQ
      IF (LWORK.LT.0) THEN
        IERR = 10
	RETURN
      ENDIF
C
C     NOTE: IREQ .EQ. 0, iff the DAE is an ODE, i.e. M=0 and NQ .EQ. N
      IF (IREQ.EQ.NQ) THEN
         CALL DLASET('N',NQ,NQ,ZERO,ONE,Z2Q,LDZ2Q)
         CALL DLASET('N',N,N,ZERO,ONE,Z1,LDZ1)
         RETURN
      END IF
C
C     Compute AH = Z2Q**T * AQ.
      CALL DGEMM('T','N',NQ-IREQ,N,NQ,ONE,Z2Q,LDZ2Q,AQ,LDAQ,
     $     ZERO,AH,LDAH)
      CALL DGESVD('A','A',NQ-IREQ,N,AH,LDAH,RWORK(OSIGMA),
     $     RWORK(OU),NQ,RWORK(OVT),NQ,RWORK(OWORK),LWORK,INFO)
C
C     Compute AH = Z2Q * U.
      CALL DGEMM('N','N',NQ,NQ-IREQ,NQ-IREQ,ONE,Z2Q,LDZ2Q,
     $     RWORK(OU),NQ,ZERO,AH,LDAH)
C
C     Set Z2Q = AH.
      CALL DLACPY('F',NQ,NQ-IREQ,AH,LDAH,Z2Q,LDZ2Q)
      IF (ID.NE.0) THEN
C
C     Compute Z1 = EQ * VT**T. NOTE: Z1 is NxN-IA <= NxN.
         CALL DGEMM('N','T',N,N-IA,N,ONE,EQ,LDEQ,RWORK(OVT+IA),NQ,ZERO,
     $              Z1,LDZ1)
         CALL DGESVD('O','N',N,N-IA,Z1,LDZ1,RWORK(OSIGMA),RWORK(OU),NQ,
     $               RWORK(OVT),NQ,RWORK(OWORK),LWORK,INFO)
      END IF
      RETURN
C *** Last line of DNFOTR ***
      END
