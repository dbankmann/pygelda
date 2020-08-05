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
C====>    date          = "17 July 1995",
C====>    time          = "19:11:16 MESZ",
C====>    filename      = "dnfpen.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "34803 275 1452 10296",
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
      SUBROUTINE DNFPEN(EDIF, ADIF, FDIF, N, T, L, WQ, LDWQ,
     $                  EQ, LDEQ, AQ, LDAQ, FQ, IPAR, RPAR, IERR)
C
C     PURPOSE
C
C     DNFPEN is part of the code DGELDA [2] which solves linear DAEs
C     with variable coefficients of the form
C
C                               .
C                           E(t)x(t) = A(t)x(t) + f(t)
C                             x(t_0) = x_0
C
C     DNFPEN forms the extended derivative array with the matricies EQ,
C     AQ and the vector FQ.
C
C     ARGUMENT LIST
C
C       USER-SUPPLIED SUBROUTINES
C
C      EDIF - User supplied SUBROUTINE.
C             This is a subroutine which the user provides to define the
C             matrix E(t) and its derivatives. It has the form
C
C             SUBROUTINE EDIF(N,T,IDIF,E,LDE,IPAR,RPAR,IERR).
C
C             The subroutine takes as input the number of equations N,
C             the time T and the integer parameter IDIF. EDIF must not
C             alter these input parameters or the leading dimension LDE
C             of the array E. As output, the subroutine produces the
C             IDIF-th derivative of E(t) at time T in the leading N by N
C             part of the array E. The integer flag IERR is always zero
C             on input and EDIF should alter IERR only if IDIF is larger
C             than the highest derivative of E(t) the subroutine
C             provides (set IERR=-2) or if another problem occurs (set
C             IERR=-1). IPAR and RPAR are integer and real arrays which
C             can be used for the communication between the user's
C             calling program and the subroutines EDIF, ADIF and FDIF.
C
C             In the calling program, EDIF must be declared as external.
C
C      ADIF - User supplied SUBROUTINE.
C             This is a subroutine which the user provides to define the
C             matrix A(t) and its derivatives. It is of the form
C
C              SUBROUTINE ADIF(N,T,IDIF,A,LDA,IPAR,RPAR,IERR)
C
C             and the input and output parameters are similar to these
C             of EDIF.
C
C      FDIF - User supplied SUBROUTINE.
C             This is a subroutine which the user provides to define the
C             vector f(t) and its derivatives. It is of the form
C
C             SUBROUTINE FDIF(N,T,IDIF,F,IPAR,RPAR,IERR)
C
C             and the input and output parameters are similar to these
C             of EDIF. Except that the first N elements of the
C             1-dimensional array F contain the IDIF-th derivative of
C             f(t) at time T. Note further, since F is a 1-dimensional
C             array no leading dimension is needed.
C
C       ARGUMENTS IN
C
C         N - INTEGER.
C             The number of equations in the DAE system.
C             N .GE. 1.
C
C         T - DOUBLE PRECISION.
C             The time T.
C
C         L - INTERGER
C             The stage of the extended derivative array to be computed.
C
C        WQ - DOUBLE PRECISION array of DIMENSION (LDWQ,*).
C             The leading N*L by N part of this array must contain
C              (j)
C             E  (T), j=0,...,L-1.
C
C      LDWQ - INTEGER.
C             The leading dimension of array WQ as declared in the
C             calling program.
C             LDWQ .GE. N.
C
C        EQ - DOUBLE PRECISION array of DIMENSION (LDEQ,*).
C             The leading N*L by N*L part of this array must contain the
C             matrix EQ of stage L-1.
C
C      LDEQ - INTEGER.
C             The leading dimension of array EQ as declared in the
C             calling program.
C             LDEQ .GE. N*(L+1).
C
C        AQ - DOUBLE PRECISION array of DIMENSION (LDAQ,*).
C             The leading N*L by N*L part of this array must contain the
C             matrix AQ of stage L-1.
C
C      LDAQ - INTEGER.
C             The leading dimension of array AQ as declared in the
C             calling program.
C             LDAQ .GE. N*(L+1).
C
C        FQ - DOUBLE PRECISION array of DIMENSION (*).
C             The first N*L components of this array must contain the
C             vector FQ of stage L-1.
C
C       ARGUMENTS OUT
C
C        WQ - DOUBLE PRECISION array of DIMENSION (LDWQ,*).
C             The leading N*(L+1) by N part of this array contains
C              (j)
C             E  (T), j=0,...,L.
C
C        EQ - DOUBLE PRECISION array of DIMENSION (LDEQ,*).
C             The leading N*(L+1) by N*(L+1) part of this array contains
C             the matrix EQ of stage L.
C
C        AQ - DOUBLE PRECISION array of DIMENSION (LDAQ,*).
C             The leading N*(L+1) by N*(L+1) part of this array contains
C             the matrix AQ of stage L.
C
C        FQ - DOUBLE PRECISION array of DIMENSION (*).
C             The first N*(L+1) components of this array contain the
C             vector FQ of stage L.
C
C      IPAR - INTEGER array of DIMENSION (*).
C             This integer array can be used for communication between
C             the calling program and the EDIF, ADIF and FDIF
C             subroutines.
C
C      RPAR - DOUBLE PRECISION array of DIMENSION (*).
C             This real array can be used for communication between the
C             calling program and the EDIF, ADIF and FDIF subroutines.
C
C             IPAR and RPAR are not altered by DNFPEN. If these arrays
C             are used, they must be dimensioned in the calling program
C             and in EDIF, ADIF and FDIF as arrays of appropriate
C             length. Otherwise, ignore these arrays by treating them as
C             dummy arrays of length one.
C
C     ERROR INDICATOR
C
C      IERR - INTEGER.
C             Unless the routine detects an error (see next section),
C             IERR contains 0 on exit.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C     IERR = -1 : An error occured in EDIF, ADIF or FDIF.
C     IERR = -2 : The input parameter IDIF, i.e. L, of EDIF, ADIF or
C                 FDIF was larger then the highest derivative the
C                 routine provides.
C
C     METHOD
C
C     If we denote the i-th derivative by a superscript (i) then EQ, AQ
C     and FQ are defined componentwise by
C
C                  /i\   (i-j)      / i \   (i-j-1)
C      (EQ)    := |   | E    (T) - |     | A      (T)  ; i,j = 0,...,L
C          i,j     \j/              \j+1/
C
C                  /   (i)
C                 |   A   (T) ; for i = 0,...,L , j=0
C      (AQ)    := |
C          i,j     \  0       ; otherwise
C
C                  (i)
C      (FQ)    := f  (T)      ; i = 0,...,L
C          i
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
C     FURTHER COMMENTS
C     The subroutines must be used iteratively to obtain EQ, AQ and FQ
C     for a specific stage L.
C
C     CONTRIBUTORS
C
C     W. Rath, J. Weickert (TU Chemnitz)
C
C     REVISIONS
C
C     1995, July 17 [Version 1.1]
C       Changed documentation to meet SLICOT standard.
C
C     1995, July 10 [Version 1.0]
C       First release. (W. Rath, J. Weickert)
C
C     ************************************************************************
C
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO = 0.0D0)
C     .. Subroutine Arguments ..
      EXTERNAL         EDIF, ADIF, FDIF
C     .. Scalar Arguments ..
      DOUBLE PRECISION T
      INTEGER          IERR, L, LDWQ, LDEQ, LDAQ, N
C     .. Array Arguments ..
      DOUBLE PRECISION AQ(LDAQ,*), EQ(LDEQ,*), FQ(*), RPAR(*),
     $                 WQ(LDWQ,*)
      INTEGER          IPAR(*)
C     .. Local Scalars ..
      DOUBLE PRECISION C1, C2
      INTEGER          I, J, JQ, LN
C     .. External Functions ..
      DOUBLE PRECISION DBINOM
      EXTERNAL         DBINOM
C     .. External Subroutines ..
      EXTERNAL         DLACPY, DLASET
C     .. Executable Statements ..
      LN = L*N
C     Set blocks in EQ.
      CALL DLASET('N',LN,N,ZERO,ZERO,EQ(1,LN+1),LDEQ)
C     Compute new derivatives.
      CALL EDIF(N,T,L,WQ(LN+1,1),LDWQ,IPAR,RPAR,IERR)
      IF (IERR.NE.0) RETURN
      CALL ADIF(N,T,L,AQ(LN+1,1),LDAQ,IPAR,RPAR,IERR)
      IF (IERR.NE.0) RETURN
      CALL FDIF(N,T,L,FQ(LN+1),IPAR,RPAR,IERR)
      IF (IERR.NE.0) RETURN
      C1 = 1.0D0
      DO 120 JQ=0,L-1
        C2 = DBINOM(L,JQ+1)
        DO 110 J=1,N
          DO 100 I=1,N
            EQ(LN+I,JQ*N+J)= C1 * WQ((L-JQ)*N+I,J)
     $                     - C2 * AQ((L-JQ-1)*N+I,J)
  100     CONTINUE
  110   CONTINUE
        C1 = C2
  120 CONTINUE
      CALL DLACPY('N',N,N,WQ,LDWQ,EQ(LN+1,LN+1),LDEQ)
      RETURN
C *** Last line of DNFPEN ***
      END
