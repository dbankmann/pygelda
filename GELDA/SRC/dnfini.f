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
C====>    time          = "10:25:46 MESZ",
C====>    filename      = "dnfini.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "43154 202 916 7144",
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
      SUBROUTINE DNFINI(N,ID,IA,IU,X,E,LDE,A,LDA,F,DELTA,TAU,
     $                  RWORK,LRWORK,METHOD,IERR)
C
C     PURPOSE
C
C     DNFRED is part of the code DGELDA [2] which solves linear DAEs
C     with variable coefficients of the form
C
C                               .
C                           E(t)x(t) = A(t)x(t) + f(t)
C                             x(t_0) = x_0
C
C     DNFINI computes consistent initial condition x_0.
C
C     ARGUMENT LIST
C
C        ARGUMENTS IN
C
C         N - INTEGER.
C             The number of equations in the DAE system.
C             N .GE. 1.
C
C        ID - INTEGER.
C             The number DMU of differential components at the previous
C             step.
C
C        IA - INTEGER.
C             The number AMU of algebraic components at the previous
C             step.
C
C        IU - INTEGER.
C             The number UAU of undetermined components at the previous
C             step.
C
C         X - DOUBLE PRECISION array of DIMENSION (*).
C             An guess for a consistent x_0.
C
C        E - DOUBLE PRECISION array of DIMENSION (LDE,*).
C            The leading N by N part of this array must contain the
C            matrix E of the strangeness free DAE at a fixed time T.
C
C      LDE - INTEGER.
C            The leading dimension of array E as declared in the
C            calling program.
C            LDE .GE. N.
C
C        A - DOUBLE PRECISION array of DIMENSION (LDA,*).
C            The leading N by N part of this array must contain the
C            matrix A of the strangeness free DAE at a fixed time T.
C
C      LDA - INTEGER.
C            The leading dimension of array A as declared in the
C            calling program.
C            LDA .GE. N.
C
C        F - DOUBLE PRECISION array of DIMENSION at least (*).
C            The leading N elements of this array must contain the
C            vector F of the strangeness free DAE at a fixed time T.
C
C     ARGUMENTS OUT
C
C        X - DOUBLE PRECISION array of DIMENSION (*).
C            The computed consistent initial value.
C
C    DELTA - DOUBLE PRECISION array of DIMENSION at least (N).
C            This array is used as workspace.
C
C      TAU - DOUBLE PRECISION array of DIMENSION at least (N).
C            This array is used as workspace.
C
C     WORK SPACE
C
C     RWORK - DOUBLE PRECISION array of DIMENSION at least (LRWORK).
C
C    LRWORK - INTEGER.
C             The dimension of RWORK.
C             LRWORK .GE. IA + N.
C
C             NOTE that for optimum performance LWORK .GE. IA + N*NB,
C             where NB is the optimum block size of LAPACK.
C
C     MODE PARAMETER
C
C     METHOD - INTEGER.
C              Indicates which method the user wishes to use for the
C              computation of consistant initial values as follows:
C
C              METHOD = 0, (All types of variables in the guess for
C                           consistent initial can be changed);
C
C              METHOD = 1, (The differential variables in the guess for
C                           consistant initial values should be fixed).
C
C     ERROR IDICATOR
C
C       IERR - INTEGER.
C              Unless the routine detects an error (see next section),
C              IERR contains 0 on exit.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C     IERR = -5 : LRWORK < IA + N.
C
C     REVISIONS
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
C     1995, July 17 [Version 1.1]
C       Changed documentation to meet SLICOT standard.
C
C     1995, July 10 [Version 1.0]
C       First release. (W. Rath, J. Weickert)
C
C     ************************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER        (ONE = 1.0D0)
C     .. Scalar Arguments ..
      INTEGER          N, ID, IA, IU, LDE, LDA, LRWORK, IERR, METHOD
C     .. Array Arguments ..
      DOUBLE PRECISION X(*), E(LDE,*), A(LDA,*), F(*),
     $                 RWORK(*), DELTA(*), TAU(*)
C     .. External Subroutines ..
      EXTERNAL         DAXPY, DCOPY, DGELQF, DGELS, DGEMV, DORMLQ
C     .. Executable Statements ..
      IF (IA.EQ.0) RETURN
      IF (LRWORK-IA-N.LT.0) THEN
        IERR = -5
	RETURN
      ENDIF
      IF (METHOD.EQ.1) THEN
         CALL DGELQF(ID,N,E,LDE,TAU,RWORK,LRWORK,IERR)
         CALL DORMLQ('R','T',N,N,ID,E,LDE,TAU,A,LDA,RWORK,LRWORK,
     $               IERR)
         CALL DORMLQ('L','N',N,1,ID,E,LDE,TAU,X,N,RWORK,LRWORK,
     $               IERR)
         CALL DCOPY(IA,F(ID+1),1,DELTA,1)
         CALL DGEMV('N',IA,N,ONE,A(ID+1,1),LDA,X,1,ONE,DELTA,1)
C
C     Compute the minimum norm solution. Note IA <= IA + IU.
         CALL DGELS('N',IA,IA+IU,1,A(ID+1,ID+1),LDA,DELTA,N,
     $              RWORK,LRWORK,IERR)
         CALL DAXPY(IA+IU,-ONE,DELTA,1,X(ID+1),1)
         CALL DORMLQ('L','T',N,1,ID,E,LDE,TAU,X,N,RWORK,LRWORK,
     $               IERR)
      ELSE
         CALL DCOPY(IA,F(ID+1),1,DELTA,1)
         CALL DGEMV('N',IA,N,ONE,A(ID+1,1),LDA,X,1,ONE,DELTA,1)
C
C     Compute the minimum norm solution. Note IA <= N.
         CALL DGELS('N',IA,N,1,A(ID+1,1),LDA,DELTA,N,RWORK,LRWORK,IERR)
         CALL DAXPY(N,-ONE,DELTA,1,X,1)
      ENDIF
      IERR=0
      RETURN
C *** Last line of DNFINI ***
      END
