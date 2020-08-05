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
C====>    date          = "15 July 1995",
C====>    time          = "17:05:06 MESZ",
C====>    filename      = "dbdtrp.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "14283 155 646 5123",
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
      SUBROUTINE DBDTRP (T, TOUT, NEQ, KOLD, PHI, LDPHI, PSI,
     $                   XOUT, XPOUT)
C
C     PURPOSE
C
C     Interpolation routine for DBDF.
C
C     ARGUMENT LIST
C
C       ARGUMENTS IN
C
C         T - DOUBLE PRECISION.
C             The current time in the integration.
C
C      TOUT - DOUBLE PRECISION.
C             The time at which the solution is desired.
C
C       NEQ - INTEGER.
C             The actual dimension of the problem.
C
C      KOLD - INTEGER.
C             Order used on last successful step.
C
C       PHI - DOUBLE PRECISION array of DIMENSION (LDPHI,*).
C             The leading NEQ by KOLD+1 part of this array must
C             contain scaled divided differences of X.
C     LDPHI - INTEGER.
C             The leading dimension of array PHI as declared in the
C             callingprogram.
C             LDPHI .GE. NEQ
C
C       PSI - DOUBLE PRECISION array of DIMENSION (*).
C             This array must contain the past stepsize history.
C
C       ARGUMENTS OUT
C
C      XOUT - DOUBLE PRECISION.
C             This variable contains the interpolated approximation
C             to X at TOUT.
C
C     XPOUT - DOUBLE PRECISION.
C             This variable contains the interpolated approximation
C             to XPRIME at TOUT.
C
C     METHOD
C
C     The methods in subroutine DBDSTP use polynomials to approximate
C     the solution. DBDTRP approximates the solution and its
C     derivative at time TOUT by evaluating one of these polynomials,
C     and its derivative, there. Information defining this polynomial
C     is passed from DBDSTP, so DBDTRP cannot be used alone.
C
C     DBDTRP is a modiefied version of the DBDTRP subroutine of DASSL
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
C         Preprint SPC 95\_8, TU Chemnitz-Zwickau, February 1995.
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
C     1995, July 15 [Version 1.1]
C       Changed documentation to meet SLICOT standard. (W. Rath)
C
C     1995, July 10 [Version 1]
C       First release. (W. Rath)
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION T, TOUT
      INTEGER          KOLD, LDPHI, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION PHI(LDPHI,*), PSI(*), XOUT(*), XPOUT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION C, D, GAMMA, TEMP1
      INTEGER          J, KOLDP1
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY
C     .. Executable Statements ..
      KOLDP1=KOLD+1
      TEMP1=TOUT-T
C
C     Set
C      XOUT(*)=PHI(*,1)
C      XPOUT(*)=0.0D0
      CALL DCOPY(NEQ,PHI,1,XOUT,1)
      DO 10 J=1,NEQ
         XPOUT(J) = ZERO
 10   CONTINUE
      C=ONE
      D=ZERO
      GAMMA=TEMP1/PSI(1)
      DO 20 J=2,KOLDP1
         D=D*GAMMA+C/PSI(J-1)
         C=C*GAMMA
         GAMMA=(TEMP1+PSI(J-1))/PSI(J)
C
C        Compute
C          XOUT(*)=XOUT(*)+C*PHI(*,J)
C          XPOUT(*)=XPOUT(*)+D*PHI(*,J)
         CALL DAXPY(NEQ,C,PHI(1,J),1,XOUT,1)
         CALL DAXPY(NEQ,D,PHI(1,J),1,XPOUT,1)
 20   CONTINUE
      RETURN
C *** Last line of DBDTRP ***
      END
