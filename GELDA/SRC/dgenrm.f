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
C====>    time          = "17:18:20 MESZ",
C====>    filename      = "dbdnrm.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "26359 119 534 4079",
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
      DOUBLE PRECISION FUNCTION DGENRM (NEQ, V, RWT)
C
C     PURPOSE
C
C     Compute the weighted root-mean-square vector norm for DGELDA [2].
C
C     ARGUMENT LIST
C
C       ARGUMENT IN
C
C       NEQ - INTEGER.
C             The actual dimension of the vectors V and RWT.
C
C         V - DOUBLE PRECISION array of DIMENSION (*).
C             This array must contain the vector V.
C
C       RWT - DOUBLE PRECISION array of DIMENSION (*).
C             This array must contain the reciprocal of weights for the
C             norm.
C
C       ARGUMENTS OUT
C
C    DGENRM - DOUBLE PRECISION array of DIMENSION (*).
C             This array contains the weighted root-mean-square norm
C             of the vector V with reciprocal weights RWT.
C
C     METHOD
C
C     This function routine computes the weighted root-mean-square
C     norm of the vector of length NEQ contained in the array V, with
C     1/weights contained in the array RWT of length NEQ.
C        DGENRM=SQRT((1/NEQ)*SUM(V(I)*RWT(I))**2)
C
C     DGENRM is a modiefied version of the DBDNRM function of DASSL [3].
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
C     1995, July 10 [Version 1.0]
C       First release. (W. Rath)
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  V(*), RWT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM, VMAX
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      DGENRM = 0.0D0
      VMAX = 0.0D0
      DO 10 I = 1,NEQ
        IF(ABS(V(I)*RWT(I)) .GT. VMAX) VMAX = ABS(V(I)*RWT(I))
 10   CONTINUE
      IF(VMAX .LE. 0.0D0) GO TO 30
      SUM = 0.0D0
      DO 20 I = 1,NEQ
         SUM = SUM + ((V(I)*RWT(I))/VMAX)**2
 20   CONTINUE
      DGENRM = VMAX*SQRT(SUM/NEQ)
 30   CONTINUE
      RETURN
C *** Last line of DGENRM ***
      END
