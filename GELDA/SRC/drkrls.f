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
C====>    time          = "13:27:45 MEZ",
C====>    filename      = "drkrls.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "36293 105 412 3647",
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
      SUBROUTINE DRKRLS(N, N3, IU, FAC1, UROUND,
     $                  AH, LDAH, W, LDW, E0, LDE0, A0, LDA0,
     $                  B, YY, IWORK, WORK, IDID)
C
C     PURPOSE
C
C     Decompose and solve the (3n*3n) real Runge-Kutta system [1]              C
C
C     REFERENCES
C
C     [1] P. Kunkel, V. Mehrmann, W. Rath and J. Weickert.
C         GELDA: A software package for the solution of general linear
C         differential algebraic equations.
C         Preprint SPC 95_8, TU Chemnitz-Zwickau, February 1995.
C
      IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      INTEGER          OIPVT, ONFA
      PARAMETER (ONFA   = 13)
      PARAMETER (OIPVT  = 21)
C     .. Scalar Arguments ..
      DOUBLE PRECISION FAC1, UROUND
      INTEGER          IDID, IU, LDAH, LDA0, LDE0, LDW, N, N3
C     .. Array Arguments ..
      DOUBLE PRECISION AH(LDAH,*), A0(LDA0,*), E0(LDE0,*), W(LDW,*),
     $                 B(*), WORK(*), YY(*)
      INTEGER          IWORK(*)
C     .. Local Scalars ..
      INTEGER          I, IER, IRANK, J
C     .. External Subroutines ..
      EXTERNAL         DCOPY, DGELSX, DGETRF, DGETRS
C     .. Executable Statements ..
C
C --- Decomposition and solution of the Runge-Kutta system
      IF (IU.EQ.0) THEN
         CALL DGETRF (N3,N3,W,LDW,IWORK(OIPVT),IER)
         IF (IER.NE.0) GOTO 98
         IWORK(ONFA)=IWORK(ONFA)+1
         CALL DGETRS ('N',N3,1,W,LDW,IWORK(OIPVT),B,N3,IER)
         IF (IER .NE. 0) GOTO 99
      ELSE
         DO 11 I=1,N3
            IWORK(OIPVT+I-1)=0
  11     CONTINUE
         CALL DGELSX(N3,N3,1,W,LDW,B,N3,IWORK(OIPVT),UROUND,IRANK,
     $               WORK,IER)
         IF (IER.NE.0) GOTO 99
         IWORK(ONFA)=IWORK(ONFA)+1
      END IF
      CALL DCOPY(N3,B,1,YY,1)
C
C --- Supply matrix of error estimation in full rank case
      IF (IU.EQ.0) THEN
         DO 13 I=1,N
            DO 12 J=1,N
               AH(I,J) = FAC1*E0(I,J) - A0(I,J)
  12        CONTINUE
  13     CONTINUE
         CALL DGETRF (N,N,AH,LDAH,IWORK(OIPVT+N3),IER)
         IF (IER.NE.0) GOTO 98
      END IF
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
C --- Last line of subroutine DRKRLS
      END
