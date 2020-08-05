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
C====>    version       = "1.0",
C====>    date          = "12 October 1995",
C====>    time          = "13:35:17 MEZ",
C====>    filename      = "dbinom.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "19104 124 447 4000",
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
      FUNCTION DBINOM(I,J)
C
C     PURPOSE
C
C     DBINOM computes
C
C             ( I )         I!
C            (     ) = ------------
C             ( J )     J! (I-J)!   .
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION DBINOM
      INTEGER          I, J
C     .. Local Arrays ..
      DOUBLE PRECISION B(10,10)
C     .. External Functions ..
      DOUBLE PRECISION DFACLN
      EXTERNAL         DFACLN
C     .. Intrinsic ..
      INTRINSIC        ANINT, EXP
C     .. Save Statements ..
      SAVE             B
C     .. Data Statements ..
      DATA             B/100*-1.0D0/
C     .. Executable Statements ..
      IF (J.LT.0 .OR. J.GT.I) THEN
         DBINOM=0.D0
      ELSE
         IF (I.LE.9 .AND. J.LE.9) THEN
            IF (B(I+1,J+1).LT.0.0D0) THEN
               B(I+1,J+1) = ANINT(EXP(DFACLN(I)-DFACLN(J)-DFACLN(I-J)))
            ENDIF
            DBINOM = B(I+1,J+1)
         ELSE
            DBINOM=ANINT(EXP(DFACLN(I)-DFACLN(J)-DFACLN(I-J)))
         ENDIF
      ENDIF
      RETURN
C *** Last line ob DBINOM
      END
      FUNCTION DFACLN(N)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      DOUBLE PRECISION DFACLN
      INTEGER          N
C     .. Local Arrays ..
      DOUBLE PRECISION A(10)
C     .. External Functions ..
      DOUBLE PRECISION DGAMLN
      EXTERNAL         DGAMLN
C     .. Save Statements ..
      SAVE             A
C     .. Data Statements ..
      DATA             A/10*-1.0D0/
C     .. Executable Statements ..
      IF (N.LT.0) PAUSE 'ERROR IN DFAC'
      IF (N.LE.9) THEN
         IF (A(N+1).LT.0.0D0) A(N+1) = DGAMLN(N+1.0D0)
         DFACLN = A(N+1)
      ELSE
         DFACLN = DGAMLN(N+1.0D0)
      ENDIF
      RETURN
C *** Last line of DFACLN ***
      END
      FUNCTION DGAMLN(XX)
      IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION STP
      PARAMETER        (STP=2.50662827465D0)
      DOUBLE PRECISION HALF, ONE, FPF
      PARAMETER        (HALF=0.5D0, ONE =1.0D0, FPF=5.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION DGAMLN, XX
C     .. Local Scalars ..
      DOUBLE PRECISION SER, TMP, X
      INTEGER          J
C     .. Local Arrays ..
      DOUBLE PRECISION COF(6)
C     .. Intrinsic Statements ..
      INTRINSIC        LOG
C     .. Data Statements
      DATA             COF / 76.18009173D0,   -86.50532033D0,
     $                       24.01409822D0,    -1.231739516D0,
     $                        0.120858003D-2,  -0.536382D-5/
C     .. Executable Statements ..
      X = XX - ONE
      TMP = X +FPF
      TMP = (X + HALF)*LOG(TMP)- TMP
      SER = ONE
      DO 10 J = 1,6
         X = X + ONE
         SER = SER + COF(J)/X
 10   CONTINUE
      DGAMLN = TMP + LOG(STP*SER)
      RETURN
      END
