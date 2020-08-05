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
C====>    time          = "18:03:50 MESZ",
C====>    filename      = "dgeini.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "38493 115 520 4425",
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
      SUBROUTINE DGEINI(EDIF, ADIF, FDIF, NEQ, TN, MXINDX, X, XPRIME,
     $                  M, ID, IA, IU, IREQ, EQ, LDEQ, AQ, LDAQ, FQ,
     $                  E, LDE, A, LDA, F, Z1, LDZ1, Z2Q, LDZ2Q, D, AH,
     $                  LDAH, V, IPAR, RPAR, RWORK, LRWORK, INIVAL,
     $                  IERR)
C
C     PURPOSE
C
C     Compute strangeness index and characteristic values of a linear DAE
C     at time TN. Then get a strangeness free matrix pair at time T.
C     If INIVAL=1 use this pair to compute consistent initial values X.
C     Finnally compute initial derivative XPRIME.
C
C     DGEINI is only for internal use and is called from DGELDA.
C
C     CONTRIBUTORS
C
C     W. Rath, J. Weickert (TU Chemnitz, Germany).
C
C     REVISIONS
C
C     1995, July 15 [Version 1.1]
C       Changed interface to meet SLICOT standard. (W. Rath)
C
C     1995, July 10 [Version 1]
C       First release. (W. Rath)
C
C     ******************************************************************
C
C     .. Subroutine Arguments ..
      EXTERNAL         ADIF, EDIF, FDIF
C     .. Scalar Arguments ..
      DOUBLE PRECISION TN
      INTEGER          IA, ID, IERR, INIVAL, IREQ, IU, LDA, LDAQ, LDEQ,
     $                 LDE, LDZ1, LDZ2Q, LDAH, LRWORK, MXINDX, M, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), AH(LDAH,*), AQ(LDAQ,*), D(*),
     $                 E(LDE,*), EQ(LDEQ,*), F(*), FQ(*),
     $                 RPAR(*), RWORK(*), V(*),  X(*), XPRIME(*),
     $                 Z1(LDZ1,*), Z2Q(LDZ2Q,*)
      INTEGER          IPAR(*)
C     .. Local Scalars
      INTEGER          IDID, NQ
C     .. External Subroutines ..
      EXTERNAL         DCOPY, DGEMV, DGELS, DLACPY, DNFIND, DNFINI,
     $                 DNFRED
C     .. Executable statements ..
C
      IDID = 0
C
C     Compute index and get reduced form at starting point.
      CALL DNFIND(EDIF, ADIF, FDIF, NEQ, TN, MXINDX, M, ID, IA, IU,
     $            IREQ, AH, LDAH, EQ, LDEQ, AQ, LDAQ, FQ, Z2Q, LDZ2Q,
     $            IPAR, RPAR, RWORK, LRWORK, IDID)
      IF (IDID.NE.0 .OR. M.LT.0) THEN
         IERR = -22
         GOTO 320
      ENDIF
      NQ = (M+1)*NEQ
      CALL DNFRED(NEQ, NQ, ID, IA, IU, IREQ, EQ, LDEQ, AQ, LDAQ, FQ,
     $            Z2Q, LDZ2Q, Z1, LDZ1, E, LDE, A, LDA, F, AH, LDAH,
     $            RWORK, LRWORK, IDID)
      IF (IDID.NE.0) THEN
         IERR = -23
         GOTO 320
      ENDIF
C
C     Compute initail value if necessary and compute initial derivative.
      IF (INIVAL .EQ. 0) THEN
         CALL DLACPY('N',NEQ,NEQ,E,LDE,EQ,LDEQ)
         CALL DLACPY('N',NEQ,NEQ,A,LDA,AQ,LDAQ)
         CALL DNFINI(NEQ,ID,IA,IU,X,EQ,LDEQ,AQ,LDAQ,F,V,D,RWORK,LRWORK,
     $               0,IDID)
         IF (IDID.NE.0) THEN
            IERR = -24
            GOTO 320
         ENDIF
      ENDIF
      CALL DCOPY(NEQ,F,1,XPRIME,1)
      CALL DGEMV('N',ID,NEQ,1.0D0,A,LDA,X,1,1.0D0,XPRIME,1)
      CALL DLACPY('N',NEQ,NEQ,E,LDE,AH,LDAH)
      CALL DGELS('N',ID,NEQ,1,AH,LDAH,XPRIME,NEQ,RWORK,LRWORK,
     $     IDID)
      IF (IDID.NE.0) THEN
         IERR = -25
      END IF
C *** Last line of DGRINI ***
 320  END
