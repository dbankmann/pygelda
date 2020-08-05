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
C====>    date          = "19 July 1995",
C====>    time          = "18:36:36 MESZ",
C====>    filename      = "dgeerm.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "39122 193 1130 8186",
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
      SUBROUTINE DGEERM (EDIF, ADIF, FDIF, NEQ, T, TOUT, X, XPRIME,
     $                   CVAL, IPAR, RPAR, IWORK, LIW, RWORK, LRW,
     $                   RTOL, ATOL, METHOD, INFO, IWARN, IERR)
C
C     PURPOSE
C
C     Prints information about an error which DGELDA [1] detected.
C
C     If IERR=-998 or IERR=-999 the run is terminated.
C
C     ARGUMENT LIST
C
C     The argument list is identicall to the argument list of DGELDA.
C
C     REFERENCES
C
C     [1] P. Kunkel, V. Mehrmann, W. Rath and J. Weickert.
C         GELDA: A software package for the solution of general linear
C         differential algebraic equations.
C         Preprint SPC 95_8, TU Chemnitz-Zwickau, February 1995.
C
C     CONTRIBUTORS
C
C     W. Rath (TU Chemnitz, Germany).
C
C     REVISIONS
C
C     1995, July 19 [Version 1.1]
C       Changed order of IPAR, RPAR, and IWORK, RWORK to meet
C       SLICOT interface standard. (W. Rath)
C
C     1995, July 10 [Version 1]
C       First release. (W. Rath)
C
C     ******************************************************************
C
C     .. Subroutine Arguments ..
      EXTERNAL         EDIF, ADIF, FDIF
C     .. Scalar Arguments ...
      DOUBLE PRECISION T, TOUT
      INTEGER          IERR, IWARN, LIW, LRW, NEQ, METHOD
C     .. Array Arguments ...
      DOUBLE PRECISION ATOL(*), RTOL(*), RPAR(*), RWORK(*), X(*),
     $                 XPRIME(*)
      INTEGER          CVAL(4), INFO(20), IPAR(*), IWORK(*)
C     .. Executable Statements ..
      IF (IERR.GE.0) RETURN
      PRINT *
      PRINT *,' *********************** ERROR IN DGELDA',
     $     ' ***********************'
      PRINT *,' ERROR NUMBER = ', IERR
      IF (IERR .EQ. -1) THEN
         PRINT *,' AT CURRENT T = ',RWORK(4)
         PRINT *,' 500 STEPS TAKEN ON THIS BEFORE REACHING TOUT'
      ELSEIF (IERR .EQ. -2) THEN
         PRINT *,' AT T = ',RWORK(4)
         PRINT *,' TOO MUCH ACCURACY REQUESTED FOR PRECISION OF',
     $        ' MACHINE.'
         PRINT *,' RTOL AND ATOL WERE INCREASED TO APPROPRIATE',
     $        ' VALUES'
      ELSEIF (IERR .EQ. -3) THEN
         PRINT *,' AT T = ',RWORK(4),' SOME ELEMENT OF WT HAS BECOME',
     $        ' .LE. 0.0'
      ELSEIF (IERR .EQ. -6) THEN
         PRINT *,' AT T = ',RWORK(4),' AND STEPSIZE H = ',RWORK(3)
         PRINT *,' THE ERROR TEST FAILED REPEATEDLY OR',
     $        ' WITH ABS(H)=HMIN'
      ELSEIF (IERR .EQ. -8) THEN
         PRINT *,' AT T = ',RWORK(4),' AND STEPSIZE H = ',RWORK(3)
         PRINT *,' THE LINEAR SYSTEM SOLVER FAILED SEVERAL TIMES'
      ELSEIF (IERR .EQ. -10) THEN
         PRINT *,' AT T = ',RWORK(4),' AND STEPSIZE H = ',RWORK(3)
         PRINT *,' THE INTEGRATION COULD NOT CONTINUE BECAUSE IERR IN'
         PRINT *,' EDIF, ADIF or FDIF WAS EQUAL TO MINUS ONE'
      ELSEIF (IERR .EQ. -21) THEN
         PRINT *,' AT T = ',RWORK(4),' AND STEPSIZE H = ',RWORK(3)
         PRINT *,' PARAMETER IDIF OF EDIF, ADIF OR FDIF  WAS EQUAL',
     $        ' TO MINUS TWO'
      ELSEIF (IERR .EQ. -22) THEN
         PRINT *,' AT T = ',RWORK(4),' AND STEPSIZE H = ',RWORK(3)
         PRINT *,' FAILED TO DETERMINE STRANGENESS INDEX'
      ELSEIF (IERR .EQ. -23) THEN
         PRINT *,' AT T = ',RWORK(4),' AND STEPSIZE H = ',RWORK(3)
         PRINT *,' FAILED TO COMPUTE EQUIVALENT STRANGENESS INDEX',
     $        ' 0 SYSTEM'
      ELSEIF (IERR .EQ. -24) THEN
         PRINT *,' AT T = ',RWORK(4),' AND STEPSIZE H = ',RWORK(3)
         PRINT *,' FAILED TO COMPUTE INITIAL X'
      ELSEIF (IERR .EQ. -25) THEN
         PRINT *,' AT T = ',RWORK(4),' AND STEPSIZE H = ',RWORK(3)
         PRINT *,' FAILED TO COMPUTE INITIAL DERIVATIVE'
      ELSEIF (IERR .EQ. -26) THEN
         PRINT *,' AT T =', RWORK(4),' AND STEPSIZE H = ',RWORK(3)
         PRINT *,
     $    ' UNABLE TO CONTINUE DUE TO CHANGE IN CHARACTERISTIC VALUES'
         PRINT *,' STRANGENESS INDEX        ',IWORK(LIW)
         PRINT *,' DIFFERENTIAL COMPONENTS  ',IWORK(LIW-1)
         PRINT *,' ALGEBRAIC    COMPONENTS  ',IWORK(LIW-2)
         PRINT *,' UNDETERMINED COMPONENTS  ',IWORK(LIW-3)
      ELSEIF (IERR .EQ. -101) THEN
         PRINT *, ' SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE'
      ELSEIF (IERR .EQ. -102) THEN
         PRINT *, ' NEQ = ',NEQ,'  .LE. 0'
      ELSEIF (IERR .EQ. -103) THEN
         PRINT *, ' MAXORD = ',IWORK(3),' NOT IN RANGE'
      ELSEIF (IERR .EQ. -104) THEN
         PRINT *, ' MXINDX = ',IWORK(4),' IS .LT. 0'
      ELSEIF (IERR .EQ. -105) THEN
         WRITE (*,*) ' RWORK LENGTH NEEDED, LENRW = ', IWORK(LIW),
     $        ' EXCEEDS LRW = ',LRW
      ELSEIF (IERR .EQ. -106) THEN
         PRINT *, ' IWORK LENGTH NEEDED, LENIW = ', IWORK(LIW),
     $        ' EXCEEDS LIW = ',LIW
      ELSEIF (IERR .EQ. -107) THEN
         PRINT *, ' SOME ELEMENT OF RTOL IS .LT. 0'
      ELSEIF (IERR .EQ. -108) THEN
         PRINT *, ' SOME ELEMENT OF ATOL IS .LT. 0'
      ELSEIF (IERR .EQ. -109) THEN
         PRINT *, ' ALL ELEMENTS OF RTOL AND ATOL ARE ZERO'
      ELSEIF (IERR .EQ. -110) THEN
         PRINT *, ' INFO(4) = 1 AND TSTOP = ',RWORK(1),' BEHIND',
     $        ' TOUT = ',TOUT
      ELSEIF (IERR .EQ. -111) THEN
         PRINT *, ' HMAX = ',RWORK(2),' .LE. 0.0'
      ELSEIF (IERR .EQ. -112) THEN
         PRINT *, ' TOUT = ',TOUT,' BEHIND T = ',T
      ELSEIF (IERR .EQ. -113) THEN
         PRINT *, ' INFO(8)=1 AND H0=0.0'
      ELSEIF (IERR .EQ. -114) THEN
         PRINT *, ' SOME ELEMENT OF WT IS .LE. 0.0'
      ELSEIF (IERR .EQ. -115) THEN
         PRINT *, ' TOUT = ',TOUT,' TOO CLOSE TO T = ',T,
     $        ' TO START INTEGRATION'
      ELSEIF (IERR .EQ. -116) THEN
         PRINT *, ' INFO(4)=1 AND TSTOP = ',RWORK(1),' BEHIND T = ',T
      ELSEIF (IERR .EQ. -119) THEN
         PRINT *, ' TOUT = T = ', T
      ELSEIF (IERR .EQ. -120) THEN
         PRINT *, ' WRONG INPUT FOR NMAX, IWORK(20) = ',IWORK(20)
      ELSEIF (IERR .EQ. -121) THEN
         PRINT *, ' CURIOUS INPUT FOR SAFE, RWORK(11) = ',RWORK(11)
      ELSEIF (IERR .EQ. -122) THEN
         PRINT *, ' CURIOUS INPUT FOR FACL OR FACR, RWORK(12) = ',
     $         RWORK(12),' RWORK(13) = ',RWORK(13)
      ELSEIF (IERR .EQ. -123) THEN
         PRINT *, ' CURIOUS INPUT FOR QUOT1 OR QUOT2, RWORK(14) = ',
     $         RWORK(14),' RWORK(15) = ',RWORK(15)
      ELSEIF (IERR .EQ. -998) THEN
         PRINT *, ' THE LAST STEP TERMINATED WITH A NEGATIVE VALUE OF',
     $        ' IERR'
         PRINT *, ' AND NO APPROPRIATE ACTION WAS TAKEN.'
         PRINT *, ' RUN TERMINATED. '
         STOP
      ELSEIF (IERR .EQ. -999) THEN
         PRINT *, ' REPEATED OCCURRENCES OF ILLEGAL INPUT$$'
         PRINT *, ' RUN TERMINATED. APPARENT INFINITE LOOP'
         STOP
      ELSE
         PRINT *,' UNKNOWN ERROR !!'
         PRINT *,' PLEASE CONTACT THE AUTHOR !!'
         STOP
      ENDIF
      PRINT *,' ***************************************',
     $     '************************'
      PRINT *
      END
