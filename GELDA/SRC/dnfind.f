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
C====>    version       = "1.1.1",
C====>    date          = "12 October 1995",
C====>    time          = "11:16:37 MEZ",
C====>    filename      = "dnfind.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "39237 348 1768 12421",
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
      SUBROUTINE DNFIND(EDIF, ADIF, FDIF, N, T, LMAX, M, IDF, IAF, IUF,
     $                  IREQ, AH, LDAH, EQ, LDEQ, AQ, LDAQ, FQ, ZQ,
     $                  LDZQ, IPAR, RPAR, RWORK, LRWORK, IERR)
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
C     DNFIND computes the strangeness index MU and the characteristic
C     values of the matrix pair (E(t), A(t)) at time T and the right
C     orthogonal factor ZQ of the SVD EQ = ZQ*S*VT. EQ is the left
C     matrix of the extended derivative array of stage MU.
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
C      LMAX - INTEGER.
C             The maximal strangeness index.
C             LMAX .GE. 1.
C
C      LDAH - INTEGER.
C             The leading dimension of array AH as declared in the
C             calling program.
C             LDAH .GE. (LMAX+1)*N.
C
C      LDEQ - INTEGER.
C             The leading dimension of array EQ as declared in the
C             calling program.
C             LDEQ .GE. (LMAX+1)*N.
C
C      LDAQ - INTEGER.
C             The leading dimension of array AQ as declared in the
C             calling program.
C             LDAQ .GE. (LMAX+1)*N.
C
C      LDZQ - INTEGER.
C             The leading dimension of array ZQ as declared in the
C             calling program.
C             LDZQ .GE. (LMAX+1)*N.
C
C       ARGUMENTS OUT
C
C         M - INTEGER.
C             The strangeness index MU of the DAE.
C
C       IDF - INTEGER.
C             The number DMU of differential components.
C
C       IAF - INTEGER.
C             The number DMU of differential components.
C
C       IUF - INTEGER.
C             The number UMU of undetermined components.
C
C      IREQ - INTEGER.
C             The rank of EQ.
C
C        AH - DOUBLE PRECISION array of DIMENSION (LDAH,*).
C             The leading N*(M+1) by N part of this array contains
C              (j)
C             E  (T), j=0,...,M.
C
C        EQ - DOUBLE PRECISION array of DIMENSION (LDEQ,*).
C             The leading N*(M+1) by N*(M+1) part of this array contains
C             the matrix EQ of stage M.
C
C        AQ - DOUBLE PRECISION array of DIMENSION (LDAQ,*).
C             The leading N*(M+1) by N*(M+1) part of this array contains
C             the matrix AQ of stage M.
C
C        FQ - DOUBLE PRECISION array of DIMENSION (*).
C             The first N*(M+1) components of this array contain the
C             vector FQ of stage L.
C
C        ZQ - DOUBLE PRECISION array of DIMENSION (LDZQ,*).
C             The leading N*(M+1) by N*(M+1) part of this array
C             contains the right orthogonal factor of the SVD
C             EQ=ZQ*S*VT.
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
C     WORK SPACE
C
C      RWORK - DOUBLE PRECISION array of DIMENSION at least (LRWORK).
C
C     LRWORK - NTEGER.
C              The length of RWORK.
C              LRWORK .GE. 3*((LMAX+1)*N)^2 + 6*(LMAX+1)*N.
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
C     IERR = -1 : An error occured in EDIF, ADIF or FDIF.
C     IERR = -2 : The input parameter IDIF of EDIF, ADIF or FDIF was
C                 larger then the highest derivative the routine
C                 provides.
C     IERR = -3 : An argument of DGESVD had an illegal value.
C     IERR = -4 : DGESVD failed to converge.
C     IERR = -5 : LRWORK < 3*((LMAX+1)*N)^2 + 6*(LMAX+1)*N.
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
C     The characteristic values of the matrix pair (EQ,AQ) of stage I are
C     defined by
C
C                RI = rank EQ
C                AI = rank (Z^* AQ T)
C                SI = rank(V^* Z^* AQ T')
C                DI = RI - SI
C                UI = N*(I+1) - RI - AI - SI
C
C     where
C
C                T basis of kernel EQ
C                Z basis of corange EQ
C                T' basis of cokernel EQ
C                V basis of corange(Z^* AQ T).
C
C     Recursive formulas are used to compute the strangeness index MU
C     and the characteristics of the original pair (E(t), A(t)) at time
C     T in terms of the above values.
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
C     1995, October 12 [Version 1.1.1]
C       Now we check if the characteristic values are greater then zero.
C       (W. Rath, J. Weickert)
C
C     1995, July 17 [Version 1.1]
C       Changed documentation to meet SLICOT standard. (W. Rath)
C
C     1995, July 10 [Version 1.0]
C       First release. (W. Rath, J. Weickert)
C
C     ******************************************************************
C
C     .. Subroutine Arguments ..
      EXTERNAL         EDIF, ADIF, FDIF
C     .. Scalar Arguments ..
      DOUBLE PRECISION T
      INTEGER          IERR, IDF, IAF, IREQ, IUF,
     $                 LDAQ, LDEQ, LDAH, LDZQ, LMAX, LRWORK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION AQ(LDAQ,*), EQ(LDEQ,*), FQ(*), RPAR(*), RWORK(*),
     $                 AH(LDAH,*), ZQ(LDZQ,*)
      INTEGER          IPAR(*)
C     .. Local Scalars
      INTEGER          L, NQ, LFIN, IRQ, IAQ, ISQ, IDQ, IUQ,
     $                 IR, IA, IS, ISM1, ID, IU, IW, IC, IV
C     .. External Subroutines ..
      EXTERNAL         DNFPEN, DNFVAL
C     .. Executable Statements ..
      L=0
      NQ=N
      CALL DNFPEN(EDIF,ADIF,FDIF,N,T,L,AH,LDAH,
     $            EQ,LDEQ,AQ,LDAQ,FQ,IPAR,RPAR,IERR)
      IF (IERR.NE.0) RETURN
      CALL DNFVAL(N,NQ,EQ,LDEQ,AQ,LDAQ,IRQ,IAQ,ISQ,IDQ,IUQ,
     $            ZQ,LDZQ,RWORK,LRWORK,IERR)
      IF (IERR.NE.0) RETURN
      IR=IRQ
      IA=IAQ
      IS=ISQ
      ID=IDQ
      IU=IUQ
      IW=IUQ
      IC=IAQ+ISQ
      IV=IAQ+ISQ
      LFIN=0
      IF (IS.EQ.0) GOTO 2000
      IF (LMAX.GT.0) THEN
        DO 1000 L=1,LMAX
          NQ=(L+1)*N
          CALL DNFPEN(EDIF,ADIF,FDIF,N,T,L,AH,LDAH,
     $                EQ,LDEQ,AQ,LDAQ,FQ,IPAR,RPAR,IERR)
          IF (IERR.NE.0) RETURN
          CALL DNFVAL(N,NQ,EQ,LDEQ,AQ,LDAQ,IRQ,IAQ,ISQ,IDQ,
     $	              IUQ,ZQ,LDZQ,RWORK,LRWORK,IERR)
          IF (IERR.NE.0) RETURN
	  ISM1 = IS
          IR   = IR - ISM1
          IS   = ISQ - IV
          IW   = ISM1 - IS - IAQ
          IU   = IW + IU
          IA   = N - IR - IS - IU
          ID   = IR - IS
          IC   = ISM1 - IW
          IV   = IV + IC
          LFIN=L
          IF (IS .LT. 0 .OR. IW .LT. 0 .OR. IA .LT. 0 .OR.
     $        ID .LT. 0 .OR. IC .LT. 0) THEN
             IERR = -100
             RETURN
          ENDIF
          IF (IS.EQ.0) GOTO 2000
 1000   CONTINUE
      END IF
      M=-1
      RETURN
 2000 CONTINUE
      M=LFIN
      IDF=ID
      IAF=IA
      IUF=IU
      IREQ = IRQ
      RETURN
C *** Last line of DNFIND ***
      END
