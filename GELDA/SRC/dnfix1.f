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
C====>    time          = "17:06:29 MESZ",
C====>    filename      = "dnfix1.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "46770 413 2189 15325",
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
      SUBROUTINE DNFIX1(EDIF, ADIF, FDIF, N, T, LMAX, M, ID, IA, IU,
     $                  IREQ, MQ, IDQ, IAQ, IUQ, E, LDE, A, LDA, F,
     $                  EQ, LDEQ, AQ, LDAQ, FQ, Z1, LDZ1, Z2Q, LDZ2Q,
     $                  AH, LDAH, IPAR, RPAR, RWORK, LRWORK, MCONST,
     $                  IERR)
C
C     PURPOSE
C
C     DNFIX1 is part of the code DGELDA [2] which solves linear DAEs
C     with variable coefficients of the form
C
C                               .
C                           E(t)x(t) = A(t)x(t) + f(t)
C                             x(t_0) = x_0
C
C     DNFIX1 performs the following steps:
C
C     1. Compute the strangeness index MU and the other characteristic
C        values of the matrix pair (E(t), A(t)) and time T.
C     2. Compute the matricies EQ, AQ und the vector FQ of the extended
C        derivative array of stage MU at time T.
C     3. Check if the strangeness index and the other characteristic values
C        have changed.
C     4. Compute the coefficients of an equivalent starngeness free DAE
C        at time T. The two DAEs are equivalent in the sense that their
C        solutions are identical.
C
C     ARGUMENT LIST
C
C        USER-SUPPLIED SUBROUTINES
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
C        ARGUMENTS IN
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
C         M - INTEGER.
C             The strangeness index MU of the DAE.
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
C      IREQ - INTEGER.
C             The rank of the matrix EQ at the previous step.
C
C       LDE - INTEGER.
C             The leading dimension of array E as declared in the
C             calling program.
C             LDE .GE. N.
C
C       LDA - INTEGER.
C             The leading dimension of array A as declared in the
C             calling program.
C             LDA .GE. N.
C
C      LDEQ - INTEGER.
C             The leading dimension of array EQ as declared in the
C             calling program.
C             LDEQ .GE. N*(LMAX+1).
C
C      LDAQ - INTEGER.
C             The leading dimension of array AQ as declared in the
C             calling program.
C             LDAQ .GE. N*(LMAX+1).
C
C      LDZ1 - INTEGER.
C             The leading dimension of array Z1 as declared in the
C             calling program.
C             LDZ1 .GE. N.
C
C    LDZ2Q - INTEGER.
C             The leading dimension of array Z2Q as declared in the
C             calling program.
C             LDZ2Q .GE. N*(LMAX+1).
C
C      LDAH - INTEGER.
C             The leading dimension of array AH as declared in the
C             calling program.
C             LDAH .GE. N*(LMAX+1)
C
C     ARGUMENTS OUT
C
C      IREQ - INTEGER.
C             The rank of EQ at time T.
C
C        MQ - INTEGER.
C             The strangeness index MU of the DAE at time T.
C
C       IDQ - INTEGER.
C             The number DMU of differential components at time T.
C
C       IAQ - INTEGER.
C             The number DMU of algebraic components at time T.
C
C       IUQ - INTEGER.
C             The number DMU of undetermined components at time T.
C
C         E - DOUBLE PRECISION array of DIMENSION (LDE,*).
C             The leading N by N part of this array contains the matrix
C             E of the strangeness free DAE at time T.
C
C         A - DOUBLE PRECISION array of DIMENSION (LDA,*).
C             The leading N by N part of this array contains the matrix
C             A of the strangeness free DAE at time T.
C
C         F - DOUBLE PRECISION array of DIMENSION (*).
C             The leading N elements of this array contain the vector F of
C             the strangeness free DAE at time T.
C
C        EQ - DOUBLE PRECISION array of DIMENSION (LDEQ,*).
C             The leading N*(MD+1) by N*(MD+1) part of this array
C             contains the matrix EQ of stage MQ at time T.
C
C        AQ - DOUBLE PRECISION array of DIMENSION (LDAQ,*).
C             The leading N*(MD+1) by N*(MD+1) part of this array
C             contains the matrix AQ of stage MQ at time T.
C
C        FQ - DOUBLE PRECISION array of DIMENSION (*).
C             The leading N*(MD+1) elemnts of this array contain vector
C             contain the vector FQ of stage MQ at time T.
C
C        Z1 - DOUBLE PRECISION array of DIMENSION (LDZ1,*).
C             The leading N by N part of this array contains the
C             transformation matrix Z1 which extracts the differential
C             part.
C
C       Z2Q - DOUBLE PRECISION array of DIMENSION (LDZ2Q,*).
C             The leading N*(MD+1) by IAQ part of this array contains
C             the transformation matrix Z2Q which extracts the algebraic
C             part.
C
C        AH - DOUBLE PRECISION array of DIMENSION (LDAH,*).
C             This array is used as workarray.
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
C     MODE PARAMETERS
C
C     MCONST - LOGICAL.
C              Indicates whether the matricies E(t) and A(t) or constant
C              as follows:
C
C              MCONST =  .TRUE., (E(t) and A(t) are both constant
C                                 matricies);
C              MCONST = .FLASE.  (E(t) or A(t) is a time dependent
C                                 matrix).
C
C     ERROR INDICATOR
C
C       IERR - INTEGER.
C              Unless the routine detects an error (see next section),
C              IERR contains 0 on exit.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C     IERR = -1  : An error occured in EDIF, ADIF or FDIF.
C     IERR = -2  : The input parameter IDIF of EDIF, ADIF or FDIF was
C                  larger then the highest derivative the routine
C                  provides.
C     IERR = -22 : Failed to determine strangeness index.
C     IERR = -23 : Failed to compute an equivalent strangeness index 0
C                  system.
C     IERR = -26 : Unable to continue due to a change in characteristic
C                  values.
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
C     All information we need to compute the strangeness free DAE is
C     hidden in the extended derivative array EQ, AQ anf FQ.
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
C     1995, July 17 [Version 1.1]
C       Changed documentation to meet SLICOT standard.
C
C     1995, July 10 [Version 1.0]
C       First release. (W. Rath, J. Weickert)
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER        (ZERO = 0.0D0, ONE = 1.0D0)
C     .. Subroutine Arguments ..
      EXTERNAL         EDIF, ADIF, FDIF
C     .. Scalar Arguments ..
      DOUBLE PRECISION T
      INTEGER          IERR, ID, IA, IU, IDQ, IAQ, IUQ, IREQ,
     $                 LDA, LDAH, LDAQ, LDE, LDEQ, LDZ1, LDZ2Q,
     $                 LMAX, LRWORK, M, MQ, N
      LOGICAL          MCONST
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), AH(LDAH,*),
     $                 AQ(LDAQ,*), E(LDE,*), EQ(LDEQ,*),
     $                 F(*), FQ(*), RPAR(*), RWORK(*),
     $                 Z1(LDZ1,*), Z2Q(LDZ2Q,*)
      INTEGER          IPAR(*)
C     .. Local Scalars ..
      INTEGER          I, NQ
C     .. External Subroutines ..
      EXTERNAL         DCOPY, DGEMV, DLACPY, DLASET
C     .. Executable Statements ..
      NQ = N*(M+1)
C
C     Compute relevant matrices.
      IERR = 0
      IF (MCONST) THEN
         IF (M.EQ.0) THEN
            CALL FDIF(N,T,0,F,IPAR,RPAR,IERR)
            IF (IERR.NE.0) RETURN
         ELSE
            DO 10 I=0,M
               CALL FDIF(N,T,I,FQ(I*N+1),IPAR,RPAR,IERR)
               IF (IERR.NE.0) RETURN
 10         CONTINUE
            CALL DGEMV('T',N,ID,ONE,Z1,LDZ1,FQ,1,ZERO,F,1)
            CALL DGEMV('T',NQ,IA+IU,ONE,Z2Q,LDZ2Q,FQ,1,ZERO,F(ID+1),1)
         ENDIF
      ELSE
         CALL DNFIND(EDIF, ADIF, FDIF, N, T, LMAX, MQ, IDQ, IAQ, IUQ,
     $               IREQ, AH, LDAH, EQ, LDEQ, AQ, LDAQ, FQ, Z2Q, LDZ2Q,
     $               IPAR, RPAR, RWORK, LRWORK, IERR)
         IF (IERR.NE.0) THEN
            IERR = -22
            RETURN
         ENDIF
         IF (MQ.NE.M .OR. IDQ.NE.ID .OR. IAQ.NE.IA .OR. IUQ.NE.IU) THEN
            IERR = -26
            RETURN
         ENDIF
         IF (M.EQ.0) THEN
            CALL DCOPY(N,FQ,1,F,1)
            CALL DLACPY('F',N,N,EQ,LDEQ,E,LDE)
            CALL DLACPY('F',N,N,AQ,LDEQ,A,LDE)
            CALL DLASET('N',NQ,NQ,ZERO,ONE,Z2Q,LDZ2Q)
         ELSE
            CALL DNFRED(N, NQ, ID, IA, IU, IREQ, EQ, LDEQ, AQ, LDAQ,
     $                  FQ,Z2Q, LDZ2Q, Z1, LDZ1, E, LDE, A, LDA, F, AH,
     $                  LDAH,RWORK, LRWORK, IERR)
            IF (IERR.NE.0) THEN
               IERR = -23
               RETURN
            ENDIF
         ENDIF
      ENDIF
      RETURN
C *** Last line of DNFIX1 ***
      END
