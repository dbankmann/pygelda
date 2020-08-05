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
C====>    date          = "14 July 1995",
C====>    time          = "17:07:03 MESZ",
C====>    filename      = "dgelda.f",
C====>    address       = "Fakultaet fuer Mathematik
C====>                     TU Chemnitz-Zwickau
C====>                     D-09107 Chemnitz
C====>                     FRG",
C====>    telephone     = "(049) (0)371-531-3953",
C====>    FAX           = "(049) (0)371-531-2657",
C====>    checksum      = "41503 1492 7595 57676",
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
      SUBROUTINE DGELDA (EDIF, ADIF, FDIF, NEQ, T, TOUT, X, XPRIME,
     $                   CVAL, IPAR, RPAR, IWORK, LIW, RWORK, LRW,
     $                   RTOL, ATOL, METHOD, INFO, IWARN, IERR)
C
C     PURPOSE
C
C     DGELDA solves linear differentiell algebraic equations (DAEs) with
C     variable coefficients of the form
C
C                           .
C                       E(t)x(t) = A(t)x(t) +f(t)
C                          x(t0) = x0.
C
C     for x in a specified range of the independent variable t.
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
C             The number of equations.
C             N .GE. 1.
C
C         T - DOUBLE PRECISION.
C             The initial point of the integration.
C             NOTE that this scalar is overwritten.
C
C      TOUT - DOUBLE PRECISION.
C             The point at which a solution is desired. Integration
C             either forward in T (TOUT > T) or backward in T (TOUT < T)
C             is allowed.
C
C             At the beginning of the integration (INFO(1)=0, see below)
C             the user can set T = TOUT. The code will then compute the
C             strangeness index and the number of differential,
C             algebraic and undetermined equations and return this
C             values in CVAL. Furthermore, if INFO(11)=0 consistent
C             initial values will be computed and stored in X.
C
C         X - DOUBLE PRECISION array of DIMENSION (N).
C             If INFO(11)=0, this array can contain a guess for the
C             initial value. A consistent initial value close (in the
C             least square sense) to this guess is then computed. If no
C             guess is available, set all elements of X to zero.
C
C             If INFO(11)=1, this array must contain consistent initial
C             values of the N solution components at the initial point
C             T.
C             NOTE that this array is overwritten.
C
C       ARGUMENTS OUT
C
C         T - DOUBLE PRECISION.
C             The solution was successfully advanced to the output value
C             of T.
C
C         X - DOUBLE PRECISION array of DIMENSION (N).
C             Contains the computed solution approximation at T.
C
C    XPRIME - DOUBLE PRECISION array of DIMENSION (N).
C             Contains the computed first derivative of the solution
C             approximation at T.
C      CVAL - INTEGER array of DIMENSION (4).
C             Contains the characteristic values of the DAE:
C             CVAL(1) contains the strangeness index MU.
C             CVAL(2) contains the number DMU of differential
C                     components.
C             CVAL(3) contains the number AMU of algebraic components.
C             CVAL(4) contains the number UMU of undetermined
C                     components.
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
C             IPAR and RPAR are not altered by DGELDA or its
C             subprograms. If these arrays are used, they must be
C             dimensioned in the calling program and in EDIF, ADIF and
C             FDIF as arrays of appropriate length. Otherwise, ignore
C             these arrays by treating them as dummy arrays of length
C             one.
C
C     RWORK, IWORK - See WORKSPACE
C             These real and integer arrays provide the workspace (see
C             below).
C             Usually the information they contain are of no interest,
C             but sometimes the following may be useful:
C             IWORK(7)  contains the order of the method to be attempted
C                       on the next step.
C             IWORK(8)  contains the order of the method used on the last
C                       step.
C             IWORK(11) contains the number of steps taken so far.
C             IWORK(12) contains the number of calls to EDIF, ADIF and
C                       FDIF so far.
C             IWORK(13) contains the number of factorizations of the
C                       system matrix so far.
C             IWORK(14) contains the total number of error test failures
C                       so far.
C             RWORK(3)  contains the stepsize H to be attempted on the
C                       next step.
C             RWORK(4)  contains the current value of the independent
C                       variable, i.e., the farthest point integration
C                       has reached. This will be different from T only
C                       when interpolation has been performed (IERR=3).
C             RWORK(7)  contains the stepsize used on the last successful
C                       step.
C
C     WORK SPACE
C
C     IWORK - INTEGER array of DIMENSION at least (LIW).
C
C       LIW - INTEGER.
C             The length of IWORK.
C             If METHOD = 1 then LIW .GE. 20+2*N.
C             If METHOD = 2 then LIW .GE. 20+6*N.
C
C      RWORK - DOUBLE PRECISION array of DIMENSION at least (LRW).
C
C       LRW - NTEGER.
C             The length of RWORK.
C             If METHOD = 1 then
C               LRW .GE. 50 + (21 + MAXORD + 14*MXINDX)*N
C                    + (6 + MXINDX + 6*(MXINDX+1)^2)*N^2.
C             If METHOD = 2 then
C               LRW .GE. 20 + (29 + 14*MXINDX)*N
C                    + (15 + MXINDX + 6*(MXINDX+1)^2)*N^2.
C             NOTE that for good performance, LRW should generally be
C             larger.
C
C     TOLERANCES
C
C      RTOL - DOUBLE PRECISION array of DIMENSION (*)
C             The relative error tolerances which the user provides to
C             indicate how accurately he wishs the solution to be
C             computed. The user may choose RTOL and ATOL to be both
C             scalars or else both vectors.
C             If RTOL and ATOL are scalars (INFO(2)=0) the user has to
C             declare this array to be RTOL(1).
C
C      ATOL - DOUBLE PRECISION array of DIMENSION (*)
C             The absolute error tolerances which the user provides.
C             If ATOL and RTOL are scalars (INFO(2)=0) the user has to
C             declare  this array to be ATOL(1).
C
C             The tolerances are used by the code in a local error test
C             at each step which requires roughly that
C             |LOCAL ERROR| .LE. RTOL*|X| + ATOL
C             for each vector component.
C
C             The true (global) error is the difference between the true
C             solution of the initial value problem and the computed
C             approximation. Practically all present day codes,
C             including this one, control the local error at each step
C             and do not even attempt to control the global error directly.
C
C             Usually, but not always, the true accuracy of the computed
C             X is comparable to the error tolerances. This code will
C             usually, but not always, deliver a more accurate solution
C             if the user reduces the tolerances and  integrate again.
C             By comparing two such solutions the user can get a fairly
C             reliable idea of the true error in the solution at the
C             bigger tolerances.
C
C             Setting ATOL=0 results in a pure relative error test on
C             that component. Setting RTOL=0 results in a pure absolute
C             error test on that component. A mixed test with non-zero
C             RTOL and ATOL corresponds roughly to a relative error test
C             when the solution component is much bigger than ATOL and
C             to an absolute error test when the solution component is
C             smaller than the threshold ATOL.
C
C             The code will not attempt to compute a solution at an
C             accuracy unreasonable for the machine being used. It will
C             advise the user if he asks for too much accuracy and
C             inform the user as to the maximum accuracy it believes
C             possible.
C
C     MODE PARAMETERS
C
C    METHOD - INTEGER.
C             Indicates which integration method should be used as
C             follows:
C             METHOD=1 the code uses the BDF solver.
C             METHOD=2 the code uses the Runge-Kutta solver.
C
C      INFO - INTEGER array of DIMENSION (20)
C             The basic task of the code is to solve the system from T
C             to TOUT and return an answer at TOUT. INFO is an integer
C             array which is used to communicate exactly how the user
C             wants this task to be carried out. The simplest use of the
C             code corresponds to setting all entries of INFO to 0 (See
C             below).
C
C     WARNING INDICATOR
C
C     IWARN - INTEGER.
C             Is always zero in this version of DGELDA.
C
C     ERROR INDICATOR
C
C      IERR - INTEGER.
C             Unless the code detects an error (see next section), IERR
C             contains a positive value on exit.
C             IERR = 1: In the BDF-solver (METHOD=1) a step was
C                       successfully taken in the intermediate-output
C                       mode. The code has not yet reached TOUT.
C             IERR = 2: The integration to TOUT was successfully completed
C                       (T=TOUT) by stepping exactly to TOUT.
C             IERR = 3: In the BDF-solver the integration to TOUT was
C                       successfully completed (T=TOUT) by stepping past
C                       TOUT. X is obtained by interpolation.
C             IERR = 4: At the first step the user set T=TOUT and the code
C                       computed successfully the strangeness index and
C                       the other characteristic quantities. Furthermore,
C                       if INFO(11)=0, consistent initial values are
C                       stored in X.
C
C     WARNINGS AND ERRORS DETECTED BY THE ROUTINE
C
C                       *** Task interrupted ***
C
C     IERR =  -1 : A large amount of work has been expended (More then
C                  NMAX steps).
C     IERR =  -2 : The error tolerances are too stringent.
C     IERR =  -3 : The local error test cannot be satisfied because the
C                  user specified a zero component in ATOL and the
C                  corresponding computed solution component is
C                  zero. Thus, a pure relative error test is impossible
C                  for this component.
C     IERR =  -6 : DGELDA had repeated error test failures on the last
C                  attempted step.
C     IERR =  -8 : The linear system solver failed several times.
C     IERR = -10 : The integration could not continue because IERR in
C                  EDIF, ADIF or FDIF was repeatedly equal to minus 1.
C     IERR = -21 : The input parameter IDIF of EDIF, ADIF or FDIF was
C                  larger then the highest derivative the routine
C                  provides.
C     IERR = -22 : Failed to determine strangeness index.
C     IERR = -23 : Failed to compute an equivalent strangeness index 0
C                  system.
C     IERR = -24 : Failed to compute an initial X.
C     IERR = -25 : Failed to compute an initial derivative.
C     IERR = -26 : Unable to continue due to change in characteristic
C                  values.
C
C                       *** Task terminated ***
C
C     IERR = -101 : Some element of INFO vector is not zero or one.
C     IERR = -102 : N .LE. 0.
C     IERR = -103 : MAXORD not in range.
C     IERR = -104 : MXINDX < 0.
C     IERR = -105 : LRW is less than the required length for RWORK.
C     IERR = -106 : LIW is less than the required length for IWORK.
C     IERR = -107 : Some element of RTOL is < 0.
C     IERR = -108 : Some element of ATOL is < 0.
C     IERR = -109 : All elements of RTOL and ATOL are zero.
C     IERR = -110 : INFO(4)=1 and TOUT is behind TSTOP.
C     IERR = -111 : HMAX .LE. 0.0
C     IERR = -112 : TOUT is behind T.
C     IERR = -113 : INFO(8)=1 and H0=0.0.
C     IERR = -114 : Some element of RTOL*|X| + ATOL is .LE. 0.0.
C     IERR = -115 : TOUT is too close to T to start integration.
C     IERR = -116 : INFO(4)=1 and T is behind TSTOP.
C     IERR = -119 : INFO(1) = 1 and TOUT = T.
C     IERR = -120 : NMAX .LE. 0.
C     IERR = -121 : SAFE .LE. 0.001 or SAFE .GE. 1.0.
C     IERR = -122 : FACL < 1.0 or FACR > 1.0.
C     IERR = -123 : QUOT1 > 1.0 or QUOT2 < 1.0.
C     IERR = -998 : The last step was terminated with a negative value
C                   of IERR larger than -100, and no appropriate action
C                   was taken.
C     IERR = -999 : The previous call was terminated because of illegal
C                   input (IERR<-100) and there is illegal input in the
C                   present call, as well. (Suspect infinite loop.)
C
C     METHOD
C
C     The most important invariant in the analysis of linear DAE's is
C     the so called STRANGENESS INDEX [4], which generalizes the
C     differential index [1] for systems with undetermined components.
C
C     The implementation of DGELDA [5] is based on the construction of
C     the discretization scheme introduced in [3], which first
C     determines all the local invariants and then transforms the system
C     into a strangeness-free DAE with the same solution set.
C
C     The strangeness-free DAE is solved by either BDF methods, which
C     were adapted from DASSL of Petzold [6], or a Runge-Kutta method,
C     which was adapted from RADAU5 of Hairer/Wanner [2].
C
C     REFERENCES
C
C     [1] K. E. Brenan, S. L. Campbell and L. R. Petzold.
C         Numerical Solution of Initial-Value Problems in Differential
C         Algebraic Equations.
C         Elsevier, North Holland, New York, N. Y., 1989.
C
C     [2] E. Hairer and G. Wanner.
C         Solving Ordinary Differential Equations II.
C         Springer-Verlag, Berlin, 1991.
C
C     [3] P. Kunkel and V. Mehrmann.
C         Canonical forms for linear differential-algebraic equations
C         with variable coefficients.
C         J. Comput. Appl. Math., 56:225--259, 1994.
C
C     [4] P. Kunkel and V. Mehrmann.
C         A new class of discretization methods for the solution of
C         linear differential-algebraic equations.
C         Materialien LXII , FSP Mathematisierung, Universitaet
C         Bielefeld,  1992.
C         To appear in SIAM J. Numer. Anal.
C
C     [5] P. Kunkel, V. Mehrmann, W. Rath and J. Weickert.
C         GELDA: A software package for the solution of general linear
C         differential algebraic equations.
C         Preprint SPC 95_8, TU Chemnitz-Zwickau, February 1995.
C
C     [6] L. R. Petzold.
C         A description of DASSL: A differential/algebraic system solver.
C         In R. S. Stepleman et al., editors, IMACS Trans. Scientific
C         Computing Vol. 1, pages 65--68. North-Holland, Amsterdam, 1983.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     FURTHER COMMENTS
C
C       SETTING UP THE INFO ARRAY BEFORE THE FIRST CALL
C
C         This array is used to give the code more details about how the
C         user wants his problem to be solved. The user must respond to
C         all of the following items, which are arranged as
C         questions. The simplest use of the code corresponds to
C         answering all questions as yes, i.e. setting all entries of
C         INFO to 0.
C
C         INFO(1) This parameter enables the code to initialize
C                 itself. The user must set it to indicate the start of
C                 every new problem.
C                 Is this the first call for this problem ...
C                 Yes - Set INFO(1)=0
C                 No  - Not applicable here.
C                       See below for continuation calls.
C
C         INFO(2) The error tolerances RTOL and ATOL are used to specify
C                 how much  accuracy the user wants. The simplest use is
C                 to take them both to be scalars. To obtain more
C                 flexibility, they can both be vectors.
C                 Are both error tolerances RTOL, ATOL scalars ...
C                 Yes - Set INFO(2) = 0
C                       and input scalars for both RTOL and ATOL
C                 No  - Set INFO(2) = 1
C                       and input arrays for both RTOL and ATOL
C
C         INFO(3) The code integrates from T in the direction of TOUT by
C                 steps. If the user wishs, it will return the computed
C                 solution at the next intermediate step (the
C                 intermediate-output mode) or TOUT, whichever comes
C                 first. If the user must have solutions at many
C                 specific TOUT points, this code will compute them
C                 efficiently.
C                 Do you want the solution only at TOUT (and not at the
C                 next intermediate step) ...
C                 Yes - Set INFO(3) = 0
C                 No  - Set INFO(3) = 1.
C
C         INFO(4) If the user uses the BDF--solver (METHOD=1), the code
C                 may integrate past TOUT and interpolate to obtain the
C                 result at TOUT, to handle solutions at many specific
C                 values TOUT efficiently. Sometimes it is not possible
C                 to integrate beyond some  point TSTOP because the
C                 equation changes there or it is not defined past
C                 TSTOP. Then the user must tell the code not to go past.
C                 Can the integration be carried out without any
C                 restrictions on the independent variable T ...
C                 Yes-  Set INFO(4)=0
C                 No  - Set INFO(4)=1
C                       and define the stopping point TSTOP by setting
C                       RWORK(1)=TSTOP,
C
C         INFO(5) The code assumes that the user wants to solve a time
C                 dependent problem. If E(t) and A(t) are both constant
C                 the user can help the code.
C                 Are E(t) or A(t) time dependent ...
C                 Yes - Set INFO(5)=0
C                 No  - Set INFO(5)=1.
C
C         INFO(7) The user can specify a maximum (absolute value of)
C                 stepsize, so that the code will avoid passing over
C                 very large regions.
C                 Do you want the code to decide on its own maximum
C                 stepsize ...
C                 Yes - Set INFO(7)=0
C                 No  - Set INFO(7)=1
C                       and define HMAX by setting RWORK(2)=HMAX.
C
C         INFO(8) Differential/algebraic problems may occasionally
C                 suffer from severe scaling difficulties on the first
C                 step. If the user knows a great deal about the scaling
C                 of his problem, he can help to alleviate this problem
C                 by specifying an initial stepsize H0.
C                 Do you want the code to define its own initial
C                 stepsize ...
C                 Yes - Set INFO(8)=0
C                 No  - Set INFO(8)=1
C                       and define H0 by setting RWORK(3)=H0.
C
C         INFO(9) If storage is a severe problem and the user use the
C                 BDF-solver (METHOD=1), he can save some memory by
C                 restricting the maximum order MAXORD. The default
C                 value is 5. For each order decrease below 5, the code
C                 requires fewer locations, however it is likely to be
C                 slower. In any case, the user must have 1.LE.MAXORD 5.
C                 Do you want the maximum order to default to 5 ...
C                 Yes - Set INFO(9)=0
C                 No  - Set INFO(9)=1
C                       and define MAXORD by setting IWORK(3)=MAXORD.
C
C        INFO(10) The code tries to calculate the strangeness index of
C                 the problem, however it needs more memory for high
C                 index problems. The default value for the maximum
C                 index MXINDX is 3. The user can decrease it below 3 to
C                 save memory (if the index of his problem is smaller
C                 than 3) or increase it to solve a higher index
C                 problem. Note, that EDIF, ADIF and FDIF must provide
C                 E(t), A(t), F(t) and (maybe) their first MXINDX
C                 derivatives. In any case, the user must have MXINDX.GE.0.
C                 Do you want the maximum index to default to 3 ...
C                 Yes - Set INFO(10)=0
C                 No  - Set INFO(10)=1
C                       and define MAXINDX by setting IWORK(4)=MAXINDX.
C
C        INFO(11) In this code it is not neccesary to provide consistent
C                 initial conditions. Using the special structure of the
C                 strangeness free DAE, the code can compute consistent
C                 initial values to start the integration (see also
C                 INFO(12)). However, often consistent initial values
C                 are known and the code should use these values.
C                 Do you want the code to compute consistent initial
C                 values ...
C                 Yes - Set INFO(11)=0
C                 No  - Set INFO(11)=1
C                       and input consistent initial values in X.
C
C        INFO(12) If INFO(11)=0, the code computes consistent initial
C                 values in the least squares sense. The default method
C                 is to  compute consistent initial values which are
C                 close (in the least squares sense) to the given
C                 X. Sometimes the user knows which are the differential
C                 variables and he wants to prescribe these variables.
C                 In this case, the user can use a different method,
C                 which keeps the differential variables fixed.
C                 Do you want the code to use the default method for
C                 computing consistent initial values ...
C                 Yes - Set INFO(12)=0
C                 No  - Set INFO(12)=1.
C
C        INFO(14) A maximum number of steps NMAX must be specified in
C                 order to prevent the code from computing infinitely
C                 further in the case of repeated step rejection. The
C                 default value for NMAX is 10 000.
C                 Do you want the maximum number of steps to default to
C                 10 000 ...
C                 Yes - Set INFO(11)=0
C                 No  - Set INFO(11)=1
C                       and define NMAX by setting IWORK(20)=NMAX.
C
C        INFO(15) For the Runge--Kutta branch, the user can choose
C                 between two step size strategies. If not specified
C                 otherwise, the code will use the modified predictive
C                 controller of Gustafsson, which seems to produce
C                 safer results. As alternative for simple problems, the
C                 user can apply the classical step size control which
C                 produces often slightly faster runs.
C                 Do you want to use the predictive controller of
C                 Gustafsson ...
C                 Yes - Set INFO(12)=0
C                 No  - Set INFO(12)=1.
C
C        INFO(16) For the Runge--Kutta branch, a safety factor SAFE in
C                 step size prediction is used in the formula for
C                 calculating the new step size in dependency of the old
C                 one and the error norm. The smaller SAFE is chosen,
C                 the more the new step size is restricted. SAFE must
C                 lie in the interval (0.001,1). The default value is
C                 SAFE=0.9. Furthermore, parameters FACL, FACR for
C                 step size selection restrict the relation between the
C                 old and the new stepsize. The new step size is chosen
C                 subject to 1/FACL.LE.HNEW/HOLD.LE.1/FACR. The default
C                 values are FACL=5.0 and FACR=0.125.
C                 Do you want SAFE, FACL, FACR to default to $0.9, 5.0$
C                 and 0.125, respectively ...
C                 Yes - Set INFO(13)=0
C                 No  - Set INFO(13)=1
C                       and define SAFE, FACL, FACR by setting
C                       RWORK(11)=SAFE,  RWORK(12)=FACL, RWORK(13)=FACR.
C
C        INFO(17) For the Runge--Kutta branch, if HNEW is not far from
C                 HOLD (QUOT1 < HNEW/HOLD < QUOT2) and the matrices E
C                 and A are constant, work can be saved by setting
C                 HNEW=HOLD and using the system matrix of the previous
C                 step, so that a new LU--decomposition is not
C                 necessary. For small systems one may have QUOT1=1.0,
C                 QUOT2=1.2, for large full systems QUOT1=0.99,
C                 QUOT2=2.0 might be good. Default values are QUOT1=1.0,
C                 QUOT2=1.2.
C                 Do you want QUOT1, QUOT2 to default to 1.0 and 1.2,
C                 respectively ...
C                 Yes - Set INFO(14)=0
C                 No  - Set INFO(14)=1
C                      and define QUOT1, QUOT2 by setting
C                      RWORK(14)=QUOT1, RWORK(15)=QUOT2.
C
C       CONTINUING THE INTEGRATION
C
C       The user must monitor the IERR parameter in order to determine
C       what to do next.
C
C       Do not alter any quantity not specifically permitted below, in
C       particular do not alter N, T, X(*), RWORK(*), IWORK(*) or the
C       differential equation in subroutines EDIF, ADIF and FDIF. Any
C       such alteration constitutes a new problem and must be treated as
C       such, i.e., the user must start afresh. The user cannot change
C       from vector to scalar error control or vice versa (INFO(2)), but
C       he can change the size of the entries of RTOL and ATOL.
C       Increasing a tolerance makes the equation easier to integrate.
C       Decreasing a tolerance will make the equation harder to
C       integrate and should generally be avoided.
C
C       The user can switch from the intermediate-output mode to the
C       interval mode (INFO(3)) or vice versa at any time.
C
C       If it has been necessary to prevent the integration from going
C       past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
C       code will not integrate to any TOUT beyond the currently
C       specified TSTOP. Once TSTOP has been reached the user must
C       change the value of TSTOP or set INFO(4)=0. The user may change
C       INFO(4) or TSTOP at any time but he must supply the value of
C       TSTOP in RWORK(1) whenever he set INFO(4)=1.
C
C       The user should not change INFO(5), IWORK(1), or IWORK(2) unless
C       he is going to restart the code.
C
C                   *** Following a completed task ***
C
C       If
C       IERR = 1 call the code again to continue the integration
C                another step in the direction of TOUT.
C       IERR = 2 or 3 define a new TOUT and call the code again. TOUT
C                must be different from T. The user cannot change the
C                direction of integration without restarting.
C       IERR = 4 define a new TOUT and call the code again. TOUT must be
C                different from T.
C
C                   *** Following an interrupted task ***
C
C       To show the code that the user realizes the task was interrupted
C       and that he wants to continue, he must take appropriate action
C       and set INFO(1)=1.
C
C       If
C       IERR = -1  DGELDA has taken  10000 steps. If the user wants to
C                  continue, set INFO(1)=1 and call the code again. An
C                  additional 10000 steps will be allowed.
C       IERR = -2  The error tolerances RTOL and ATOL have been increased
C                  to values the code estimates appropriate for
C                  continuing. The user may want to change them himself.
C                  If the user is sure he wants to continue with relaxed
C                  error tolerances, set INFO(1)=1 and call the code
C                  again.
C       IERR = -3  A solution component is zero and the user sets the
C                  corresponding component of ATOL to zero. If the user
C                  is sure he want to continue, he must first alter the
C                  error criterion to use positive values for those
C                  components of ATOL corresponding to zero solution
C                  components, then set INFO(1)=1 and call the code
C                  again.
C       IERR = -6  Repeated error test failures occurred on the last
C                  attempted step in DGELDA. A singularity in the
C                  solution may be present. If the user is absolutely
C                  certain he wants to continue, he should restart the
C                  integration.
C       IERR = -8  The linear system solver routines failed several
C                  times. It is possible that your problem is ill-posed,
C                  and cannot be solved using this code.
C       IERR = -10 IERR was repeatedly equal to minus one in EDIF, ADIF
C                  or FDIF. If you are absolutely certain you want to
C                  continue, you should restart the integration.
C       IERR = -21 IERR=-2 was encountered, and control is being
C                  returned to the calling program.
C       IERR = -22 DGELDA could not  determine the strangeness index. It
C                  is possible that your problem is ill-posed, and cannot
C                  be solved using this code.
C       IERR = -23 DGELDA could not compute an equivalent strangeness
C                  index zero system. It is possible that your problem
C                  is ill-posed, and cannot be solved using this code.
C       IERR = -24 DGELDA could not compute an initial X. It is
C                  possible that your problem is ill-posed, and cannot
C                  be solved using this code.
C       IERR = -25 DGELDA could not compute an initial derivative. It
C                  is possible that your problem is ill-posed, and
C                  cannot be solved using this code.
C       IERR = -26 The strangeness index or the characteristic values
C                  changed on the last step. It is possible that a
C                  solution to your problem either does not exist.
C
C                   *** Following a terminated task ***
C
C        IF IERR<-100 the user cannot continue the solution of this
C        problem. An attempt to do so will result in the users run being
C        terminated.
C
C     CONTRIBUTORS
C
C     W. Rath, J. Weickert (TU Chemnitz, Germany).
C
C     REVISIONS
C
C     1995, July 15 [Version 1.1]
C       Changed order of IPAR, RPAR, and IWORK, RWORK to meet
C       SLICOT interface standard.
C       Added equilibration of linear systems for the BDF solver,
C       therefore we need additional space in RWORK (2*NEQ). (W. Rath)
C
C     1995, July 10 [Version 1]
C       First release. (W. Rath)
C
C     ******************************************************************
C
C     .. Parameters ..
C     Set pointers into CVAL
      INTEGER          OIA, OID, OIU, OM
      PARAMETER        (OM=1, OID=2, OIA=3, OIU=4)
C
C     Set pointers into IWORK
      INTEGER          OIWM, OMXORD, OMXINDX, OPHASE, OK, OKOLD,
     $                 ONS, ONSTL, ONST, ONEV, ONFA, OETF, OCTF, ONSING,
     $                 OREQU, OCEQU, OIREQ, ONMAX
      PARAMETER        (OIWM=1, OMXORD=3, OMXINDX=4, OPHASE=6,
     $                  OK=7, OKOLD=8, ONS=9, ONSTL=10, ONST=11,
     $                  ONEV=12, ONFA=13, OETF=14, OCTF=15, ONSING=16,
     $                  OREQU=17, OCEQU=18, OIREQ=19, ONMAX=20)
C
C     Set pointers into RWORK
      INTEGER          OTSTOP, OHMAX, OH, OTN, OCJ, OCJOLD, OHOLD,
     $                 OERRACC, OROUND, OSAFE, OFACL, OFACR, OQUOT1,
     $                 OQUOT2, OE
      PARAMETER        (OTSTOP=1, OHMAX=2, OH=3, OTN=4, OCJ=5, OCJOLD=6,
     $                  OHOLD=7, OERRACC=8, OROUND=9, OSAFE=11,
     $                  OFACL=12, OFACR=13, OQUOT1=14, OQUOT2=15, OE=21)
C     .. Subroutine Arguments ..
      EXTERNAL         EDIF, ADIF, FDIF
C     .. Scalar Arguments ...
      DOUBLE PRECISION T, TOUT
      INTEGER          IERR, IWARN, LIW, LRW, NEQ, METHOD
C     .. Array Arguments ...
      DOUBLE PRECISION ATOL(*), RTOL(*), RPAR(*), RWORK(*), X(*),
     $                 XPRIME(*)
      INTEGER          CVAL(4), INFO(20), IPAR(*), IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION ATOLI, FACL, FACR, H, HMAX, HMIN, HO, QUOT1,
     $                 QUOT2, R, RH, RTOLI, SAFE, TDIST, TN, TNEXT,
     $                 TSTOP, UROUND, XPNORM
      INTEGER          I, ITEMP, LDA, LDA0, LDAH, LDAQ, LDCONT, LDE,
     $                 LDE0, LDEQ, LDPHI, LDW, LDZ1, LDZ2Q, LENIW,
     $                 LENRW, LWORK, MXINDX, MXORD, N3, NMAX, NZFLG,
     $                 OA, OA0, OALPHA, OAH, OAQ, OB, OBETA, OCS,
     $                 OCONT, OD, ODELTA, OE0, OEQ, OERRV, OF, OFQ,
     $                 OGAMMA, OPHI, OPSI, ORS, OSIGMA, OV, OW,
     $                 OWORK, OWT, OYY, OZ1, OZ2Q, OZR1, OZR2, OZR3
      LOGICAL          DONE, MCONST, PRED
C     .. External Functions ..
      DOUBLE PRECISION  DGENRM, DLAMCH
      EXTERNAL          DGENRM, DLAMCH
C     .. External Subroutines ..
      EXTERNAL          DBDSTP, DBDTRP, DCOPY, DGEINI, DSCAL, DRKSTP
C     .. Executable Statements ..
      IF (INFO(1).NE.0) GOTO 100
C
C----------------------------------------------------------------------
C     This block is executed for the initial call only.
C     It contains checking of inputs and initializations.
C----------------------------------------------------------------------
C
C     First check INFO array to make sure all elements of INFO
C     are either zero or one.
      DO 10 I=2,20
         IF(INFO(I).NE.0. AND. INFO(I).NE.1) GOTO 701
10    CONTINUE
C
C     Check if NEQ is greater than zero.
      IF (NEQ.LE.0) GOTO 702
C
C     Check and compute maximum order.
      MXORD=5
      IF(INFO(9).NE.0) THEN
         MXORD=IWORK(OMXORD)
         IF (MXORD.LT.1.OR.MXORD.GT.5) GOTO 703
      ENDIF
      IWORK(OMXORD)=MXORD
C
C     Check and compute maximum index.
      MXINDX=3
      IF (INFO(10).NE.0) THEN
         MXINDX=IWORK(OMXINDX)
         IF (MXINDX.LT.0) THEN
            IERR = -104
            GOTO 704
         ENDIF
      ENDIF
      IWORK(OMXINDX)=MXINDX
C
C     Compute LENRW and LENIW.
      IF (METHOD .EQ. 1) THEN
         LENRW=50+(21+IWORK(OMXORD)+14*IWORK(OMXINDX))*NEQ
     $        + (6+IWORK(OMXINDX)+6*(IWORK(OMXINDX)+1)**2)*NEQ*NEQ
         LENIW=20+2*NEQ
      ELSE
         LENRW=20 + (15+14*(IWORK(OMXINDX)+1))*NEQ
     $        +(14+(IWORK(OMXINDX)+1)+6*(IWORK(OMXINDX)+1)**2)*NEQ*NEQ
         LENIW=20+6*NEQ
      END IF
C
C     Check lengths of RWORK and IWORK.
      IF (LRW.LT.LENRW) GOTO 705
      IF (LIW.LT.LENIW) GOTO 706
C
C     Check HMAX.
      IF (INFO(7).NE.0) THEN
         HMAX=RWORK(OHMAX)
         IF(HMAX.LE.0.0D0) GOTO 711
      ENDIF
C
C     Check NMAX , the maximal number of steps.
      IF (INFO(14).EQ.0) THEN
         NMAX=10000
      ELSE
         NMAX=IWORK(20)
         IF (NMAX.LE.0) GOTO 720
      END IF
      IWORK(ONMAX)=NMAX
C
C     Check the safety factor and parameters in step size prediction.
      IF (INFO(16).EQ.0) THEN
         SAFE=0.9D0
         FACL=5.0D0
         FACR=0.125D0
         RWORK(OSAFE)=SAFE
         RWORK(OFACL)=FACL
         RWORK(OFACR)=FACR
      ELSE
         SAFE=RWORK(OSAFE)
         IF (SAFE.LE.0.001D0.OR.SAFE.GE.1.0D0) GOTO 721
         FACL=RWORK(OFACL)
         FACR=RWORK(OFACR)
         IF (FACL.LT.1.0D0.OR.FACR.GT.1.0D0) GOTO 722
      END IF
C
C     QUOT1 and QUOT2: If QUOT1 < HNEW/HOLD < QUOT2, stepsize = const.
      IF (INFO(17).EQ.0) THEN
         QUOT1=1.0D0
         QUOT2=1.2D0
         RWORK(OQUOT1)=QUOT1
         RWORK(OQUOT2)=QUOT2
      ELSE
         QUOT1=RWORK(OQUOT1)
         QUOT2=RWORK(OQUOT2)
         IF (QUOT1.GT.1.0D0 .OR. QUOT2.LT.1.0D0) GOTO 723
      END IF

C
C     Initialize counters
      IWORK(ONST)=0
      IWORK(ONEV)=0
      IWORK(ONFA)=0
      IWORK(OETF)=0
      IWORK(OCTF)=0
      IWORK(ONSTL)=0
      IERR=1
      GO TO 200
C
C----------------------------------------------------------------------
C     This block is for continuation calls
C     only. Here we check INFO(1), and if the
C     last step was interrupted we check whether
C     appropriate action was taken.
C----------------------------------------------------------------------
C
100   CONTINUE
      IF(INFO(1).GE.1) GOTO 110
      IF(INFO(1).NE.-1) THEN
         IERR = -101
         GOTO 701
      ENDIF
C
C     If we are here, the last step was interrupted
C     by an error condition from DBDSTP, and
C     appropriate action was not taken. This
C     is a fatal error.
      PRINT *,' THE LAST STEP TERMINATED WITH A NEGATIVE VALUE OF',
     $     ' IERR =',IERR
      PRINT *,' AND NO APPROPRIATE ACTION WAS TAKEN.'
      PRINT *,' RUN TERMINATED'
      IERR = -998
      STOP
110   CONTINUE
      IWORK(ONSTL)=IWORK(ONST)
C
C----------------------------------------------------------------------
C     This block is executed on all calls.
C     The error tolerance parameters are
C     checked, and the work array pointers
C     are set.
C----------------------------------------------------------------------
C
200   CONTINUE
C
C     Check if E(t) and A(t) are constant.
      IF (INFO(5).EQ.0) THEN
        MCONST = .FALSE.
      ELSE
        MCONST = .TRUE.
      ENDIF
C
C     Check whether predictive step size control shall be used.
      IF (INFO(15).EQ.0) THEN
         PRED = .TRUE.
      ELSE
         PRED = .FALSE.
      END IF
C
C     Check  RTOL and ATOL.
      NZFLG=0
      IF(INFO(2).EQ.0) THEN
        RTOLI=RTOL(1)
        ATOLI=ATOL(1)
	IF(RTOLI.GT.0.0D0.OR.ATOLI.GT.0.0D0) NZFLG=1
        IF(RTOLI.LT.0.0D0) GOTO 707
        IF(ATOLI.LT.0.0D0) GOTO 708
      ELSE
        DO 210 I=1,NEQ
          RTOLI=RTOL(I)
          ATOLI=ATOL(I)
          IF(RTOLI.GT.0.0D0.OR.ATOLI.GT.0.0D0) NZFLG=1
          IF(RTOLI.LT.0.0D0) GOTO 707
          IF(ATOLI.LT.0.0D0) GOTO 708
210     CONTINUE
      ENDIF
      IF (NZFLG.EQ.0) GOTO 709
C
C     Set up leading dimensions for RWORK matricies.
      IF (METHOD .EQ. 1) THEN
         LDPHI   = NEQ
         LDW     = NEQ
      ELSE
         N3      = 3*NEQ
         LDCONT  = NEQ
         LDW     = N3
         LDE0    = NEQ
         LDA0    = NEQ
      END IF
      LDE     = NEQ
      LDA     = NEQ
      LDZ1    = NEQ
      LDZ2Q   = (IWORK(OMXINDX)+1)*NEQ
      LDEQ    = (IWORK(OMXINDX)+1)*NEQ
      LDAQ    = (IWORK(OMXINDX)+1)*NEQ
      LDAH    = (IWORK(OMXINDX)+1)*NEQ
C
C     Set up RWORK storage. IWORK is fixed in DATA SEGMENT.
      OA      = OE + LDE*LDE
      OF      = OA + LDA*LDA
      OEQ     = OF + NEQ
      OAQ     = OEQ + LDEQ*LDEQ
      OFQ     = OAQ + LDAQ*NEQ
      OZ1     = OFQ + (IWORK(OMXINDX)+1)*NEQ
      OZ2Q    = OZ1 + LDZ1*LDZ1
      OAH     = OZ2Q + LDZ2Q*LDZ2Q
      OD      = OAH + LDAH*LDAH
      OV      = OD + (IWORK(OMXINDX)+1)*NEQ
      OWT     = OV +  (IWORK(OMXINDX)+1)*NEQ
      IF (METHOD .EQ. 1) THEN
         OALPHA  = OWT + NEQ
         OBETA   = OALPHA + 6
         OGAMMA  = OBETA + 6
         OPSI    = OGAMMA + 6
         OSIGMA  = OPSI + 6
         ODELTA  = OSIGMA + 6
         OERRV   = ODELTA + NEQ
         OPHI    = OERRV + NEQ
         ORS     = OPHI + (IWORK(OMXORD)+1)*LDPHI
         OCS     = ORS + NEQ
         OW      = OCS + NEQ
      ELSE
         OZR1    = OWT + NEQ
         OZR2    = OZR1 + NEQ
         OZR3    = OZR2 + NEQ
         OYY     = OZR3 + NEQ
         OB      = OYY + N3
         OCONT   = OB + N3
         OE0     = OCONT + LDCONT*4
         OA0     = OE0 + LDE0**2
         OW      = OA0 + LDA0**2
      END IF
      OWORK   = OW + LDW*LDW
      LWORK   = LRW - OWORK
C
      IF(INFO(1).EQ.1) GOTO 400
C
C----------------------------------------------------------------------
C     This block is executed on the initial call
C     only.
C     Set the initial step size, and
C     the error weight vector, and PHI.
C     Compute initial XPRIME and X, if necessary.
C----------------------------------------------------------------------
C
      TN=T
      IERR=1
C
C     Compute unit roundoff.
      UROUND = DLAMCH('P')
      RWORK(OROUND) = UROUND
      IF (INFO(1).EQ.0) THEN
C
C     Computel characteristic values and consistent initial condidions if
C     necessary.
         CALL DGEINI(EDIF, ADIF, FDIF, NEQ, TN, IWORK(OMXINDX), X,
     $               XPRIME, CVAL(OM), CVAL(OID), CVAL(OIA),
     $               CVAL(OIU), IWORK(OIREQ), RWORK(OEQ), LDEQ,
     $               RWORK(OAQ), LDAQ, RWORK(OFQ), RWORK(OE), LDE,
     $               RWORK(OA), LDA, RWORK(OF), RWORK(OZ1), LDZ1,
     $               RWORK(OZ2Q), LDZ2Q, RWORK(OD), RWORK(OAH), LDAH,
     $               RWORK(OV), IPAR, RPAR, RWORK(OWORK), LWORK,
     $               INFO(11),  IERR)
         IF (IERR.LT.0) THEN
            GOTO 340
         ENDIF
         IF (TOUT .EQ. T) THEN
            INFO(1) = 2
            IERR = 4
            GOTO 580
         ENDIF
      ENDIF
C
C     Check to see that TOUT is different from T.
      IF (TOUT .EQ. T) GOTO 719
C
C     Compute HMIN.
      HMIN = 4.0D0*UROUND*MAX(ABS(T),ABS(TOUT))
C
C     Check initial interval to see that it is long enough.
      TDIST = ABS(TOUT - T)
      IF (TDIST .LT. HMIN) GOTO 715
C
C     Set error weight vector WT.
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 310 I=1,NEQ
        IF (INFO(2) .NE.0) THEN
          RTOLI=RTOL(I)
          ATOLI=ATOL(I)
        ENDIF
	RWORK(OWT+I-1)=RTOLI*ABS(X(I))+ATOLI
        IF (RWORK(OWT+I-1).LE.0.0D0) GOTO 714
 310  CONTINUE
C
C     Invert the WT vector to minimize number of division operations.
      DO 320 I = 1,NEQ
         RWORK(OWT+I-1) = 1.D0/RWORK(OWT+I-1)
 320  CONTINUE
C
C     Check HO, if this was input.
      IF (INFO(8) .NE. 0) THEN
         HO = RWORK(OH)
         IF ((TOUT - T)*HO .LT. 0.0D0) GOTO 712
         IF (HO .EQ. 0.0D0) GOTO 713
      ELSE
C
C     Compute initial stepsize, to be used by DBDSTP and DRKSTP.
        HO = 0.001D0*TDIST
        XPNORM = DGENRM(NEQ,XPRIME,RWORK(OWT))
        IF (XPNORM .GT. 0.5D0/HO) HO = 0.5D0/XPNORM
        HO = SIGN(HO,TOUT-T)
      ENDIF
C
C     Adjust HO if necessary to meet HMAX bound.
      IF (INFO(7) .NE. 0) THEN
         RH = ABS(HO)/RWORK(OHMAX)
         IF (RH .GT. 1.0D0) HO = HO/RH
      ENDIF
C
C     Compute TSTOP, if applicable.
      IF (INFO(4) .NE. 0) THEN
         TSTOP = RWORK(OTSTOP)
         IF ((TSTOP - T)*HO .LT. 0.0D0) GOTO 716
         IF ((T + HO - TSTOP)*HO .GT. 0.0D0) HO = TSTOP - T
         IF ((TSTOP - TOUT)*HO .LT. 0.0D0) GOTO 710
      ENDIF
C
C     Load H with HO. store H in RWORK(OH).
      H = HO
      RWORK(OH) = H
C
C     Load X and H*XPRIME into PHI(*,1) and PHI(*,2).
      IF (METHOD .EQ. 1) THEN
         CALL DCOPY(NEQ,X,1,RWORK(OPHI),1)
         ITEMP = OPHI + NEQ
         DO 330 I = 1,NEQ
            RWORK(ITEMP + I - 1) = H*XPRIME(I)
 330     CONTINUE
      END IF
C
 340  IF (IERR.LT.0) GOTO 600
      GOTO 500
C
C-------------------------------------------------------
C     This block is for continuation calls only. Its
C     purpose is to check stop conditions before
C     taking a step.
C     Adjust H if necessary to meet HMAX bound.
C-------------------------------------------------------
C
400   CONTINUE
      UROUND=RWORK(OROUND)
      DONE = .FALSE.
      TN=RWORK(OTN)
      H=RWORK(OH)
      IF (INFO(7) .NE. 0) THEN
        RH = ABS(H)/RWORK(OHMAX)
        IF (RH .GT. 1.0D0) H = H/RH
      ENDIF
      IF(T .EQ. TOUT) GOTO 719
      IF((T - TOUT)*H .GT. 0.0D0) GOTO 712
      IF(INFO(4) .EQ. 1) THEN
         TSTOP=RWORK(OTSTOP)
         IF ((TN-TSTOP)*H.GT.0.0D0) GOTO 716
         IF ((TSTOP-TOUT)*H.LT.0.0D0) GOTO 710
      END IF
      IF (METHOD .EQ. 1) THEN
C
C        Backward-differentiation branch.
         IF (INFO(4) .NE. 1) THEN
           IF (INFO(3) .NE. 1) THEN
             IF ((TN-TOUT)*H.LT.0.0D0) GOTO 420
             CALL DBDTRP(TN,TOUT,NEQ,IWORK(OKOLD),
     $            RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
             T=TOUT
             IERR = 3
           ELSE
             IF ((TN-T)*H .LE. 0.0D0) GOTO 420
             IF ((TN - TOUT)*H .LE. 0.0D0) THEN
                CALL DBDTRP(TN,TN,NEQ,IWORK(OKOLD),
     $               RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
                T = TN
                IERR = 1
             ELSE
                CALL DBDTRP(TN,TOUT,NEQ,IWORK(OKOLD),
     $               RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
                T = TOUT
                IERR = 3
             ENDIF
           ENDIF
         ELSEIF (INFO(3) .NE. 1) THEN
           IF ((TN-TOUT)*H.LT.0.0D0) GOTO 410
           CALL DBDTRP(TN,TOUT,NEQ,IWORK(OKOLD),
     $            RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
           T=TOUT
           IERR = 3
         ELSE
           IF ((TN-T)*H .LE. 0.0D0) GOTO 410
           IF ((TN - TOUT)*H .LE. 0.0D0) THEN
              CALL DBDTRP(TN,TN,NEQ,IWORK(OKOLD),
     $             RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
              T = TN
              IERR = 1
           ELSE
              CALL DBDTRP(TN,TOUT,NEQ,IWORK(OKOLD),
     $             RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
              T = TOUT
              IERR = 3
           ENDIF
         ENDIF
         DONE = .TRUE.
         GOTO 420
 410     CONTINUE
C
C        Check whether we are within roundoff of tstop.
         IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*
     $      (ABS(TN)+ABS(H))) THEN
            CALL DBDTRP(TN,TSTOP,NEQ,IWORK(OKOLD),
     $           RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
           IERR=2
           T=TSTOP
           DONE = .TRUE.
           GOTO 420
         ELSE
 460        TNEXT=TN+H
           IF ((TNEXT-TSTOP)*H.LE.0.0D0) GOTO 420
           H=TSTOP-TN
           RWORK(OH)=H
         ENDIF
      END IF
C
420   IF (DONE) GO TO 580
C
C-------------------------------------------------------
C     The next block contains the call to the
C     one-step integrators DBDSTP and DRKSTP.
C     This is a looping point for the integration steps.
C     Check for too many steps.
C     Update WT.
C     Check for too much accuracy requested.
C     Compute minimum stepsize.
C-------------------------------------------------------
C
 500  CONTINUE
      DONE=.FALSE.
C
C     Check for too many steps.
      IF ((IWORK(ONST)-IWORK(ONSTL)).GE.IWORK(ONMAX)) THEN
        IERR=-1
        GOTO 540
      ENDIF
C
C     Update WT.
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 520 I=1,NEQ
         IF (INFO(2) .NE.0) THEN
            RTOLI=RTOL(I)
            ATOLI=ATOL(I)
         ENDIF
         RWORK(OWT+I-1)=RTOLI*ABS(X(I))+ATOLI
         IF (RWORK(I+OWT-1).LE.0.0D0) THEN
            IERR=-3
            GOTO 540
         ENDIF
 520  CONTINUE
C
C     Invert the WT vector to minimize number of division operations.
      DO 530 I = 1,NEQ
         RWORK(OWT+I-1) = 1.D0/RWORK(OWT+I-1)
 530  CONTINUE
C
C     Test for too much accuracy requested.
      R=DGENRM(NEQ,X,RWORK(OWT))*100.0D0*UROUND
      IF(R.GT.1.0D0) THEN
C
C       Multiply RTOL and ATOL by R and return.
        IF(INFO(2).NE.1) THEN
           RTOL(1)=R*RTOL(1)
           ATOL(1)=R*ATOL(1)
	ELSE
	  CALL DSCAL(NEQ,R,RTOL,1)
	  CALL DSCAL(NEQ,R,ATOL,1)
	ENDIF
	IERR=-2
	GOTO 540
      ENDIF
C
C     Compute minimum stepsize.
      HMIN=4.0D0*UROUND*MAX(ABS(TN),ABS(TOUT))
C
C     Test H vs. HMAX.
      IF (INFO(7) .NE. 0) THEN
         RH = ABS(H)/RWORK(OHMAX)
         IF (RH .GT. 1.0D0) H = H/RH
      ENDIF
C
C     Runge-kutta integrator should not step past TOUT.
      IF (METHOD .EQ. 2 .AND.
     $    (TN+H*1.0001D0-TOUT)*H .GE. 0.0D0) THEN
         HO=H
         H=TOUT-TN
         DONE=.TRUE.
      END IF
C
      IF (METHOD .EQ. 1) THEN
         CALL DBDSTP(EDIF, ADIF, FDIF, NEQ, TN, H, RWORK(OHOLD), HMIN,
     $               RWORK(OCJ),RWORK(OCJOLD), IWORK(OK), IWORK(OKOLD),
     $               CVAL(OM), CVAL(OID), CVAL(OIA), CVAL(OIU),
     $               IWORK(OIREQ), IWORK(OMXINDX), IWORK(ONS),
     $               IWORK(OPHASE), IWORK(OREQU), IWORK(OCEQU), INFO(1),
     $               X, XPRIME, RWORK(OE), LDE,
     $               RWORK(OA), LDA, RWORK(OZ1), LDZ1, RWORK(OZ2Q),
     $               LDZ2Q, RWORK(OWT), RWORK(OW), LDW,
     $               RWORK(OPHI), LDPHI, RWORK(OALPHA),
     $               RWORK(OBETA), RWORK(OGAMMA), RWORK(OPSI),
     $               RWORK(OSIGMA), RWORK(ODELTA), RWORK(OERRV),
     $               RWORK(ORS), RWORK(OCS), RWORK(OF), RWORK(OEQ),
     $               LDEQ, RWORK(OAQ), LDAQ, RWORK(OFQ), RWORK(OAH),
     $               LDAH, IPAR, RPAR, IWORK(OIWM), RWORK(OWORK), LWORK,
     $               MCONST, IWARN, IERR)
      ELSE
         CALL DRKSTP(EDIF, ADIF, FDIF, NEQ, N3, TN, H, RWORK(OHOLD),
     $               HMIN, RWORK(OERRACC), RWORK(OROUND), RWORK(OSAFE),
     $               RWORK(OFACL), RWORK(OFACR), RWORK(OQUOT1),
     $               RWORK(OQUOT2), CVAL(OM), CVAL(OID), CVAL(OIA),
     $               CVAL(OIU), IWORK(OIREQ), IWORK(OMXINDX),
     $               IWORK(ONSING), INFO(1), DONE, X, XPRIME,
     $               RWORK(OE0), LDE0, RWORK(OA0), LDA0, RWORK(OE), LDE,
     $               RWORK(OA), LDA, RWORK(OEQ), LDEQ, RWORK(OAQ), LDAQ,
     $               RWORK(OZ1), LDZ1, RWORK(OZ2Q), LDZ2Q, RWORK(OAH),
     $               LDAH, RWORK(OW), LDW, RWORK(OCONT), LDCONT,
     $               RWORK(OWT), RWORK(OZR1), RWORK(OZR2), RWORK(OZR3),
     $               RWORK(OYY), RWORK(OB), RWORK(OF), RWORK(OFQ), IPAR,
     $               RPAR, IWORK(OIWM), RWORK(OWORK), LWORK, PRED,
     $               MCONST, IWARN, IERR)
      ENDIF
540   IF(IERR.LT.0)GO TO 600
C
C--------------------------------------------------------
C     This block handles the case of a successful return
C     from DBDSTP or DRKSTP (IERR=1).
C     Test for stop conditions.
C--------------------------------------------------------
C
C     Runge-Kutta integrator reached exactly TOUT.
      IF (DONE) THEN
         H=HO
         T=TN
         IERR=3
         GO TO 580
      END IF
C
C     Backward-differentiation branch.
      IF (METHOD .EQ. 1) THEN
         IF (INFO(4).EQ.0) THEN
           IF (INFO(3).EQ.0) THEN
             IF ((TN-TOUT)*H.LT.0.0D0)GO TO 500
             CALL DBDTRP(TN,TOUT,NEQ,IWORK(OKOLD),
     $             RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
             IERR=3
             T=TOUT
           ELSEIF ((TN-TOUT)*H.LT.0.0D0) THEN
             T=TN
             IERR=1
           ELSE
             CALL DBDTRP(TN,TOUT,NEQ,IWORK(OKOLD),
     $             RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
             IERR=3
             T=TOUT
	   ENDIF
         ELSEIF (INFO(3).EQ.0) THEN
           IF ((TN-TOUT)*H.GE.0.0D0) THEN
             CALL DBDTRP(TN,TOUT,NEQ,IWORK(OKOLD),
     $             RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
             T=TOUT
             IERR=3
    	   ELSEIF (ABS(TN-TSTOP).GT.100.0D0*UROUND*
     $        (ABS(TN)+ABS(H))) THEN
             TNEXT=TN+H
             IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 500
             H=TSTOP-TN
             GO TO 500
	   ELSE
             CALL DBDTRP(TN,TSTOP,NEQ,IWORK(OKOLD),
     $             RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
             IERR=2
             T=TSTOP
   	   ENDIF
         ELSEIF ((TN-TOUT)*H.LT.0.0D0) THEN
           IF (ABS(TN-TSTOP).GT.100.0D0*UROUND*(ABS(TN)+ABS(H))) THEN
             T=TN
             IERR=1
	   ELSE
             CALL DBDTRP(TN,TSTOP,NEQ,IWORK(OKOLD),
     $             RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
             IERR=2
             T=TSTOP
	   ENDIF
         ELSE
           CALL DBDTRP(TN,TOUT,NEQ,IWORK(OKOLD),
     $           RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME)
           T=TOUT
           IERR=3
         ENDIF
C
C     Runge-Kutta branch.
      ELSE
         IF (INFO(3) .EQ. 0) THEN
            GOTO 500
         ELSE
            T=TN
            IERR=1
         END IF
      END IF
C
C--------------------------------------------------------
C     All successful returns from DGELDA are made from
C     this block.
C--------------------------------------------------------
C
580   CONTINUE
      RWORK(OTN)=TN
      RWORK(OH)=H
      RETURN
C
C----------------------------------------------------------------------
C     This block handles all unsuccessful
C     returns other than for illegal input.
C----------------------------------------------------------------------
C
600   CONTINUE
      INFO(1)=-1
      T=TN
      RWORK(OTN)=TN
      RWORK(OH)=H
      RETURN
C
C----------------------------------------------------------------------
C     This block handles all error returns due
C     to illegal inputl, as detected before calling
C     DBDSTP or DRKSTP.
C     If this happens twice in succession,
C     execution is terminated.
C----------------------------------------------------------------------
C
701   CONTINUE
      IERR = -101
      GO TO 750
C
 702  CONTINUE
      IERR = -102
      GO TO 750
C
703   CONTINUE
      IERR = -103
      GO TO 750
C
 704  CONTINUE
      IERR = -104
      GO TO 750
C
 705  CONTINUE
      IWORK(LIW) = LENRW
      IERR = -105
      GO TO 750
C
 706  CONTINUE
      IWORK(LIW) = LENIW
      IERR = -106
      GO TO 750
C
 707  CONTINUE
      IERR = -107
      GO TO 750
C
 708  CONTINUE
      IERR = -108
      GO TO 750
C
 709  CONTINUE
      IERR = -109
      GO TO 750
C
 710  CONTINUE
      IERR = -110
      GO TO 750
C
 711  CONTINUE
      IERR = -111
      GO TO 750
C
 712  CONTINUE
      IERR = -112
      GO TO 750
C
 713  CONTINUE
      IERR = -113
      GO TO 750
C
 714  CONTINUE
      IERR = -114
      GO TO 750
C
 715  CONTINUE
      IERR = -115
      GO TO 750
C
 716  CONTINUE
      IERR = -116
      GO TO 750
C
 719  CONTINUE
      IERR = -119
      GO TO 750
C
 720  CONTINUE
      IERR = -120
      GO TO 750
C
 721  CONTINUE
      IERR = -121
      GO TO 750
C
 722  CONTINUE
      IERR = -122
      GO TO 750
C
 723  CONTINUE
      IERR = -123
      GO TO 750
C
 750  CONTINUE
      IF(INFO(1).EQ.-1) THEN
         IERR = -999
      ENDIF
C
      INFO(1)=-1
      RETURN
C *** Last line of DGELDA ***
      END
