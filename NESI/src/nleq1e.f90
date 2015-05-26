module nleq_driver
use nleq_mod

contains

!-----------------------------------------------------------------
! This is just to provide an argument for the call to nleq1
!-----------------------------------------------------------------
subroutine dummy
end subroutine

!-----------------------------------------------------------------
!
!*  Title
!
!     Numerical solution of nonlinear (NL) equations (EQ)
!     especially designed for numerically sensitive problems.
!     (E)asy-to-use driver routine for NLEQ1.
!
!*  Written by        U. Nowak, L. Weimann
!*  Purpose           Solution of systems of highly nonlinear equations
!*  Method            Damped affine invariant Newton method
!                     (see references below)
!*  Category          F2a. - Systems of nonlinear equations
!*  Keywords          Nonlinear equations, Newton methods
!*  Version           2.3
!*  Revision          November 1991
!*  Latest Change     November 1991
!*  Library           CodeLib
!*  Code              Fortran 77, Double Precision
!*  Environment       Standard Fortran 77 environment on PC's,
!                     workstations and hosts.
!*  Copyright     (c) Konrad-Zuse-Zentrum fuer
!                     Informationstechnik Berlin (ZIB)
!                     Takustrasse 7, D-14195 Berlin-Dahlem
!                     phone : + 49/30/84185-0
!                     fax   : + 49/30/84185-125
!*  Contact           Lutz Weimann
!                     ZIB, Division Scientific Computing,
!                          Department Scientific Software
!                     phone : + 49/30/84185-185
!                     fax   : + 49/30/84185-107
!                     e-mail: weimann@zib.de
!
!*    References:
!
!     /1/ P. Deuflhard:
!         Newton Techniques for Highly Nonlinear Problems -
!         Theory and Algorithms.
!         Academic press Inc. (To be published)
!
!     /2/ U. Nowak, L. Weimann:
!         A Family of Newton Codes for Systems of Highly Nonlinear
!         Equations - Algorithm, Implementation, Application.
!         ZIB, Technical Report TR 90-10 (December 1990)
!
!  ---------------------------------------------------------------
!
!* Licence
!    You may use or modify this code for your own non commercial
!    purposes for an unlimited time.
!    In any case you should not deliver this code without a special
!    permission of ZIB.
!    In case you intend to use the code commercially, we oblige you
!    to sign an according licence agreement with ZIB.
!
!* Warranty
!    This code has been tested up to a certain level. Defects and
!    weaknesses, which may be included in the code, do not establish
!    any warranties by ZIB. ZIB does not take over any liabilities
!    which may follow from aquisition or application of this code.
!
!* Software status
!    This code is under care of ZIB and belongs to ZIB software class 1.
!
!     ------------------------------------------------------------
!
!*    Summary:
!     ========
!     Damped Newton-algorithm for systems of highly nonlinear
!     equations - damping strategy due to Ref. (1).
!
!     (The iteration is done by subroutine N1INT actually. NLEQ1E
!      calls the standard interface driver NLEQ1, which itself does
!      some house keeping and builds up workspace.)
!
!     Jacobian approximation by numerical differences.
!
!     The numerical solution of the arising linear equations is
!     done by means of the subroutines *GEFA and *GESL ( Gauss-
!     algorithm with column-pivoting and row-interchange;
!     replace '*' by 'S' or 'D' for single or double precision
!     version respectively).
!
!     ------------------------------------------------------------
!
!*    External subroutine to be supplied by the user
!     ==============================================
!
!     FCN(N,X,F,IFAIL) Ext    Function subroutine - the problem function
!                             (The routine must be named exactly FCN !)
!       N              Int    Number of vector components (input)
!                             Must not be altered!
!       X(N)           Dble   Vector of unknowns (input)
!                             Must not be altered!
!       F(N)           Dble   Vector of function values (output)
!       IFAIL          Int    FCN evaluation-failure indicator. (output)
!                             On input:  Has always value 0 (zero).
!                             On output: Indicates failure of FCN eval-
!                                uation, if having a nonzero value.
!                             If <0: NLEQ1E will be terminated with
!                                    IFAIL returned via IERR.
!                             If =1: A new trial Newton iterate will be
!                                    computed, with the damping factor
!                                    reduced to it's half.
!                             If =2: A new trial Newton iterate will
!                                    computed, with the damping factor
!                                    reduced by a reduction factor,
!                                    which must be output through F(1)
!                                    by FCN, and it's value must be
!                                    >0 and < 1.
!
!*    Parameters list description
!     ===========================
!
!     Name   Type    In/Out Meaning
!
!     N      Int     In     Number of unknowns ( N .LE. 50 )
!     X(N)   Dble    In     Initial estimate of the solution
!                    Out    Solution values ( or final values,
!                           respectively )
!     RTOL   Dble    In     Required relative precision of
!                           solution components -
!                           RTOL.GE.EPMACH*TEN*N
!                    Out    Finally achieved accuracy
!     IERR   Int     Out    Return value parameter
!                           < 0 Termination forced by user function FCN
!                               due to IFAIL<0 on output, IERR is set
!                               to IFAIL
!                           = 0 successfull completion of the iteration,
!                               solution has been computed
!                           > 0 see list of error/warning messages below
!
!*   Error and warning messages:
!    ===========================
!
!      1    Termination, since Jacobian matrix became singular
!      2    Termination after 100 iterations
!      3    Termination, since damping factor became too small
!      4    Warning: Superlinear or quadratic convergence slowed down
!           near the solution.
!           Iteration has been stopped therefore with an approximation
!           of the solution not such accurate as requested by RTOL,
!           because possibly the RTOL requirement may be too stringent
!           (i.e. the nonlinear problem is ill-conditioned)
!      5    Warning: Iteration stopped with termination criterion
!           (using RTOL as requested precision) satisfied, but no
!           superlinear or quadratic convergence has been indicated yet.
!           Therefore, possibly the error estimate for the solution may
!           not match good enough the really achieved accuracy.
!     20    Bad input value to parameter N; 1 .LE. N .LE. 50 required
!     21    Nonpositive value for RTOL supplied
!     82    Termination, because user routine FCN returned with IFAIL>0
!
!     Note   : in case of failure:
!        -    use better initial guess
!        -    or refine model
!        -    or use non-standard options and/or analytical Jacobian
!             via the standard interface NLEQ1
!
!*    Machine dependent constants used:
!     =================================
!
!     DOUBLE PRECISION EPMACH  in  N1PCHK, N1INT
!     DOUBLE PRECISION GREAT   in  N1PCHK
!     DOUBLE PRECISION SMALL   in  N1PCHK, N1INT, N1SCAL
!
!*    Subroutines called: NLEQ1
!
!-----------------------------------------------------------------
SUBROUTINE NLEQ1E(N,X,RTOL,nfout,IERR)
INTEGER N, nfout
DOUBLE PRECISION X(N), RTOL
INTEGER IERR
INTEGER NMAX, LIOPT, LIWK, LRWK, LUPRT
PARAMETER ( NMAX=50 )
PARAMETER ( LIOPT=50, LIWK=NMAX+50, LRWK=(NMAX+13)*NMAX+60 )
PARAMETER ( LUPRT=6 )
EXTERNAL FCN
DOUBLE PRECISION XSCAL(NMAX)
INTEGER IOPT(LIOPT), IWK(LIWK)
DOUBLE PRECISION RWK(LRWK)
INTEGER I, NIW, NRW
CHARACTER CHGDAT*20, PRODCT*8

!
!     Version: 2.3               Latest change:
!     -----------------------------------------
!
DATA      CHGDAT      /'November 15, 1991   '/
DATA      PRODCT      /'NLEQ1E  '/
!*    Begin
NIW = N+50
NRW = (N+13)*N+60
!     Checking dimensional parameter N
IF ( N.LT.1 .OR. N.GT.NMAX ) THEN
    WRITE(*,1001) NMAX, N
    IERR = 20
    RETURN
ENDIF
1001    FORMAT(/,' Error: Bad input to parameter N supplied',/,8X,'choose 1 .LE. N .LE. ',I3,' , your input is: N = ',I5)
DO I=1,LIOPT
    IOPT(I) = 0
enddo
DO I=1,NIW
    IWK(I) = 0
enddo
DO I=1,NRW
    RWK(I) = 0.0D0
enddo
DO I=1,N
    XSCAL(I) = 0.0D0
enddo
!     Print errors, warnings, monitor and time monitor
!     to standard output
IOPT(11) = 3
IOPT(12) = nfout	!LUPRT
IOPT(13) = 3
IOPT(14) = nfout	!LUPRT
IOPT(19) = 1
IOPT(20) = nfout	!LUPRT
!     Maximum number of Newton iterations
IWK(31) = 100
CALL NLEQ1(N,FCN,DUMMY,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
IF (IERR.EQ.82 .AND. IWK(23).LT.0) IERR=IWK(23)
!     End of subroutine NLEQ1E
RETURN
END subroutine

end module
