module deform

use global
use nleq_driver

implicit none

contains

!---------------------------------------------------------------------------------
! To determine the parameters of a squashed and truncated sphere, conserving volume.
! The scaling is a > 1 in the xy plane, and b < 1 in the z axis.
! The bottom cR is removed.
! The squashed sphere is parametrised by (r,theta):
! r = aRsin(theta)
! z = bR(1-cos(theta))
! The angle theta0 corresponding to the slice at z=cR is given by:
! bR(1-cos(theta0)) = cR
! ==> cos(theta0) = (b-c)/b
! For the volume of the squashed sphere with the bottom cR sliced off to equal
! the volume of the original sphere:
! a^2.b(q - (1/3).q^3 + 2/3) = 4/3  where q = (b-c)/b
!
! The parameters we use to define the squashing are:
!    alpha_shape = (contact diameter)/diameter = aRsin(theta0)/(aR) = sin(theta0) = (1-q^2)^(1/2)
! and
!    beta_shape = height/diameter = (2b-c)R/(2aR) = (2b-c)/a
! 
! Therefore, for given (alpha, beta), the three simultaneous non-linear equations to solve are:
!	a^2.b(q - (1/3).q^3 + 2/3) - 4/3 = 0
!	(2b-c)/a - beta_shape = 0
!	(1-q^2)^(1/2) - alpha_shape = 0
!
! We set x(1) = a, x(2) = b, x(3) = c.
!---------------------------------------------------------------------------------
subroutine squasher
implicit none
!use deform
!use nleq_driver
INTEGER NN
PARAMETER ( NN=3 )
INTEGER N,I,N1
DOUBLE PRECISION RTOL
INTEGER IERR
DOUBLE PRECISION X(NN)
REAL STIME,ETIME,CPTIME
double precision F(NN)
     
N = 3
RTOL = 1.0D-6
! Initial guess
X(1) = 1.1
X(2) = 0.9
X(3) = 0.05

CALL NLEQ1E(N,X,RTOL,nfout,IERR)
write(*,*) 'Returned: ierr: ', ierr
write(*,'(3f12.6)') (X(i),i=1,N)	    
write(*,*) 'Check:'
call FCN(N,X,F,IERR)
write(*,'(3f12.6)') (F(i),i=1,N)
adrop = x(1)
bdrop = x(2)
cdrop = x(3)    
end subroutine
   
end module

!---------------------------------------------------------------------------------
! x(:) = (a,b,c)
!---------------------------------------------------------------------------------
subroutine FCN(N,X,F,IFAIL)
use deform
implicit none
INTEGER N, IFAIL
DOUBLE PRECISION X(N),F(N)
integer i, is
DOUBLE PRECISION q

q = (x(2) - x(3))/x(2)
F(1) = x(1)*x(1)*x(2)*(q - (q**3)/3 + 2./3.) - 4./3.
F(2) = (2*x(2) - x(3))/(2*x(1)) - beta_shape
F(3) = (1 - q*q) - alpha_shape*alpha_shape
!write(*,'(3f10.4)') x
!write(*,'(3f10.4)') F
!stop
IFAIL = 0
end subroutine
