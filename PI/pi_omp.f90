! compute pi by calculating 4*arctan(1)
! arctan(1) = Integral of (1+x^2)^-1 from 0 to 1

program pi

! include this module to use the OpenMP run-time library
! conditional OpenMP compilation in f90 is designated with !$ 

!$ use omp_lib

implicit none
real(8) :: pi_int,x
integer(4), parameter :: n=100000000  !number of intervals between 0 and 1
real(8), parameter :: h=1/real(n)  !width of each interval
integer(4) :: i

pi_int=0.0

! spawn threads with the following data scoping
! shared -- all threads see same memory address
! private -- each thread has an automatic copy of the variable
! reduction -- threads get a private copy that is reduced 
!              at the end of the loop and returned to the parent process

!$omp parallel default(shared) private(i,x) reduction(+:pi_int)
!$ print *,'thread',omp_get_thread_num(),'started'
!$omp do
do i=1,n
  x=h*(dble(i)-0.5)    !calculate at center of interval
  pi_int=pi_int+4.0/(1.0+x**2)
enddo
!$omp end do
!$omp end parallel

print *,'pi =', h*pi_int 

end program pi
