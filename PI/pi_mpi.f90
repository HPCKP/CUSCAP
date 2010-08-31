! compute pi by calculating 4*arctan(1)
! arctan(1) = Integral of (1+x^2)^-1 from 0 to 1

program pi_mpi

! include this module to use the OpenMP run-time library
! conditional OpenMP compilation in f90 is designated with !$

!$ use omp_lib

implicit none

! include this header file to access MPI parameters

include 'mpif.h'

integer(4), parameter :: n=10000000  !number of intervals between 0 and 1
real(8), parameter :: h=1/real(n)  !width of each interval

real(8) :: mypi,sumpi,sum,x
integer(4) ::myrank,numprocs,i,ierr

! this call is necessary - starts MPI communications between processes

call mpi_init(ierr)

! find out the "rank" of each process

call mpi_comm_rank(mpi_comm_world,myrank,ierr)

! find out how many processes there are in the "communicator" mpi_comm_world 

call mpi_comm_size(mpi_comm_world,numprocs,ierr)

sum=0.0

! only the "root" process returns diagnostic information

if (myrank==0) then 
  print *, 'Calculating PI with',numprocs,'processes'
endif

! OpenMP:
! spawn threads with the following data scoping
! shared -- all threads see same memory address
! private -- each thread has an automatic copy of the variable
! reduction -- threads get a private copy that is reduced
!              at the end of the loop and returned to the parent process

! MPI:
! each process computes a fraction of the iterations based on it's rank

!$omp parallel default(shared) private(i,x) reduction(+:sum)
!$ if (myrank==0.and.omp_get_thread_num()==0) print *,omp_get_num_threads(),'threads used per process'
!$ print *,'thread',omp_get_thread_num(),'of',omp_get_num_threads(),'started on rank',myrank
!$omp do
do i=myrank+1,n,numprocs
  x=h*(dble(i)-0.5)
  sum=sum+4.0/(1.0+x**2)
enddo
!$omp end do
!$omp end parallel

mypi=h*sum

! each of the MPI processes calls this MPI routine to reduce the sum 
! after this call, rank 0 will contain the reduced value in sumpi

call MPI_REDUCE(mypi,sumpi,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

if (myrank==0) print *, sumpi

! this is called to stop MPI communications

call mpi_finalize(ierr)
end program pi_mpi
