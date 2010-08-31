// compute pi by calculating 4*arctan(1)
// arctan(1) = Integral of (1+x^2)^-1 from 0 to 1
//
// the number of threads can be set at runtime using the OMP_NUM_THREADS
// environment variable

#include <stdio.h>
#include <math.h>

// have to include omp.h to use OpenMP functions

#include <omp.h>

// have to include mpi.h to define MPI parameters

#include <mpi.h>

int main(int argc, char *argv[])
{

  double mypi,sumpi,h,sum,x;
  int n,myrank,numprocs,i;

// this call is necessary - starts MPI communications between processes

  MPI_Init(&argc, &argv);

// find out the "rank" assigned to this process

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

// find out how many processes there are in the "communicator" mpi_comm_world

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  n=100000;
  h=1.0/n;
  sum=0.0;

// root process returns diagnostic information
//
// check to see how many threads we will launch. since we only want this to compile in 
// the parallel version, we use the _OPENMP macro name

  if (myrank == 0) {
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
    printf("Calculating PI with %d processes and %d threads\n", numprocs, omp_get_num_threads());
#else
    printf("Calculating PI with %d processes\n", numprocs);
#endif
    }

// check to see how many threads have started using OpenMP functions
// since we only want this to compile in the parallel version, we use the _OPENMP
// macro name

#pragma omp parallel default(shared) private(i,x) reduction(+:sum)
#ifdef _OPENMP
  printf("process %d thread %d started\n", myrank, omp_get_thread_num());
#endif

// this pragma indicates that the following loop will be done in parallel
//
// each team of threads uses the rank of the parent process to determine
// which intervals to integrate 

#pragma omp for
  for (i = myrank+1; i <= n; i+= numprocs) {
    x = h * ( i - 0.5 ) ;    //calculate at center of interval
    sum += 4.0 / ( 1.0 + pow(x,2) ) ;
    }

// now the threads go idle, and the parent process resumes with the reduced value of sum 

  mypi = h * sum;

// each of the MPI processes calls this MPI routine to reduce mypi 
//
// after this call, rank 0 will contain the reduced value in sumpi

  MPI_Reduce(&mypi,&sumpi,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD) ;

  if (myrank == 0) {
    printf("%f\n",sumpi) ;
    }

// this is called to stop MPI communications

  MPI_Finalize();
  return 0;
}
