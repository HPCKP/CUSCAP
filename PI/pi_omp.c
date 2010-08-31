// compute pi by calculating 4*arctan(1)
// arctan(1) = Integral of (1+x^2)^-1 from 0 to 1
//
// the number of threads can be set at runtime using the OMP_NUM_THREADS 
// environment variable

#include <stdio.h>
#include <math.h>

// have to include omp.h to use OpenMP functions

#include <omp.h>

int main(int argc, char *argv[])
{
  double pi_int = 0.0;  // where we store value of pi that is accumulated
  int n = 100000000 ;  //  number of intervals between 0 and 1
  double h = 1.0 / n ;  // integration interval width 
  double x;
  int i;

// this pragma indicates that the following omp for loop will be 
// executed such that:
//
// -all threads share the same memory location for every variable by default
// -each thread has a "private" copy of the iteration counter `i` and the 
//  intermediate value `x`.   
// -each thread has a private copy of `pi_int`, that is reduced to one value
//  by summation and stored in the pi_int address of the parent process at 
//  the end of the for loop
//
// these distintinctions are commonly called "variable scoping"

#pragma omp parallel default(shared) private(i,x) reduction(+:pi_int)
  {
// check to see how many threads have started using OpenMP functions
// since we only want this to compile in the parallel version, we use the _OPENMP
// macro name

#ifdef _OPENMP  
  printf("thread %d of %d started\n", omp_get_thread_num(), omp_get_num_threads());
#endif

// this pragma indicates that the following loop will be done in parallel

#pragma omp for
  for (i = 1; i <= n; ++i) {
      x = h * ( i - 0.5 ) ;    //calculate at center of interval
      pi_int += 4.0 / ( 1.0 + pow(x,2) ) ; 
    }
  
  }
// now the threads go idle, and the original process resumes with the reduced value of pi_int

  printf("PI = %f\n", h*pi_int);
  return 0;
}
