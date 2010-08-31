#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[])
{
	long num_steps;
	int i;
	double x, pi, step, sum = 0.0;

	if (argc != 2) {
		fprintf (stderr, "usage: pi num_steps\n");
		exit(-1);
	}
	num_steps = atol (argv[1]);

	step = 1.0/(double) num_steps;
	for (i=1; i<= num_steps; i++) {
		x = (i - 0.5)*step;
		sum = sum + 4.0 /(1.0 + x*x);
	}
	pi = step * sum;
	printf ("Pi is about %.16f\n", pi);
	return 0;
} 
