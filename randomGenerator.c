#include <stdlib.h> // using `exit`
#include <stdio.h>  // using `puts`
#include <math.h>   // using `fmod`

#define AA    471
#define BB    1586
#define CC    6988
#define DD    9689
#define MM    16383
#define RIMAX 2147483648.0  // 2^31
#define RandomInteger (++nd, ra[nd & MM] = ra[(nd-AA) & MM] ^ ra[(nd-BB) & MM] ^ ra[(nd-CC) & MM] ^ ra[(nd-DD) & MM])

void seed(long seed);
static long ra[MM+1], nd;

void seed(long seed) {
  int i;
	if(seed < 0) {puts("SEED error."); exit(1);}
	ra[0] = (long)fmod(16807.0*(double)seed, 2147483647.0);
	for(i = 1; i <= MM; i++) ra[i] = (long)fmod(16807.0 * (double)ra[i-1], 2147483647.0);
}

long randl(long n) { // random long integer between 0 and n-1; n must be larger than 0
  return(RandomInteger % n);
}

double randd(void) { // random double between 0 and 1
  return((double) RandomInteger / RIMAX);
}

double randf(double m) { // random double between 0 and m
  return(m*((double) RandomInteger / RIMAX));
}

unsigned long randul(unsigned long n) { // random unsigned long integer between 0 and n-1; n must be larger than 0
// NOTE: No extensive testing was done.
	return((unsigned long)(randd()*n));
}

int randomBinomial(int n, double p) { // random variable from a binomial distribution of choice probabilty `p` (between 0 and `n`)
  int i, c = 0;
  for (i = 0; i < n; i++) if(randd() < p) c++;
	return(c);
}

// Source: http://stackoverflow.com/a/311716/712498
void randomSample(int n, int k, int *sample) { // samples k from n without replacement
    int t = 0, m = 0; // t = total input records dealt with; m = number of items selected so far
    double r;
    while (m < k) {
        r = randd(); // uniform (0,1) real
        if ((n-t)*r >= k-m) t++;
        else {
					sample[m] = t;
          t++; m++;
        }
    }
}

int randomChoice(double *v, double max, int n) { // returns an integer `i` from the range (0, n-1) with probability proportional to `v[i]`
	int i = 0;
	double sum = 0, r;
	//for(i = 0; i < n; i++) max +=v[i]; // if `max` is not supplied and has to calculate locally
	r = randd()*max;
	while((sum < r) && (i < n)) {
		sum += v[i];
		i++;	
	}
	return(i-1);
}
