#include "bits.h"
#include <time.h>
#include <math.h>
/*
    Returns random number from [0, 1]
*/
double rand_0_1() {
    return fabs((double)rand() / (double)RAND_MAX);
}

/*
	m - array of bytes
	n - number of bytes
	p - probability of changing one bit.
*/
void add_noise(void *m, ull n, double p) {
	ull i = 0, j = 0;
	unsigned char *bytes = (unsigned char *) m;
	//srand (time(NULL));
	
	// for (j = 0; j < 3; ++j) {
	// 	reverse_bit(bytes + i, j);
	// }

	for (i = 0; i < n; ++i) {
		for (j = 0; j < 8; ++j) {
			if (rand_0_1() < p)
				reverse_bit(bytes + i, j);
		}
	}
}