#include "bits.h"

void print(ull n, ull m, unsigned char *A) {
	ull i = 0, j = 0;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < m; ++j)
			printf("%lld", get(n, m, A, i, j, m));
		printf("\n");
	}
}

void set(ull n, ull m, unsigned char *A, ull i, ull j, ull val, ull width) {
    ull byte_idx = (i * width + j) / 8;
    ull idx = 7 - (i * width + j) % 8;
    set_bit(A + byte_idx, idx, val);
}

ull get(ull n, ull m, unsigned char *A, ull i, ull j, ull width) {
    ull byte_idx = (i * width + j) / 8;
    ull idx = 7 - (i * width + j) % 8;
    return get_bit(A + byte_idx, idx);
}

void set_bit(unsigned char *byte, ull idx, ull val) {
	unsigned char mask = 1;
	
	mask = mask << idx;
	if (val == 0) {
		mask = ~mask;
		*byte = *byte & mask;
	}
	else
		*byte = *byte | mask;
}

ull get_bit(unsigned char *byte, ull idx) {
	unsigned char mask = 1;

	if (idx > 7 || idx < 0)
		return 0;
	
	if ( (*byte & (mask << idx)) != 0 )
		return 1;

	return 0; 
}

void reverse_bit(unsigned char *byte, ull idx) {
	if (get_bit(byte + idx / 8, idx % 8) == 1)
		set_bit(byte + idx / 8, idx % 8, 0);
	else
		set_bit(byte + idx / 8, idx % 8, 1);
}


int weight(unsigned char *bytes, size_t n) {
	int i = 0, j = 0;
	unsigned char mask = 1;
	int count = 0;

	for (i = 0; i < n; ++i) {
		mask = 1;
		for (j = 0; j < 8; ++j) {
			if ( (*(bytes + i) & mask)  != 0 )
				++count;
			mask = mask << 1;
		}
	}

	return count;
}