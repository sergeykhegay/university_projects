#ifndef _BITS_H_
#define _BITS_H_

#include <stdlib.h>
#include <stdio.h>

#define ull unsigned long long int

// bits matrix
void set(ull n, ull m, unsigned char *A, ull i, ull j, ull val, ull width);
ull get(ull n, ull m, unsigned char *A, ull i, ull j, ull width);
void print(ull n, ull m, unsigned char *A);

// bytes
void set_bit(unsigned char *byte, ull idx, ull val);
ull get_bit(unsigned char *byte, ull idx);
void reverse_bit(unsigned char *byte, ull idx);

int weight(unsigned char *bytes, size_t n);

#endif /* _BITS_H_ */