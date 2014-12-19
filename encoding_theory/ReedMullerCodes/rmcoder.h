#ifndef _RMCODER_H_
#define _RMCODER_H_

#include "bits.h"

ull rm_dim(ull r, ull m);
/*
	RM code generative matrix size in bytes.
*/
ull gen_matrix_mem_size(ull r, ull m);
void generate_matrix(ull r, ull m, unsigned char *A);

ull rm_dataword_size_in_bits(ull r, ull m);
ull rm_codeword_size_in_bits(ull r, ull m);
ull rm_dataword_size_in_bytes(ull r, ull m);
ull rm_codeword_size_in_bytes(ull r, ull m);
ull rm_wordsize(ull r, ull m);
ull rm_blocks(ull n, ull r, ull m);
ull encode_mem_size(ull n, ull r, ull m);



void encode(ull r, ull m, unsigned char *G, ull data_n, 
	unsigned char *data, unsigned char *C, unsigned short *masks);
void decode(ull r, ull m, unsigned char *G, ull data_n, unsigned char *data,
	unsigned char *C, unsigned short *masks, 
	unsigned char *code_tmp, unsigned char *data_block_tmp);


void encode_block(ull r, ull m, ull dim, ull wordsize, unsigned char *code_block, 
	unsigned char *data_block, unsigned short *masks, unsigned char *G);

void decode_block(ull r, ull m, ull dim, ull wordsize, unsigned char *code_block, 
	unsigned char *data_block, unsigned short *masks, unsigned char *G, ull width);

int match(ull n, ull idx, ull block_size, unsigned char *A, unsigned char *B);
double accuracy(ull n, ull block_size, unsigned char *A, unsigned char *B);

ull masks_mem_size(int max_num_of_args);
void generate_masks(int n, unsigned short *masks);

#endif /* _RMCODER_H_ */