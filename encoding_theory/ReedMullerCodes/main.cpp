#include <stdlib.h>
#include <stdio.h>

#include <memory.h>
#include <time.h>

#include "rmcoder.h"
#include "bits.h"
#include "noise.h"

#define FREE_IF_ALLOCATED(mem) if (mem != NULL) {free(mem);mem=NULL;}

#define DEBUG 0

void rm_test_0(ull r, ull m, int k, int max_iter, FILE *fout);

int main(int argc, char** argv) {
	int i = 0, j = 0;
	ull gen_size = 0, code_size = 0, data_size = 0;
	ull blocks = 0, wordsize = 0, dim = 0;
	unsigned char *G = NULL, *D = NULL, *DD = NULL, *C = NULL, *NC = NULL;
	unsigned short *masks = NULL;
	unsigned char *code_tmp = NULL, *data_block_tmp = NULL;
	double accumulator = 0;
	char buffer[100];

	FILE *fout = NULL;
	
	srand (time(NULL));

	ull r = 0, m = 16, k = 50, max_iter = 700;
	data_size = 4096;

	char **ptr = NULL;
	
	if (argc == 3) {
		r = strtol(argv[1], ptr, 10);
		m = strtol(argv[2], ptr, 10);
	}

	//rm_test_0(r, m, k, max_iter, stdout);

	for (i = r; i < m; ++i) {
		sprintf(buffer, "output_%2d", i);
		fout = fopen(buffer, "w");
		fprintf(fout, "%d\n", i);
		for (j = 0; j <= i; ++j) {
			printf("%d %d\n", i, j);
			fprintf(fout, "%d %d\n", i, j);
			rm_test_0(j, i, k, max_iter, fout);
			fprintf(fout, "\n");
		}
		fclose(fout);
	}

	//rm_test_0(r, m, k, max_iter);
	// //for (m = 12; m < 16; ++m) 
	// {
	// 	gen_size = generation_mem_size(m, m);
	// 	G = (unsigned char *) malloc(gen_size);
	// 	generate_matrix(m, m, G);

	// 	masks = (unsigned short *) malloc((1 << m) * sizeof(unsigned short));
	// 	generate_masks(m, masks);

	// 	// for (r = 0; r <= m; ++r) 
	// 	{
				
	// 		D =	(unsigned char *) malloc(data_size);
	// 		DD = (unsigned char *) malloc(data_size);

	// 		for (int i = 0; i < data_size; ++i)
	// 			D[i] = (unsigned char) rand();
		     			

	// 		code_size = encode_mem_size(data_size, r, m);
	// 		C = (unsigned char *) malloc(code_size);
	// 		NC = (unsigned char *) malloc(code_size);


	// 		encode(r, m, G, data_size, D, C, masks);
			
	// 		blocks = rm_blocks(data_size, r, m);
	// 		wordsize = rm_wordsize(r, m);
	// 		dim = rm_dim(r, m);


	// 		code_tmp = (unsigned char *) malloc(wordsize / 8 + (wordsize % 8 == 0 ? 0 : 1));
	// 		data_block_tmp = (unsigned char *) malloc(dim / 8 + (dim % 8 == 0 ? 0 : 1));

	// 		decode(r, m, G, data_size, DD, C, masks, code_tmp, data_block_tmp);

	// 		// print(1, 8 * data_size, D);//printf("\n");
	// 		// print(1, 8 * data_size, DD);printf("\n\n");


	// 		//printf("r = %d, m = %d\n", r, m);
	// 		//printf("\nAccuracy %1.9lf:\n\n", accuracy(data_size, D, DD));


	// 		for (ull i = 0; i <= k; i++) {
	// 			accumulator = 0;
	// 			printf("%1.9lf ", (double) i / (double) k);

	// 			for (ull j = 0; j < max_iter; ++j) {
	// 				memcpy(NC, C, code_size);
	// 				add_noise(NC, code_size, (double) i / (double) k);	
	// 				decode(r, m, G, data_size, DD, NC, masks, code_tmp, data_block_tmp);
	// 				accumulator += accuracy(data_size, dim, D, DD);
	// 			}

	// 			printf("%1.9lf \n", accumulator / max_iter);
	// 		}

	// 		FREE_IF_ALLOCATED(C);
	// 		FREE_IF_ALLOCATED(D);
	// 		FREE_IF_ALLOCATED(DD);
			
	// 		FREE_IF_ALLOCATED(NC);

	// 		FREE_IF_ALLOCATED(code_tmp);
	// 		FREE_IF_ALLOCATED(data_block_tmp);
	// 	}

	// 	FREE_IF_ALLOCATED(masks);
	// 	FREE_IF_ALLOCATED(G);
	// }
	return 0;
}

void rm_test_0(ull r, ull m, int k, int max_iter, FILE *fout) {
	ull wordsize = 0, dim = 0, gen_size = 0;
	unsigned char *G = NULL, *D = NULL, *DD = NULL, *C = NULL, *NC = NULL;
	unsigned short *masks = NULL, width = 0;
	ull accumulator = 0, n = 0;

	ull data_block_size_bytes = 0;
	ull code_block_size_bytes = 0;

	dim = rm_dim(r, m);
	wordsize = rm_wordsize(r, m);
	width = rm_codeword_size_in_bytes(r, m);

	gen_size = gen_matrix_mem_size(m, m);
	G = (unsigned char *) malloc(gen_size);
	generate_matrix(m, m, G);
	
	//print(dim, 8 * width, G);printf("\n");

	masks = (unsigned short *) malloc(wordsize * sizeof(unsigned short));
	generate_masks(m, masks);

	data_block_size_bytes = rm_dataword_size_in_bytes(r, m);
	code_block_size_bytes = rm_codeword_size_in_bytes(r, m);

	D =	(unsigned char *) malloc(data_block_size_bytes);
	DD = (unsigned char *) malloc(data_block_size_bytes);

	for (int i = 0; i < data_block_size_bytes; ++i)
		D[i] = (unsigned char) rand();

	C = (unsigned char *) malloc(code_block_size_bytes);
	NC = (unsigned char *) malloc(code_block_size_bytes);

	encode_block(r, m, dim, wordsize, C, D, masks, G);
	//printf("Data:\n");print(1, dim, D);printf("\n");
	
	memcpy(NC, C, width);
	// print(1, wordsize, C);printf("\n");
	// decode_block(r, m, dim, wordsize, NC, DD, masks, G, width);
	// print(1, dim, DD);
	srand (time(NULL));
	
	n = 2 * k / 3 + 1;
	fprintf(fout, "%lld\n", n);
	for (ull i = 0; i < n; i++) {
		accumulator = 0;
		fprintf(fout, "%1.9lf ", (double) i / (double) k);

		for (ull j = 0; j < max_iter; ++j) {
			memcpy(NC, C, code_block_size_bytes);
	
			add_noise(NC, code_block_size_bytes, (double) i / (double) k);				
			decode_block(r, m, dim, wordsize, NC, DD, masks, G, width);
			accumulator += match(data_block_size_bytes, 0, dim, D, DD);
		}

		fprintf(fout, "%1.9lf\n", (double) accumulator / (double) max_iter);
	}
	FREE_IF_ALLOCATED(C);
	FREE_IF_ALLOCATED(NC);
	FREE_IF_ALLOCATED(D);
	FREE_IF_ALLOCATED(DD);
	FREE_IF_ALLOCATED(masks);
	FREE_IF_ALLOCATED(G);
}