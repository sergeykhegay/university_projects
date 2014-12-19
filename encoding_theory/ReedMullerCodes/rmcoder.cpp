#include "bits.h"
#include <memory.h>

#define DEBUG 0

ull rm_dim(ull r, ull m) {
	ull res = 1;
	ull tmp = 1;

	if (m < 0)
		res = 0;
	else if (r >= m) 
		res = res << m;
	else if (r < 0)
		res = 0;
	else {
		for (ull i = 1; i <= r; ++i) {
			tmp = (m - i + 1) * tmp / i;
			res += tmp;
		}
	}

	return res;
}

ull rm_codeword_size_in_bits(ull r, ull m) {
	return 1 << m;
}

ull rm_codeword_size_in_bytes(ull r, ull m) {
	ull tmp = rm_codeword_size_in_bits(r, m);
	return tmp / 8 + (tmp % 8 == 0 ? 0 : 1);
}

ull rm_dataword_size_in_bits(ull r, ull m) {
	return rm_dim(r, m);
}

ull rm_dataword_size_in_bytes(ull r, ull m) {
	ull tmp = rm_dataword_size_in_bits(r, m);
	return tmp / 8 + (tmp % 8 == 0 ? 0 : 1);
}

ull rm_wordsize(ull r, ull m) {
	return 1 << m;
}

ull rm_data_blocks_num(ull data_n, ull r, ull m) {
	ull info_bits = 8 * data_n;
	ull dim = rm_dim(r, m);

	return info_bits / dim + ((info_bits % dim) == 0 ? 0 : 1);
}

/*
	RM code generative matrix size in bytes.
*/
ull gen_matrix_mem_size(ull r, ull m) {
	ull bytes = rm_dim(r, m) * rm_codeword_size_in_bytes(r, m);
	return bytes;
}

// from -> to
void copy(unsigned char *A, ull n, ull m, 
	ull height, ull width, ull r_0, ull c_0, ull r_1, ull c_1, ull width_b) {
	ull i = 0, j = 0, val = 0;

	for (i = 0; i < height; ++i)
		for (j = 0; j < width; ++j) {
			val = get(n, m, A, r_0 + i, c_0 + j, width_b);
			set(n, m, A, r_1 + i, c_1 + j, val, width_b);
		}
}

void G(ull r, ull m, unsigned char *A, ull N, ull M, ull idx, ull jdx, ull width_b) {
	ull width = 0, height = 0;
	ull j = 0;

	height = rm_dim(r, m);
	width = 1 << m;

	if (m < 0 || height == 0) return;

	
	if (height == 1) {
		for (j = 0; j < width; ++j)
			set(N, M, A, idx, jdx + j, 1, width_b);
	}
	else if (width == 1)
		set(N, M, A, idx, jdx, 1, width_b);
	else {
		height = rm_dim(r, m - 1);
		width = 1 << ((m - 1 > 0) ? m - 1 : 0);

		G(r, m - 1, A, N, M, idx, jdx, width_b);
		copy(A, N, M, height, width, idx, jdx, idx, jdx + width, width_b);
		G(r - 1, m - 1, A, N, M, idx + height, jdx + width, width_b);
	}

}

void generate_matrix(ull r, ull m, unsigned char *A) {
	ull N, M, i = 0, s = 0, width = 0;

	memset(A, 0, gen_matrix_mem_size(r, m));
	N = rm_dim(r, m);
	M = rm_wordsize(r, m);
	width = 8 * rm_codeword_size_in_bytes(r, m);
	G(r, m, A, N, M, 0, 0, width);
}

/*
	memory needed for encoded info
*/
ull encode_mem_size(ull n, ull r, ull m) {
	ull blocks, bits, wordsize;
	ull bytes;

	wordsize = rm_wordsize(r, m);
	blocks = rm_data_blocks_num(n, r, m);
	bits = blocks * wordsize;

	bytes = bits / 8 + ((bits % 8) == 0 ? 0 : 1);

	return bytes;
}

void copy_vec(ull len, unsigned char *from, ull from_i, unsigned char *to, ull to_i) {
	ull i = 0;
	int tmp = 0;
	for (i = 0; i < len; ++i) {
		tmp = get_bit(from + (from_i + i) / 8, 7 - (from_i + i) % 8);
		set_bit(to + (to_i + i) / 8, 7 - (to_i + i) % 8, tmp);
	}
} 

void multiply(unsigned char *V, ull V_m, ull V_i, unsigned char *G, ull n, ull m, 
	unsigned char *C, ull C_m, ull C_i, unsigned short *masks) {
	ull i = 0, j = 0;
	ull tmp = 0, a = 0, g = 0;

	for (j = 0; j < m; ++j) {
		tmp = 0;
		for (i = 0; i < n; ++i) {
			a = get(1, V_m, V, 0, V_i + i, V_m);
			g = get(n, m, G, masks[i], j, m);
			tmp = tmp + a * g;
		}
		set(1, C_m, C, 0, C_i + j, tmp % 2, C_m);
	}
}

void multiply_b(unsigned char *V, ull V_m, ull V_i, unsigned char *G, ull n, ull m, 
	unsigned char *code_tmp, ull width, unsigned short *masks) {
	ull i = 0, j = 0;
	ull tmp = 0, data_bit = 0, g = 0, row = 0;
	// ull width_in_bytes = width / 8;

	memset(code_tmp, 0, width);
	
	for (i = 0; i < n; ++i) {
		data_bit = get(1, V_m, V, 0, V_i + i, V_m);
		if (data_bit == 1) {
			row = masks[i];
			for (j = 0; j < width; ++j) {
				code_tmp[j] ^= *(G + row * width + j);
			}
		}	
		//printf("dd   "); print(1, 8 * width, code_tmp);
	}
}

void encode_block(ull r, ull m, ull dim, ull wordsize, unsigned char *code_block, 
	unsigned char *data_block, unsigned short *masks, unsigned char *G) {
	ull width = rm_codeword_size_in_bytes(r, m);
	multiply_b(data_block, dim, 0, G, dim, wordsize, code_block, width, masks);
}

void encode(ull r, ull m, unsigned char *G, ull data_n, unsigned char *data, 
	unsigned char *C, unsigned short *masks,
	unsigned char *code_tmp, unsigned char *data_block_tmp) {

	ull i = 0, j = 0;
	ull dim = rm_dim(r, m) , wordsize = rm_wordsize(r, m);
	ull blocks = rm_data_blocks_num(data_n, r, m);
	ull data_bits = data_n * 8;
	ull width = rm_codeword_size_in_bytes(r, m);

	for (i = 0; i < blocks; ++i) {

		if ((i + 1) * dim <= data_bits) {
			multiply_b(data, data_bits, i * dim, G, dim, wordsize, 
				code_tmp, width, masks);
		}
		else {
			multiply_b(data, data_bits, i * dim, G, data_bits % dim + 1, wordsize, 
				code_tmp, width, masks);
		}
		copy_vec(wordsize, code_tmp, 0, C, i * wordsize);
	}
}

void short_set(unsigned short *num, int idx, int val) {
	unsigned short mask = 1;
	
	mask = mask << idx;
	if (val == 0)
		*num &= ~mask;
	else
		*num |= mask;
}

int short_get(unsigned short num, int idx) {
	unsigned short mask = 1;
	return (num & (mask << idx)) != 0 ? 1 : 0; 
}





unsigned short build_index(int m, int const_n, unsigned short const_val, 
	int var_n, unsigned short var_val, unsigned short mask) {
	int i = 0, tmp = 0;
	unsigned short res = 0;

	for (i = m - 1; i >= 0; --i) {
		if (short_get(mask, i) == 1) {
			--var_n;
			tmp = short_get(var_val, var_n);
		}
		else {
			--const_n;
			tmp = short_get(const_val, const_n);
		}
		short_set(&res, i, tmp);	
	}
	return res;
}

int short_weight(unsigned short val) {
	int i = 0;
	unsigned short mask = 1;
	int count = 0;

	for (i = 0; i < 16; ++i) {
		if ( (val & mask) != 0 )
			++count;

		mask = mask << 1;
	}

	return count;
}

int get_coef(int i, int r, int m, unsigned char *codeblock, unsigned short mask) {
	int cout_ones = 0, count_zeros = 0, sum = 0;
	int var_n, const_n, val;
	unsigned short tmp = 0, var_val = 0, const_val = 0;

	var_n = short_weight(mask);
	const_n = m - var_n;
	
	if (DEBUG) {
		printf("\nmask:        ");
		print(1, 16, (unsigned char *) &mask);
	}

	for (const_val = 0; const_val < (1 << const_n); ++const_val) {
		sum = 0;

		if (DEBUG) {
			printf("const_val\n");
			print(1, 16, (unsigned char *) &const_val);
		}

		for (var_val = 0; var_val < (1 << var_n); ++var_val) {
			tmp = build_index(m, const_n, const_val, var_n, var_val, mask);
			
			if (DEBUG) {
				printf("      ");
				print(1, 16, (unsigned char *) &tmp);
			}

			val = get_bit(codeblock + tmp / 8, 7 - tmp % 8);
			if (DEBUG) {
				printf("		index : %d  val = %d\n", tmp, val);
			}
			sum += val;
		}

		sum %= 2;
		if (sum == 1)
			++cout_ones;
		else
			++count_zeros;
	}

	if (DEBUG) {
		printf("			result %d\n", cout_ones >= count_zeros);
	}
	return cout_ones >= count_zeros;
}

void substract_monomial(int row, ull dim, ull wordsize, unsigned char *G, 
	unsigned char *codeblock, ull width) {
	int i = 0, tmp = 0;

	if (DEBUG) {
		printf("\n SUBTRACT \n codeblock\n");
		print(1, wordsize, codeblock);
		// printf("val\n");
		// print(wordsize, wordsize, G);
		// printf("\n");
	}
	
	for (i = 0; i < width; ++i) {
		*(codeblock + i) ^= *(G + row * width + i);
	}
	//print(1, wordsize, codeblock);printf("\n");
}

void decode_block(ull r, ull m, ull dim, ull wordsize, unsigned char *code_block, 
	unsigned char *data_block, unsigned short *masks, unsigned char *G, ull width) {
	int i = 0, coef = 0;

	for (i = dim - 1; i >= 0; --i) {
		coef =  get_coef(i, r, m, code_block, masks[i]);
		if (coef == 1) substract_monomial(masks[i], dim, wordsize, G, code_block, width);
		set_bit(data_block + i / 8, 7 - i % 8, coef);

		if (DEBUG) {
			printf("data block:\n");
			print(1, dim, data_block);
		}
	}
}

void decode(ull r, ull m, unsigned char *G, ull data_n, unsigned char *data,
	unsigned char *C, unsigned short *masks, 
	unsigned char *code_tmp, unsigned char *data_block_tmp) {

	int i = 0, j = 0;
	ull blocks = 0, wordsize = 0, dim = 0, width = 0;

	blocks = rm_data_blocks_num(data_n, r, m);
	wordsize = rm_wordsize(r, m);
	dim = rm_dim(r, m);
	width = rm_codeword_size_in_bytes(r, m);

	for (i = 0; i < blocks; ++i) {
		//printf("i = %d, blocks = %d\n", i, blocks);
		copy_vec(wordsize, C, i * wordsize, code_tmp, 0);
		//printf("    copied\n");
		decode_block(r, m, dim, wordsize, code_tmp, data_block_tmp, masks, G, width);
		//printf("    decoded\n");

		copy_vec( (dim < data_n * 8 - i * dim ? dim:  data_n * 8 - i * dim), data_block_tmp, 0, data, i * dim); // TODO FIX FIRST ARG
		//printf("    copied\n");
	}

}



int match(ull n, ull idx, ull block_size, unsigned char *A, unsigned char *B) {
	int val_A = 0, val_B  = 0;
	int i = 0;
	for (i = 0; i < block_size; ++i) {
		val_A = get(1, n, A, 0, idx + i, n);
		val_B = get(1, n, B, 0, idx + i, n);
		if (val_A != val_B)
			return 0;
	}
	return 1;
}
/*
	n - bytes.
*/
double accuracy(ull n, ull block_size, unsigned char *A, unsigned char *B) {
	ull count = 0, i = 0;
	int N = 8 * n / block_size;
	
	for (i = 0; i < N; ++i) {
		if (match(n, i * block_size, block_size, A, B) != 1)
			++count;
	}
	return 1 - (double) count / (double) (N);
}


int compare_ushorts(const void *a, const void *b) {
  	int val_a = weight((unsigned char*) a, 2), val_b = weight((unsigned char*) b, 2);

  	if (val_a > val_b)
  		return 1;
  	else if (val_a < val_b)
    	return -1;
  	else {
  		if (a == b)
    		return 0;
    	else if (a < b)
    		return -1;
    	else
    		return 1;
    }
}

ull masks_mem_size(int max_num_of_args) {
	return 1 << max_num_of_args;
}

void generate_masks(int n, unsigned short *masks) {
	unsigned short i = 0;

	for (i = 0; i < (1 << n); ++i)
		*(masks + i) = i;
	
	qsort(masks, 1 << n, sizeof(unsigned short), compare_ushorts);
}	