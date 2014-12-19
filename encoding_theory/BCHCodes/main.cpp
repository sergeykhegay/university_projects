#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

#include <vector>
#include <map>
#include <bitset>
#include <iostream>
#include <set>

#include "noise.h"

using namespace std;

#define lli long long int // __int128_t //
#define ull unsigned long long int

// ELEM begin
lli polynomials[] = {
	0b00000000000000000011, // 1
	0b00000000000000000111, // 2 
	0b00000000000000001011, // 3
	0b00000000000000010011, // 4
	0b00000000000000100101, // 5
	0b00000000000001000011, // 6
	0b00000000000010000011, // 7
	0b00000000000101100001, // 8
	0b00000000001000010001, // 9
	0b00000000010000001001, // 10
	0b00000000100000000101, // 11
	0b00000001000010011001, // 12
	0b00000010000000011011, // 13
	0b00000101100000000011, // 14
	0b00001000000000000011, // 15
	0b00010000000000101101, // 16
};

class Elem {
	lli _pow;
	
public:
	static ull _num;
	static ull _q;
	static ull _m;
	static lli *_pow2pol;
	static lli *_pol2pow;


	Elem(lli pow) {
		if (pow == -1)
			_pow = -1;
		else {
			_pow = pow % _q;
			if (_pow < 0)
				_pow += _q;
		}
	}

	Elem(): _pow(1) {
		;
	}
	
	static Elem primitive() {
		return Elem(1);
	}

	static Elem zero() {
		return Elem(-1);
	}

	bool is_zero() const {
		return (_pow == -1);
	}
	
	lli pow() {
		return _pow;
	}

	friend const Elem operator+(const Elem& left, const Elem& right);
	
	friend const Elem operator-(const Elem& left, const Elem& right);

	friend const Elem operator*(const Elem& left, const Elem& right);

	friend const Elem operator/(const Elem& left, const Elem& right);

	friend bool operator<(const Elem& left, const Elem& right);

   	friend bool operator==(const Elem& left, const Elem& right);
	

	static Elem rev(const Elem& e) {
		lli pow = (-e._pow) % _q;
		if (pow < 0)
			pow += _q;
		return Elem(pow);
	}

	static void print_bits(lli n) {
		bitset<18> temp(n);
    	cout << temp << endl;
	}

	static void gen_pow2pol() {
		lli mask = 1 << _m;
		lli tmp = 1, polinom = polynomials[_m - 1];
		//print_bits(tmp);
		_pow2pol[0] = tmp;

		for (ull i = 1; i < _q; ++i) {
			tmp = tmp << 1;
			if ((tmp & mask) != 0)
				tmp ^= polinom;
			_pow2pol[i] = tmp;
		}

		for (int i = 0; i < _q; ++i)
			print_bits(_pow2pol[i]);

		cout << endl;
	}

	static void gen_pol2pow() {
		_pol2pow[0] = -1;
		for (ull i = 0; i < _q; ++i) {
			_pol2pow[_pow2pol[i]] = i;
		}
		
		for (int i = 1; i < _num; ++i)
			print_bits(_pol2pow[i]);
	}

	/*
		whole class represents a field 
		of 2 ^ m elements
	*/
   	static void initialize(int m) {
   		_m = m;
		_num = 1 << m;
		_q = _num - 1;

		_pow2pol = (lli *) malloc(_num * sizeof(lli));
		_pol2pow = (lli *) malloc(_num * sizeof(lli));

		gen_pow2pol();
		gen_pol2pow();
   	}

   	static void destruct() {
		free(_pow2pol);
		free(_pol2pow);
   	}

   	void print() {
   		printf("%lld\n", _pow);
   	}
};

ull Elem::_m = 0;
ull Elem::_q = 0;
ull Elem::_num = 0;
lli* Elem::_pow2pol = NULL;
lli* Elem::_pol2pow = NULL;

const Elem operator+(const Elem& left, const Elem& right) {
	if (left.is_zero())
		return right;
	if (right.is_zero())
		return left;

	lli polin = Elem::_pow2pol[left._pow] ^ Elem::_pow2pol[right._pow];
	return Elem(Elem::_pol2pow[polin]);
}
	
const Elem operator-(const Elem& left, const Elem& right) {
	return (left + right);
}

const Elem operator*(const Elem& left, const Elem& right) {
	if (left.is_zero() || right.is_zero())
		return Elem::zero();

	return (left._pow + right._pow) % Elem::_q;
}

const Elem operator/(const Elem& left, const Elem& right) {
	if (right.is_zero())
		throw "Division by zero element";
	return left * Elem::rev(right);
}

bool operator<(const Elem& left, const Elem& right) {
      return left._pow < right._pow;
}

bool operator==(const Elem& left, const Elem& right) {
	return left._pow == right._pow;
}

// ELEM END


class Coder {
	lli _num;		// total number of elements in the field
	lli _n;			// number of elements in multiplicative group
	lli _t;			// number of errors possible to correct
	lli _k;			// number of data bits
	lli _code_pol;	// code generative polynomial  00010011 -> x ^ 4 + x + 1
					// maximum polynomial degree = 63
	lli _deg;		// code polynomial degree
	lli _mask;
	lli _shift;     // _n - _k

protected:
	set<Elem> compute_cyclotomic_class(Elem elem) {
		set<Elem> s;
		Elem initial = elem;

		s.insert(elem);
		while (true) {
			elem = elem * elem;
			if (elem == initial)
				break;

			s.insert(elem);
		}

		return s;
	}

	set<Elem> compute_roots() {
		set<Elem> roots;
		set<Elem> tmp;

		// TODO: Why for all odds?? :'()
		for (ull i = 1; i <= 2 * _t; i += 2) { 
			tmp = compute_cyclotomic_class(Elem(i));
			roots.insert(tmp.begin(), tmp.end());
		}
		return roots;
	}
	
	void multiply(vector<Elem> &v, vector<Elem> &tmp, Elem elem, ull count) {
		tmp = v;
		printf("Elem :"); elem.print();

		for (ull i = count; i > 0; --i)  // v(x) * x
			v[i] = v[i-1];
		
		v[0] = Elem::zero();

		for (ull i = 0; i < count; ++i)  // v(x) * elem
			tmp[i] = tmp[i] * elem;


		for (ull i = 0; i < count; ++i)  // v'(x) = v(x) * (x + elem)
			v[i] = v[i] + tmp[i];
	}

	lli generate_code_polynomial(lli &pol_degree) {
		lli polynomial = 0, one = 1;
		
		set<Elem> roots = compute_roots();
		set<Elem>::iterator it;
		ull num_of_roots = roots.size();

		printf("number of roots: %lld\n", num_of_roots);
		ull count = 1;

		vector<Elem> v(num_of_roots + 1, Elem::zero());
		vector<Elem> tmp;

		v[0] = Elem(0);	

		for (it = roots.begin(); it != roots.end(); it++) {
			multiply(v, tmp, *it, count);
			++count;
		}
		
		for (ull i = 0; i < count; ++i) {
			if (v[i].pow() == 0) 
				polynomial |= (one << i);
		}
		
		printf("g(x) = ");
		bitset<33> temp(polynomial);
    	cout << temp << endl;
		
		pol_degree = num_of_roots;
		return polynomial;
	}

	lli generate_mask() {
		lli res = 0;
		lli tmp = 1;
		for (int i = 0; i < _k; ++i)
			res |= tmp << i;

		return res;
	}

	lli compute_reminder(lli pol) {
		ull pow = _n - 1;
		lli divider = _code_pol;
		lli mask = 1;

		while (pow >= _deg) {
			if ((pol & (mask << pow)) != 0)
				pol ^= (divider << (pow - _deg));
			pow--;
		}

		return pol;
	}	

public:
	/*
		m - power of 2
		t - number of errors
		d - data length bits
	*/
	Coder(int m, int t): _t(t) {
		Elem::initialize(m);

		_num = 1 << m;
		_n = _num - 1;
		_code_pol = generate_code_polynomial(_deg);
		_mask = generate_mask();
		_k = _n - _deg;
		_shift = _n - _k;
	}

	~Coder() {
		Elem::destruct();
	}

	lli encode(lli data) {
		lli reminder = 0;
		lli data_polynomial = data & _mask; // nullify higher bits from _k's 
											// position and higher.
		lli higher_bits = data_polynomial << _shift;
		
		reminder = compute_reminder(higher_bits);

		return (reminder | higher_bits);
	}
	
	Elem compute_codeword(lli codeword, Elem elem) {
		Elem res = Elem::zero();
		Elem mult = Elem(0);
		lli mask = 1;

		for (ull i = 0; i < _num; i++) {
			if ((codeword & (mask << i)) != 0)
				res = res + mult;
			mult = mult * elem;
		}

		return res;
	}

	bool has_errors(lli codeword) {
		Elem primitive = Elem::primitive();
		Elem tmp = Elem::primitive();
		
		for (ull i = 1; i <= 2 * _t; i++) {
			if (! (compute_codeword(codeword, tmp) == Elem::zero()) )
				return true;
			tmp = tmp * primitive;
		}

		return false;
	}

	ull num_of_data_bits() {
		return _n - _deg;
	}

	void whoiam() {
		printf("\n");
		printf("BCH(%lld, %lld)\n", _n - _deg, _n);
		printf("data bits number: %lld\n", _n - _deg);
		printf("deg(g(x)): %lld\n", _deg);
		printf("codeword length: %lld\n", _n);
		printf("\n");
	}
};

void test0(int m, int t);

int main(int argc, char** argv) {
	srand (time(NULL));
	// Coder coder(3, 1);
	// coder.whoiam();
	// lli code = coder.encode(1);
	
	// code ^= 8; 
	// if (coder.has_errors(code))
	// 	cout << "Has error" << endl;
	// else
	// 	cout << "No errors" << endl;

	// bitset<64> temp(code);
 //    cout << temp << endl;
	
	char **ptr = NULL;
	int m = 0, t = 0;
	if (argc == 3) {
		m = strtol(argv[1], ptr, 10);
		t = strtol(argv[2], ptr, 10);
	}
	if (t <= 0 || m <= 0)
		return 0;
	test0(m, t);

	return 0;
}

double power(double base, int n) {
	double res = 0;
	if (n == 1)
		return base;
	
	res = power(base, n/2);

	if (n % 2 == 1)
		return res * res * base;
	else
		return res * res;
}

void test0(int m, int t) {
	FILE *fout = stdout;
	ull k = 100;
	int max_iter = 4000;
	double accumulator = 0;
	lli C, NC;
	ull n = (1 << m) - 1;
	int code_block_size_bytes = (int) sizeof(lli);

	Coder coder(m, t);
	
	coder.whoiam();

	
	C = coder.encode(0b11111);
	for (ull i = 0; i <= k; ++i) {
		accumulator = 0;
		fprintf(fout, "%1.9lf ", (double) i / (double) k);

		for (ull j = 0; j < max_iter; ++j) {
			memcpy(&NC, &C, code_block_size_bytes);
	
			add_noise(&NC, code_block_size_bytes, (double) i / (double) k);						
			if (! coder.has_errors(NC))
				accumulator++;
		}

		fprintf(fout, "%1.9lf ", (double) accumulator / (double) max_iter);
		fprintf(fout, "%1.9lf\n", power(1 - (double) i / (double) k, n));

	}
}