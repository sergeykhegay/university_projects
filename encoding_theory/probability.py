import math
import sys

def newton_coef(n, m):
	return math.factorial(n) / (math.factorial(n - m) * math.factorial(m))

def newton_coef_list(n, m):
	# print n, m
	return [newton_coef(n, i) for i in xrange(m + 1)]

def probability(r, m, p):
	res = 0
	total = 2 ** m
	corrupted = int(2 ** (m - r - 1))

	C = newton_coef_list(total, corrupted)
	P = [p ** i * (1 - p) ** (total - i) for i in xrange(corrupted + 1)]
	for a, b in zip(C, P):
		res += a * b

	return res

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "Wrong usage: python probability.py <r> <m> <prob>"
	else:
		r = int(sys.argv[1])
		m = int(sys.argv[2])
		p = float(sys.argv[3])
		print probability(r, m, p)