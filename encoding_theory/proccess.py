import math
import sys
import probability


def foo(filename):
	with open(filename, "r") as fin:
		fout = open("{0}_out".format(filename), "w")
		m = int(fin.readline())
		
		fout.write("{0}\n".format(m))
		for i in xrange(m + 1):
			r = int(fin.readline().split()[1].strip())
			n = int(fin.readline().strip())
			# print r, m, n
			fout.write("{0} {1}\n{2}\n".format(m, r, n))
			for j in xrange(n):
				arg, val = (float(word) for word in fin.readline().split())
				# print r, m, arg
				p = probability.probability(r, m, arg)
				fout.write("{0:0.9f} {1:0.9f} {2:0.9f}\n".format(arg, val, p))
			fin.readline()
			fout.write("\n")

		fout.close()

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "Wrong usage: python proccess.py <filename>"
	else:
		filename = sys.argv[1]
		foo(filename)