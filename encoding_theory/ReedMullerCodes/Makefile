all: coder clean

coder: main.o bits.o rmcoder.o noise.o
	g++ -o coder main.o bits.o rmcoder.o noise.o

main.o: main.cpp bits.h rmcoder.h noise.h
	g++ -c main.cpp

bits.o: bits.cpp bits.h
	g++ -c bits.cpp

rmcoder.o: rmcoder.cpp rmcoder.h bits.h
	g++ -c rmcoder.cpp

noise.o: noise.cpp noise.h bits.h
	g++ -c noise.cpp

clean:
	rm -rf *.o