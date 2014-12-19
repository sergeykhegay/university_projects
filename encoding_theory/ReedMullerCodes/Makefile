all: coder clean

coder: main.o bits.o rmcoder.o noise.o
	g++ -o ./build/coder ./build/main.o ./build/bits.o ./build/rmcoder.o ./build/noise.o

main.o: main.cpp bits.h rmcoder.h noise.h
	g++ -c main.cpp -o ./build/main.o

bits.o: bits.cpp bits.h
	g++ -c bits.cpp -o ./build/bits.o

rmcoder.o: rmcoder.cpp rmcoder.h bits.h
	g++ -c rmcoder.cpp -o ./build/rmcoder.o

noise.o: noise.cpp noise.h bits.h
	g++ -c noise.cpp -o ./build/noise.o

clean:
	rm -rf ./build/*.o