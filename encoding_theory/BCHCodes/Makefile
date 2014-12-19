all: coder clean

coder: main.o bits.o noise.o
	g++ -o ./build/coder ./build/main.o ./build/bits.o ./build/noise.o

main.o: main.cpp bits.h noise.h
	g++ -c main.cpp -o ./build/main.o

bits.o: bits.cpp bits.h
	g++ -c bits.cpp -o ./build/bits.o

noise.o: noise.cpp noise.h bits.h
	g++ -c noise.cpp -o ./build/noise.o

clean:
	rm -rf ./build/*.o
