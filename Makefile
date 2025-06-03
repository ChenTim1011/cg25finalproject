CC = g++
CFLAGS = -std=c++11 -Wall
LIBS = `pkg-config --cflags --libs opencv4`

all: noise_generator.out

noise_generator.out: main.o gabor.o perlin.o
	$(CC) $(CFLAGS) -o noise_generator.out main.o gabor.o perlin.o $(LIBS)

main.o: main.cpp noise.h prng.h
	$(CC) $(CFLAGS) -c main.cpp $(LIBS)

gabor.o: gabor.cpp noise.h prng.h
	$(CC) $(CFLAGS) -c gabor.cpp

perlin.o: perlin.cpp noise.h prng.h
	$(CC) $(CFLAGS) -c perlin.cpp

clean:
	rm -f *.o noise_generator