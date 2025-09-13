CC = mpic++
CFLAGS = -O3 --std=c++17
AFLAGS = -larmadillo
all: simulator.exe
simulator.exe : system.o NSL_SIMULATOR.o chromosome.o random.o
	$(CC) system.o chromosome.o random.o NSL_SIMULATOR.o -o simulator.exe $(AFLAGS)
NSL_SIMULATOR.o : NSL_SIMULATOR.cpp system.h chromosome.h
	$(CC) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o $(CFLAGS)
system.o : system.cpp system.h random.h chromosome.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)
chromosome.o : chromosome.cpp chromosome.h
	$(CC) -c chromosome.cpp -o chromosome.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o simulator.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
