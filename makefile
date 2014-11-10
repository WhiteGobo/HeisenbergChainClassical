CC=gcc
CCFLAGS=-lm  -D_REENTRANT -lpthread

PROGRAM=main

all: program

program: timechain.o maincycleall.o plotter.o reader.o
	$(CC) $(CCFLAGS) -I src/ maincycleall.o timechain.o plotter.o reader.o -o program2

timechain.o: src/timechain/timechain.c src/plotter/plotter.h src/reader/reader.h
	$(CC) $(CCFLAGS) -c src/timechain/timechain.c src/plotter/plotter.h src/reader/reader.h

plotter.o: src/plotter/plotter.c
	$(CC) $(CCFLAGS) -c src/plotter/plotter.c

reader.o: src/reader/reader.c
	$(CC) $(CCFLAGS) -c src/reader/reader.c

main.o: src/main.c src/timechain/timechain.h
	$(CC) $(CCFLAGS) -c src/main.c src/timechain/timechain.h

maincycleall.o: src/maincycleall.c src/timechain/timechain.h
	$(CC) $(CCFLAGS) -c src/maincycleall.c src/timechain/timechain.h

clean:
	rm *.o
