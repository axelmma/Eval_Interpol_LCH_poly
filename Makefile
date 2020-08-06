CC=gcc 
LIBS= -lflint -lgmp -lmpfr

exe: main.o LCH.o 
	$(CC) $(LIBS) -o exe LCH.o main.o 

main.o: main.c LCH.h
	$(CC) -c -DDEBUG main.c

LCH.o: LCH.c
	$(CC) -c -DDEBUG LCH.c 

clean: 
	rm -rf *.o

