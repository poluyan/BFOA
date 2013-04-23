TARGET = bfoa
LIBS = -lm
CC = gcc
CFLAGS = -std=c90 -Wall -Wextra

all:
	$(CC) $(CFLAGS) -c main.c
	$(CC) main.o -o $(TARGET) $(LIBS)

clean:
	rm -f *.o
	rm -f $(TARGET)

.PHONY: all clean