CC		:= g++
ARCH	:= x86_64
CFLAGS	:= -g0 -O2 -Wall -fopenmp

.PHONY: all clean

all:
	@echo "Building main..."
	@$(CC) $(CFLAGS) -o main main.cpp
	@echo "done!"

clean:
	@rm -f main *.o
