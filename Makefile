CC=gcc
CFLAGS=-Wall -O3
LIBS=-lm

SRC=main.c data_utils.c cmd_parser.c ode_solvers.c vects.c
OUT=simulation

all: $(OUT) $(OUT)

$(OUT): $(SRC)
	$(CC) $(CFLAGS) -o $(OUT) $(SRC) $(LIBS)

clean:
	rm -f $(OUT) *.png out/*

.PHONY: all clean
