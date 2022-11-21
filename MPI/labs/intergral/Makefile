all: main.out

main.out: main.c
	gcc -pthread -Wall -O3 $^ -lm -o $*.out

clean:
	rm -rf main.out

.PHONY: clean