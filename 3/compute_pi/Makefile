all:
	gcc-12 -fopenmp -c -O3 compute_pi.c
	gcc-12 -fopenmp main.c compute_pi.o -o a.out
	rm *.o
	./a.out