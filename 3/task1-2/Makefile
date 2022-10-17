CC=gcc-12
all: task_2 

task_2:
	# $(CC) omp_bugreduction.c  -fopenmp -o omp_bugreduction.exe
	# ./omp_bugreduction.exe
	$(CC) omp_bugparfoc.c  -fopenmp -o omp_bugparfor.exe
	# ./omp_bugparfor.exe

demonstrated:
	$(CC) omp_hello.c -fopenmp -o omp_hello.exe 
	$(CC) omp_outmes.cpp -fopenmp -o omp_outmes.exe 
	$(CC) omp_privateshared.c -fopenmp -o omp_privateshared.exe 
	$(CC) omp_parsec.c -fopenmp -o omp_parsec.exe 
	$(CC) omp_sumarray.c -fopenmp -o omp_sumarray.exe 
	rm -rf *.o