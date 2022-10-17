all: static_lib clean

static_lib: static_lib.a main.c
	gcc main.c static_lib.a -o result.exe
	
static_lib.a: matmul.o helper.o
	ar rc static_lib.a matmul.o helper.o

matmul.o: matmul.c 
	gcc -c matmul.c -O3

helper.o: helper.c
	gcc -c helper.c -g
	
clean:
	rm -rf *.o hello
	
run:
	./result.exe 