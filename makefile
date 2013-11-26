fft: main.c fft.c fft.h util.c util.h
	gcc -o fft main.c fft.c util.c -lm --std=c99

dft: main-d.c fft.c fft.h util.c util.h
	gcc -o dft main-d.c fft.c util.c -lm --std=c99

test: tests/main.c fft.c fft.h util.c util.h
	gcc -g -o ftest tests/main.c fft.c util.c -lm --std=c99
