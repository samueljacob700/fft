#ifndef UTIL_H_
#define UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { double real; double imag; } Cnum; // num[0] is "real", num[1] "imaginary"

double magn(Cnum *);
double mag2(Cnum *);

Cnum *read_file_in(char *path);

#endif
