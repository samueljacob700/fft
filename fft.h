#ifndef FFT_H_
#define FFT_H_
#include "util.h"

void dft(const int, Cnum*, Cnum*);
void sin_wave(int, double, double, Cnum*);
void cos_wave(int, double, double, Cnum*);
void sqr_wave(int, double, double, Cnum*);

Cnum* fft(const int, Cnum* /*, int*/);
Cnum* fft_combine(const int, Cnum*, Cnum* /*, int*/);


#endif
