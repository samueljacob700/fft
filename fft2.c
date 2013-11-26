#include "fft.h"

// assumes N is array length of all args
// 
void dft(const int N, Cnum *in, Cnum *out) {
    for (int k = 0; k < N; k++) {
        double r_sum = 0;
        double i_sum = 0;
        
        for (int n = 0; n < N; n++) {
            Cnum *innum = &in[n];
            r_sum +=  innum->real * cos( 2 * M_PI * n * k / N ) +
                innum->imag * sin( 2 * M_PI * n * k / N );
            i_sum += -innum->real * sin( 2 * M_PI * n * k / N ) +
                innum->imag * cos( 2 * M_PI * n * k / N );
        }
        out[k] = (Cnum){r_sum, i_sum};
    }
}

Cnum *fft_combine(const int N, Cnum *E, Cnum *O) {
    Cnum *X = (Cnum*)malloc(sizeof(Cnum)*N);

    for (int k = 0; k < N/2; k++) {
        Cnum t = (Cnum) { ( O[k].real * cos(2*M_PI*k/N)
                           +O[k].imag * sin( 2*M_PI*k/N)),
                            
                          (-O[k].real * sin(2*M_PI*k/N)
                           +O[k].imag * cos(2*M_PI*k/N)) };
        X[k    ] = (Cnum) { E[k].real + t.real,
                            E[k].imag + t.imag };
        X[k+N/2] = (Cnum) { E[k].real - t.real,
                            E[k].imag - t.imag };
    }
    return X;
}

// adapted from en.literateprograms.org/Cooley-Tukey_FFT_algorithm_(C)
// don't need stride, because new arrays are created -- not memory optimized
Cnum *fft(const int N, Cnum *in/*, int step*/) { 
    Cnum *e_in, *o_in, *E, *O;

    if (N % 2)
    {
        Cnum *X = (Cnum*)malloc(sizeof(Cnum)*N);
        dft(N, in, X);
        /*X[0] = in[0];*/
        return X;
    }

    e_in = (Cnum*)malloc(sizeof(Cnum)*N/2);
    o_in = (Cnum*)malloc(sizeof(Cnum)*N/2);

    for (int k = 0; k < N/2; k++)
    {
        e_in[k] = in[2*k];
        o_in[k] = in[2*k + 1];
    }

    E = fft(N/2, e_in);
    O = fft(N/2, o_in);

    Cnum *X = fft_combine(N, E, O);
    // don't free input, or memory will be double-freed after recursion
    free(E);
    free(O);

    return X;
}

void sin_wave(int samples, double cycles, double phase, Cnum *out) {
        double period = samples / cycles;
        for (int i = 0; i < samples; i++) {
            double s = sin(phase + 2*M_PI * i / period);
            out[i] = (Cnum){s, 0};
        }
}
void cos_wave(int samples, double cycles, double phase, Cnum *out) {
        double period = samples / cycles;
        for (int i = 0; i < samples; i++) {
            double s = cos(phase + 2*M_PI * i / period);
            out[i] = (Cnum){s, 0};
        }
}
void sqr_wave(int samples, double cycles, double phase, Cnum *out) {
        double period = samples / cycles;
        for (int i = 0; i < samples; i++) {
            double s = phase + (i% (int)period < period/2 ?
                                      i % (int)period
                                    : period/2 - (i % (int)period)
                                );
            out[i] = (Cnum){s, 0};
        }
}
