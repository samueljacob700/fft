#include "fft.h"
#include "util.h"

int main(int argc, char **argv) {
    // init args, format :  fft -f file_name num_entries
    //                      fft samples cycles [phase] <- sin, phase in
    //                          increments of pi
    //                      fft -c samples cycles [phase] <- cos
    //                      fft -s samples cycles [phase] <- sqr
    Cnum *input;
    if (argc > 2) {
        int samples = 0;
        if (argv[1][0] == '-') {
            if (!argv[1][1]) goto usage;
            if (argv[1][1] == 'f') {
                if (argc == 4) {
                    sscanf(argv[3], "%d", &samples);
                    input = (Cnum *) malloc(sizeof(Cnum)*samples);
                    input = read_file_in(argv[2]);
                } else goto usage;
            } else {
                double cycles = 0;
                sscanf(argv[2], "%d", &samples);
                sscanf(argv[3], "%lf", &cycles);
                double phase = 0;
                input = (Cnum *) malloc(sizeof(Cnum)*samples);
                if (argc == 5) sscanf(argv[4], "%lf", &phase);
                if (argv[1][1] == 'c') cos_wave(samples, cycles, phase*M_PI, input);
                if (argv[1][1] == 's') sqr_wave(samples, cycles, phase*M_PI, input);
                else goto usage;
            }
        } else {
            if (argc != 3 && argc != 4) goto usage;
            double cycles = 0;
            sscanf(argv[1], "%d", &samples);
            sscanf(argv[2], "%lf", &cycles);
            input = (Cnum *) malloc(sizeof(Cnum)*samples);
            double phase = 0;
            if (argc == 4) sscanf(argv[3], "%lf", &phase);
            sin_wave(samples, cycles, phase*M_PI, input);
        }
        Cnum *output = (Cnum *) malloc(sizeof(Cnum)*samples);
        /*Cnum *dout = (Cnum *) malloc(sizeof(Cnum)*samples);*/

        Cnum *even = (Cnum *) malloc(sizeof(Cnum)*samples/2);
        Cnum *odd  = (Cnum *) malloc(sizeof(Cnum)*samples/2);
        for (int i = 0; i < samples/2; i++) {
            even[i] = input[2*i];
            odd[i]  = input[2*i+1];
        }
        // output = fft(samples, input);
        output = fft_combine(samples, fft(samples/2, even), fft(samples/2, odd));
        Cnum *x;
        printf("Test example for two nodes:\n");
        for (int i = 0; i < samples; i++) {
            /*double y = magn(&output[i])/samples;*/
            x = &output[i];
            printf("%f + %f i\n", x->real, x->imag);
        }
        printf("\n");
        Cnum *even_even = (Cnum *) malloc(sizeof(Cnum)*samples/4);
        Cnum *even_odd  = (Cnum *) malloc(sizeof(Cnum)*samples/4);
        for (int i = 0; i < samples/4; i++) {
            even_even[i] = even[2*i];
            even_odd[i]  = even[2*i+1];
        }
        Cnum *odd_even = (Cnum *) malloc(sizeof(Cnum)*samples/4);
        Cnum *odd_odd  = (Cnum *) malloc(sizeof(Cnum)*samples/4);
        for (int i = 0; i < samples/4; i++) {
            odd_even[i] = odd[2*i];
            odd_odd[i]  = odd[2*i+1];
        }
        printf("Test example for four nodes:\n");
        output = fft_combine(samples, fft_combine(samples/2, fft(samples/4, even_even), fft(samples/4, even_odd)), fft_combine(samples/2, fft(samples/4, odd_even), fft(samples/4, odd_odd)));
        for (int i = 0; i < samples; i++) {
            /*double y = magn(&output[i])/samples;*/
            x = &output[i];
            printf("%f + %f i\n", x->real, x->imag);
        }
        return 0;
    }


    usage:
        free (input);
        puts ("Usage : fft -[cs] samples cycles [phase]");
        puts ("        fft -f file_name samples");
        exit (1);
}
