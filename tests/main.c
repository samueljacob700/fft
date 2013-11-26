#include "../fft.h"
#include "../util.h"

int
compare (int N, Cnum *d, Cnum *f)
{
     for (int i = 0; i < N; i++) {
         double x = mag2(&d[i]);
         double y = mag2(&f[i]);
         /*printf("%f  ==  %f\n", x, y);*/
         if (fabs(x - y) > 0.000001) {
             printf("%f\n", fabs(x-y));     
             return 0;
         }
     }
     return 1;
}

int
main (int argc, char **argv)
{
    int N = 500;
    double cycles = 25;
    double phase = M_PI / 2;

    Cnum *test_file_in = read_file_in("testdata.txt");
    Cnum d_file_out[8];
    Cnum *f_file_out = (Cnum *) malloc(sizeof(Cnum)*8);

    Cnum *temp_in = (Cnum *) malloc(sizeof(Cnum)*N);

    Cnum *d_out = (Cnum *) malloc(sizeof(Cnum)*N);
    Cnum *f_out = (Cnum *) malloc(sizeof(Cnum)*N);

    int passed = 0;
    int temp;

    dft(8, test_file_in, d_file_out);
    f_file_out = fft(8, test_file_in);

    printf("Testing file input... ");
    // check file output
    if ((temp = compare(N, d_file_out, f_file_out)))
        printf("passing!\n");
    else
        printf("FAILED\n");

    passed += temp;

    sin_wave(N, cycles, phase, temp_in);

    dft(N, temp_in, d_out);
    f_out = fft(N, temp_in);

    printf("Testing sin wave input... ");
    // check sin output
    if ((temp = compare(N, d_out, f_out)))
        printf("passing!\n");
    else
        printf("FAILED\n");

    passed += temp;

    cos_wave(N, cycles, phase, temp_in);
    dft(N, temp_in, d_out);
    f_out = fft(N, temp_in);

    printf("Testing cos wave input... ");
    // check cos output
    if ((temp = compare(N, d_out, f_out)))
        printf("passing!\n");
    else
        printf("FAILED\n");

    passed += temp;

    sqr_wave(N, cycles, phase, temp_in);
    dft(N, temp_in, d_out);
    f_out = fft(N, temp_in);

    printf("Testing square wave input... ");
    // check sqr output
    if ((temp = compare(N, d_out, f_out)))
        printf("passing!\n");
    else
        printf("FAILED\n");

    passed += temp;

    printf("\n%d / 4 tests passed.\n", passed);

    free(d_out);
    free(f_out);
    free(test_file_in);
    free(f_file_out);
    free(temp_in);
}
