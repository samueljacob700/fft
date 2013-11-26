#include "util.h"

double mag2 (Cnum *n) {
    return (n->real*n->real + n->imag*n->imag);
}
double magn (Cnum *n) {
    return sqrt(mag2(n));
}

Cnum *read_file_in(char *path) {

    FILE *fr;
    int count = 0;
    char line[80];
    Cnum holder;
    Cnum* out  = NULL; 
    Cnum* outx = NULL;

    fr = fopen (path, "rt");
    while(fgets(line, 80, fr) != NULL) {
        double r, i;
        sscanf(line, "%lf, %lf\n", &r, &i);
        count++;

        outx = (Cnum*) realloc (out, count * sizeof(Cnum));

        if (outx != NULL) {
            out = outx;
            out[count-1] = (Cnum) {r, i};
        } else {
            free(out);
            puts("Error (re)allocating memory");
            exit(1);
        }
    }
    fclose(fr);
    free(out);

    return out;
}
