#include "fft.h"
#include "util.h"
#include <stddef.h>
#include <mpi.h>

#define MASTER 0

void improper_usage_abort() {
    // Print out the proper usage of the program and exit
    puts ("Usage : fft -[cs] samples cycles [phase]");
    puts ("        fft -f file_name samples");
    exit (1);
}

Cnum *parse_args(int argc, char **argv) {
    Cnum *input;
    if (argc > 2) {
        if (argv[1][0] == '-') {
            if (!argv[1][1]) improper_usage_abort();
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
                else improper_usage_abort();
            }
        } else {
            if (argc != 3 && argc != 4) improper_usage_abort();
            double cycles = 0;
            sscanf(argv[1], "%d", &samples);
            sscanf(argv[2], "%lf", &cycles);
            input = (Cnum *) malloc(sizeof(Cnum)*samples);
            double phase = 0;
            if (argc == 4) sscanf(argv[3], "%lf", &phase);
            sin_wave(samples, cycles, phase*M_PI, input);
        }
    }
}

int main(int argc, char **argv) {
    int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //Create MPI_Datatype for complex numbers
    const int nitems = 2;
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype mpi_Cnum;
    MPI_Aint offsets[2];
    offsets[0] = offsetof(Cnum, real);
    offsets[1] = offsetof(Cnum, imag);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_Cnum);
    MPI_Type_commit(&mpi_Cnum);

    int depth = 0; // Tells how many times we have split the data.
    int length; // In master this is samples, the slaves must be told the length.

    if (myrank == MASTER) {
        // TAKING INPUT
        // init args, format :  fft -f file_name num_entries
        //                      fft samples cycles [phase] <- sin, phase in
        //                          increments of pi
        //                      fft -c samples cycles [phase] <- cos
        //                      fft -s samples cycles [phase] <- sqr
        Cnum *input;
        int samples = 0;
        if (argc > 2) {
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
        }
        else {
            usage:
                free (input);
                puts ("Usage : fft -[cs] samples cycles [phase]");
                puts ("        fft -f file_name samples");
                exit (1);
        }

        //DOING STUFF
        length = samples;
        depth = 1;
        Cnum *data = input;
        int target = myrank + pow(2, depth-1);
        while (target < nprocs && // Still have processors to send to
               length % 2 == 0 && length > 2) { // Still have work to distribute
            length /= 2;
            Cnum *evenData = malloc(sizeof(Cnum)*length);
            Cnum *oddData = malloc(sizeof(Cnum)*length);
            for (int i = 0; i < length; ++i) {
                evenData[i] = data[2*i];
                oddData[i] = data[2*i + 1];
            }
            data = evenData;
            // consider changing tag
            printf("depth: %d, master, target is %d\n", depth, target);
            MPI_Send(&length, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
            MPI_Send(&depth, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
            MPI_Send(oddData, length, mpi_Cnum, target, 0, MPI_COMM_WORLD);
            depth += 1;
            target = myrank + pow(2, depth-1);
        }
        Cnum *output = (Cnum *) malloc(sizeof(Cnum)*length);
        output = fft(length, data);
        Cnum *x;
        for (int i = 0; i < length; i++) {
            x = &output[i];
            printf("%f + %f i, 0\n", x->real, x->imag);
        }

        // wait for all finished
        int dummy_data;
        MPI_Status status;
        for (int i = 1; i < nprocs; ++i) {
            MPI_Recv(&dummy_data, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            printf("Master received signal from slave %d\n", i);
        }
        printf("Master finished\n");
    }
    else { //We are slaves!
        printf("Slave %d starting\n", myrank);
        MPI_Status status;
        // consider changing tag
        MPI_Recv(&length, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("Slave %d received length\n", myrank);
        MPI_Recv(&depth, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("Slave %d received depth\n", myrank);
        Cnum *data = (Cnum *) malloc(sizeof(Cnum)*length);
        MPI_Recv(data, length, mpi_Cnum, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("Slave %d received data\n", myrank);
        depth += 1;
        int target = myrank + pow(2, depth-1);
        while (target < nprocs && // Still have processors to send to
               length % 2 == 0 && length > 2) { // Still have work to distribute
            length /= 2;
            Cnum *evenData = malloc(sizeof(Cnum)*length);
            Cnum *oddData = malloc(sizeof(Cnum)*length);
            for (int i = 0; i < length; ++i) {
                evenData[i] = data[2*i];
                oddData[i] = data[2*i + 1];
            }
            data = evenData;
            // consider changing tag
            printf("depth: %d, slave %d, target is %d\n", depth, myrank, target);
            MPI_Send(&length, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
            MPI_Send(&depth, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
            MPI_Send(oddData, length, mpi_Cnum, target, 0, MPI_COMM_WORLD);
            depth += 1;
            target = myrank + pow(2, depth-1);
        }
        Cnum *output = (Cnum *) malloc(sizeof(Cnum)*length);
        output = fft(length, data);
        Cnum *x;
        for (int i = 0; i < length; i++) {
            x = &output[i];
            printf("%f + %f i, %d\n", x->real, x->imag, myrank);
        }
        printf("Slave %d finished\n", myrank);
        int dummy_data = 1;
        MPI_Send(&dummy_data, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
