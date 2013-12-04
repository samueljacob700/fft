#include "fft.h"
#include "util.h"
#include <stddef.h>
#include <mpi.h>
#include <string.h>

#define MASTER 0

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
        int length = samples;
        int depth = 1;
        Cnum *data = input;
        int target = myrank + nprocs / pow(2, depth);
        Cnum *evenData;
        Cnum *oddData;
        while (target > myrank && // Still have processors to send to
               length % 2 == 0 && length > 2) { // Still have work to distribute
            length /= 2;
            evenData = (Cnum *) malloc(sizeof(Cnum)*length);
            oddData = (Cnum *) malloc(sizeof(Cnum)*length);
            for (int i = 0; i < length; ++i) {
                evenData[i] = data[2*i];
                oddData[i] = data[2*i + 1];
            }
            data = evenData;
            // consider changing tag
            // printf("depth: %d, master, target is %d\n", depth, target);
            MPI_Send(&length, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
            MPI_Send(&depth, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
            MPI_Send(oddData, length, mpi_Cnum, target, 0, MPI_COMM_WORLD);
            depth += 1;
            target = myrank + nprocs / pow(2, depth);
        }
        Cnum *output = (Cnum *) malloc(sizeof(Cnum)*length);
        output = fft(length, data);
        // printf("master completed fft\n");

        MPI_Status status;

        int available = myrank + nprocs / pow(2, depth-1);
        Cnum *combined;
        free(evenData);
        free(oddData);
        // printf(" ------ master : available is %d\n", available);
        while (available < nprocs) {
            oddData = (Cnum *) malloc(sizeof(Cnum)*length);
            // printf(" ------ master : about to receive %d data points from slave %d\n", length, (int)(myrank + nprocs / pow(2, depth-1)));
            MPI_Recv(oddData, length, mpi_Cnum, (int)(myrank + nprocs / pow(2, depth-1)), MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            // printf(" ------ master : received successfully\n");

            combined = fft_combine(2*length, output, oddData);
            // printf(" ------ master : successful fft_combine\n");
            length *= 2;
            free(output);
            free(oddData);
            output = combined;
            // printf(" ------ master : free'd output and oddData memory\n");

            depth -= 1;
            available = nprocs / pow(2, depth-1);
            // printf(" ------ master : available is %d\n", available);
        }
        // printf(" ------ master : done with while loop\n");
        printf("Master finished\n");

        Cnum *x;
        for (int i = 0; i < length; i++) {
            x = &output[i];
            printf("%.3f + %.3f i, 0\n", x->real, x->imag);
        }
        free(output);

    }
    else { //We are slaves!
        int length, depth;

        // printf("Slave %d starting\n", myrank);
        MPI_Status status;
        // consider changing tag
        MPI_Recv(&length, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        // printf("Slave %d received length: %d\n", myrank, length);
        MPI_Recv(&depth, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        // printf("Slave %d received depth\n", myrank);
        Cnum *data = (Cnum *) malloc(sizeof(Cnum)*length);
        MPI_Recv(data, length, mpi_Cnum, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        // printf("Slave %d received data\n", myrank);
        depth += 1;
        int target = myrank + nprocs / pow(2, depth);
        // printf("Slave %d, target is %d\n", myrank, target);
        Cnum *evenData;
        Cnum *oddData;
        while (target > myrank && // Still have processors to send to
               length % 2 == 0 && length > 2) { // Still have work to distribute
            length /= 2;
            evenData = (Cnum *) malloc(sizeof(Cnum)*length);
            oddData = (Cnum *) malloc(sizeof(Cnum)*length);
            for (int i = 0; i < length; ++i) {
                evenData[i] = data[2*i];
                oddData[i] = data[2*i + 1];
            }
            free(data); 
            data = (Cnum *) malloc(sizeof(Cnum)*length);
            memcpy(data, evenData, sizeof(Cnum)*length);
            // printf("depth: %d, slave %d, target is %d\n", depth, myrank, target);
            MPI_Send(&length, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
            MPI_Send(&depth, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
            MPI_Send(oddData, length, mpi_Cnum, target, 0, MPI_COMM_WORLD);
            depth += 1;
            target = myrank + nprocs / pow(2, depth);

            free(oddData);
            free(evenData);
        }
        // printf("Slave %d after while loop(sending)\n", myrank);
        Cnum *output = (Cnum *) malloc(sizeof(Cnum)*length);
        output = fft(length, data);
        // printf("Slave %d after first fft\n", myrank);

        // Cnum *x;
        // for (int i = 0; i < length; i++) {
        //     x = &output[i];
        //     printf("%f + %f i, %d\n", x->real, x->imag, myrank);
        // }

        int available = nprocs / pow(2, depth-1);
        int source = status.MPI_SOURCE;
        Cnum *combined = output;
        // printf(" ------ slave %d : about to enter while loop with available as %d\n", myrank, available);
        while (myrank - available != source) {
            // printf(" ------ slave %d : inside while loop\n", myrank);
            // combined = (Cnum *) malloc(sizeof(Cnum)*length*2);
            oddData = ((Cnum *) malloc(sizeof(Cnum)*length));
            // printf(" ------ slave %d : about to receive %d data points from %d\n", myrank, length, (int)(myrank + nprocs / pow(2, depth-1)));
            MPI_Recv(oddData, length, mpi_Cnum, (int)(myrank + nprocs / pow(2, depth-1)), MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            // printf(" ------ slave %d : successfully received data\n", myrank);
            combined = fft_combine(2*length, output, oddData);
            // printf(" ------ slave %d : successful fft_combine\n", myrank);
            depth -= 1;
            available = nprocs / pow(2, depth-1);
            length *= 2;
            free(output);
            free(oddData);
            output = combined;
        }
        // printf(" ------ slave %d : outside of while loop\n", myrank);
        // printf(" ------ slave %d : about to send %d data points to source %d\n", myrank, length, source);

        // should combined be output? supposedly does not matter
        MPI_Send(combined, length, mpi_Cnum, source, 0, MPI_COMM_WORLD);
        // printf(" ------ slave %d : sent to source %d\n", myrank, source);
        free(output);

        // printf("Slave %d finished\n", myrank);
    }

    MPI_Finalize();
    return 0;
}
