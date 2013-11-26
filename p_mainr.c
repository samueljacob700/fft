#include "fft.h"
#include "util.h"
#include <stddef.h>
#include <mpi.h>

#define MASTER 0

Cnum* parallel_fft(int procs_available, int data_length, Cnum* data, int source, int myrank, MPI_Datatype type) {
    printf("Proc %d entering parallel_fft\n", myrank);
    if (procs_available < 2 || data_length <= 2) {
        if (source != myrank) {
            MPI_Send(fft(data_length, data), data_length, type, source, 0, MPI_COMM_WORLD);
            return NULL;
        }
        return fft(data_length, data);
    }
    int half_data_length = data_length / 2;
    Cnum *even_data = malloc(sizeof(Cnum)*half_data_length);
    Cnum *odd_data = malloc(sizeof(Cnum)*half_data_length);
    for (int i = 0; i < half_data_length; ++i) {
        even_data[i] = data[2*i];
        odd_data[i] = data[2*i + 1];
    }
    procs_available /= 2;
    int target = myrank + procs_available;
    // consider changing tag
    MPI_Send(&half_data_length, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
    MPI_Send(&procs_available, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
    MPI_Send(odd_data, half_data_length, type, target, 0, MPI_COMM_WORLD);
    // Slave will call the second recursive case
    even_data = parallel_fft(procs_available, half_data_length, even_data, source, myrank, type); // call recursively on this machine
    // Slave sends results back
    MPI_Status status;
    MPI_Recv(&odd_data, data_length, type, target, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    if (source != myrank) {
        printf("Proc %d sending back to %dn", myrank, source);
        MPI_Send(fft_combine(data_length, even_data, odd_data), data_length, type, source, 0, MPI_COMM_WORLD);
        return NULL;
    }
    return fft_combine(data_length, even_data, odd_data);
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
        int procs_available = nprocs; // Assume this is a power of 2
        int target;
        Cnum *data = input;
        Cnum *output = (Cnum *) malloc(sizeof(Cnum)*length);
        output = parallel_fft(procs_available, length, data, myrank, myrank, mpi_Cnum);
        Cnum *x;
        for (int i = 0; i < length; i++) {
            x = &output[i];
            printf("%f + %f i, 0\n", x->real, x->imag);
        }
    }
    else { //We are slaves!
        int length;
        int procs_available; // Assume this is a power of 2
        int target;
        MPI_Status status;
        // consider changing tag
        MPI_Recv(&length, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&procs_available, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        Cnum *data = (Cnum *) malloc(sizeof(Cnum)*length);
        MPI_Recv(data, length, mpi_Cnum, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        parallel_fft(procs_available, length, data, status.MPI_SOURCE, myrank, mpi_Cnum);
    }
    MPI_Finalize();
    return 0;
}
