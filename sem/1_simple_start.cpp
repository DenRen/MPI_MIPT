#include <stdio.h>
#include <mpi.h>

#include "other_func.hpp"

int main (int argc, char* argv[]) {
    print_interest_args (argc, argv);

    printf ("Before MPI_INIT\n");
    CHECK_MPI_ERR (MPI_Init (&argc, &argv));

    printf ("Parallel sect\n");
    CHECK_MPI_ERR (MPI_Finalize ());

    printf ("After MPI_FINALIZE\n");
}