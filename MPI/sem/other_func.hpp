#pragma once

void
print_error_line (const char strerr[],
                  const char name_file[],
                  unsigned line);

#define CHECK_MPI_ERR(err_code)                                 \
    do {                                                        \
        if (err_code != MPI_SUCCESS) {                          \
            print_error_line (#err_code, __FILE__, __LINE__);   \
            return -1;                                          \
        }                                                       \
    } while (0)


void print_interest_args (int argc,
                           char* argv[]);