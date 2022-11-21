#include "other_func.hpp"

#include <stdio.h>
#include <string.h>
#include <errno.h>

void
print_error_line (const char strerr[],
                  const char name_file[],
                  unsigned line) {
    fprintf (stderr, "Failed in %s, error: %s\nLINE: %s: %d\n",
                     strerr, strerror (errno), name_file, line);
}

void print_interest_args (int argc,
                           char* argv[]) {
    if (argc != 1) {
        printf ("argc: %d\n", argc);
        for (int i = 0; i < argc; ++i) {
            printf ("argv[%d]: %s\n", i, argv[i]);
        }
    }
}