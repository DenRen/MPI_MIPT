#include <stdio.h>
#include <threads.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>
#include <pthread.h>
#include <string.h>

// #define NDEBUG

#define PRINT_ERR(err_msg) \
    do {\
        fprintf (stderr, "Error in %s:%s:%d: %s\n", __FILE__, __PRETTY_FUNCTION__, \
                                                   __LINE__, ((err_msg)));\
        if (errno) {\
            fprintf (stderr, "errno: %s\n", strerror (errno));\
        }\
    } while (0)

#ifndef NDEBUG
    #define ASSERT(cond)\
        do {\
            if (!(cond)) {\
                PRINT_ERR (#cond);\
            }\
        } while (0)
#elif
    #define ASSERT(cond)
#endif

unsigned
read_uint (char* str, int* err) {
    ASSERT (str != NULL);
    ASSERT (err != NULL);

    char* endptr = NULL;
    errno = 0;
    long val = strtol (str, &endptr, 10);

    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
        || (errno != 0 && val == 0)) {
    #ifndef NDEBUG
        PRINT_ERR ("strtol");
    #endif
        *err = 1;
        return 0;
    }

    if (endptr == str || val < 0 || val >= UINT_MAX) {
    #ifndef NDEBUG
        if (endptr == str) {
            PRINT_ERR ("Error format");
        } else {
            PRINT_ERR ("Value is not unsigned");
        }
        *err = 1;
        return 0;
    #endif
    }

    return val;
}

double
read_double (char* str, int* err) {
    ASSERT (str != NULL);
    ASSERT (err != NULL);

    char* endptr = NULL;
    errno = 0;
    double val = strtod (str, &endptr);

    if ((errno == ERANGE && (val == HUGE_VAL || val == -HUGE_VAL))
        || (errno != 0 && val == 0)) {
    #ifndef NDEBUG
        PRINT_ERR ("strtol");
    #endif
        *err = 1;
        return 0;
    }

    if (endptr == str) {
    #ifndef NDEBUG
        if (endptr == str) {
            PRINT_ERR ("Error format");
        }
        *err = 1;
        return 0;
    #endif
    }

    return val;
}

unsigned
read_number_threads (char* str) {
    int err = 0;
    unsigned num_threads = read_uint (str, &err);
    if (err != 0 || num_threads == 0) {
        printf ("Please, enter the positive number of threads\n");
        exit (EXIT_FAILURE);
    }

    return num_threads;
}

unsigned
read_number_sections (char* str) {
    int err = 0;
    unsigned num_sections = read_uint (str, &err);
    if (err != 0 || num_sections == 0) {
        printf ("Please, enter the positive number of sections\n");
        exit (EXIT_FAILURE);
    }

    return num_sections;
}

double
read_eps (char* str) {
    int err = 0;
    double eps = read_double (str, &err);
    if (err != 0 || eps <= 0) {
        printf ("Please, enter the positive eps\n");
        exit (EXIT_FAILURE);
    }

    return eps;
}

double
calc_integral (unsigned num_threads,
               unsigned num_sections,
               double a,
               double b,
               int* err);

int main (int argc, char* argv[]) {
    unsigned num_threads = 0;
    double eps = 1e-8;
    double a = 1e-4, b = 1;
    unsigned num_sections = 1u << 31;

    if (argc >= 2) {
        num_threads = read_number_threads (argv[1]);
    }
    if (argc == 3) {
        eps = read_eps (argv[2]);
        // num_sections = read_number_sections (argv[2]);
    }

    if (argc == 1 || argc > 3) {
        printf ("Please, enter the number of threads and "
                "(optional) number of secontions for integrate\n");
        exit (EXIT_SUCCESS);
    }

    int err = 0;
    num_sections = (b - a) / eps;
    double res = calc_integral (num_threads, num_sections, a, b, &err);
    if (err) {
        PRINT_ERR ("calc_integral");
        exit (EXIT_FAILURE);
    }

    printf ("%.16g\n", res);

    return 0;
}

struct task_t {
    double a, dx;
    unsigned num_sections;
    double res;
};

void
print_task (const struct task_t* task) {
    printf ("a: %g\n", task->a);
    printf ("dx: %g\n", task->dx);
    printf ("num_sec: %u\n", task->num_sections);
}

void*
kernel_calc_integral (struct task_t* task) {
    double x = task->a, dx = task->dx, acc = 0;
    unsigned num_sections = task->num_sections;

    while (num_sections--) {
        acc += sin (x) * dx / (x * x) ;
        // acc += sin (1 / x) * dx;
        x += dx;
    }

    task->res = acc;

    return NULL;
}

double
calc_integral (unsigned num_threads,
               unsigned num_sections,
               double a,
               double b,
               int* err_) {
    ASSERT (a > 0);
    ASSERT (b >= a);
    ASSERT (num_threads > 0);
    ASSERT (num_sections > 0);
    ASSERT (err_ != NULL);

    // printf ("%g %g\n", a, b);
    double tmp = a;
    a = 1 / b;
    b = 1 / tmp;

    // printf ("%g %g\n", a, b);

    unsigned num_sect_per_thread = num_sections / num_threads +
                                   (num_sections % num_threads != 0);

    struct task_t tasks[num_threads];
    tasks[0].a = a;
    tasks[0].dx = (b - a) / num_sections;
    tasks[0].num_sections = num_sect_per_thread;

    pthread_t ids[num_threads - 1];
    for (unsigned i = 0; i + 1 < num_threads; ++i) {
        // print_task (&tasks[i]);
        int err = pthread_create (&ids[i], NULL,
                                  (void * (*)(void *)) kernel_calc_integral, &tasks[i]);
        if (err) {
            PRINT_ERR ("pthread_create");
            *err_ = 1;
            return 0;
        }

        tasks[i + 1].dx = tasks[i].dx;
        tasks[i + 1].a = tasks[i].a + tasks[i + 1].dx * num_sect_per_thread;
        tasks[i + 1].num_sections = tasks[i].num_sections;
    }

    tasks[num_threads - 1].num_sections = num_sections - (num_threads - 1) * num_sect_per_thread;

    // print_task (&tasks[num_threads - 1]);
    kernel_calc_integral (&tasks[num_threads - 1]);

    double res = tasks[num_threads - 1].res;
    // printf ("%d: %g\n", num_threads - 1, tasks[num_threads - 1].res);
    for (unsigned i = 0; i + 1 < num_threads; ++i) {
        int err = pthread_join (ids[i], NULL);
        if (err) {
            PRINT_ERR ("pthread_join");
            *err_ = 1;
            return 0;
        }

        // printf ("%d: %g\n", i, tasks[i].res);

        res += tasks[i].res;
    }

    return res;
}