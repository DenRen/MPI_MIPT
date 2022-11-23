#include <omp.h>
#include <cstdio>
#include <unistd.h>
#include <vector>
#include <random>
#include <algorithm>

void swap(int *lhs, int *rhs)
{
    int tmp = *rhs;
    *rhs = *lhs;
    *lhs = tmp;
}

// [first, last]
void Qsort(int *arr, int first, int last)
{
    if (first >= last)
        return;

    const int ref_pos = (last + first) / 2;
    const int ref_elem = arr[ref_pos];

    int l = first, r = last;
    do
    {
        while (arr[l] < ref_elem)
            ++l;
        while (arr[r] > ref_elem)
            --r;
        if (l <= r)
            swap(arr + l++, arr + r--);
    } while (l <= r);

    #pragma omp task shared(arr) if (last - first > 1000)
        Qsort(arr, first, r);

    #pragma omp task shared(arr) if (last - first > 1000)
        Qsort(arr, l, last);
}

// #define RANDOM_ARRAY
// #define STDIO_ARRAY

// #define ENABLE_RES_CHECK

static int main_random_array()
{
    const int size = 100'000;
    const int num_repeats = 1000;

#ifdef ENABLE_RES_CHECK
    std::mt19937_64 mers{15 + 0 * std::random_device{}()};
    std::uniform_int_distribution<int> dist(0, 1000);
#endif

    omp_set_dynamic(0);
    std::vector<int> vec(size);

    double begin = omp_get_wtime();
    for (int i_repeat = 0; i_repeat < num_repeats; ++i_repeat)
    {
    #ifdef ENABLE_RES_CHECK
        for (auto &val : vec)
            val = dist(mers);

        auto copy_vec = vec;

        #pragma omp parallel
        {
            #pragma omp single
            {
                Qsort(copy_vec.data(), 0, copy_vec.size() - 1);
            }
        }
        std::sort(vec.begin(), vec.end());

        if (vec != copy_vec)
        {
            printf("Error:\n");
            for (int j = 0; j < vec.size(); ++j)
                if (vec[j] != copy_vec[j])
                    printf("%3d: %4d %4d\n", j, vec[j], copy_vec[j]);

            return -1;
        }
    #else
        #pragma omp parallel
        {
            #pragma omp single
            {
                Qsort(vec.data(), 0, vec.size() - 1);
            }
        }
    #endif
    }
    double end = omp_get_wtime();
    printf("time: %g\n", 1e3 * (end - begin) / num_repeats);

    return 0;
}

// Format:
// Input:  "5 9 8 7 6 5"
// Output:   "5 6 7 8 9"
static int main_stdio()
{
    int size = 0;
    if (scanf("%d", &size) != 1)
    {
        perror("Incorrect size");
        return -1;
    }

    std::vector<int> vec(size);
    for (int i = 0; i < size; ++i)
    {
        if (scanf("%d", &vec[i]) != 1)
        {
            perror("Incorrect input data");
            return -1;
        }
    }

    #pragma omp parallel
    {
        #pragma omp single
        {
            Qsort(vec.data(), 0, vec.size() - 1);
        }
    }

    for (const auto value : vec)
        printf("%d ", value);

    return 0;
}

int main()
{
#ifdef RANDOM_ARRAY
    main_random_array();
#endif

#ifdef STDIO_ARRAY
    main_stdio();
#endif

}
