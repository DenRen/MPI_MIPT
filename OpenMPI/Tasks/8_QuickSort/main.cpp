#include <omp.h>
#include <cstdio>
#include <unistd.h>
#include <vector>
#include <random>
#include <algorithm>

void swap(int* lhs, int* rhs) noexcept
{
    int tmp = *rhs;
    *rhs = *lhs;
    *lhs = tmp;
}

// [first, last]
void Qsort(int* arr, int first, int last)
{
    if (first >= last)
        return;

    const int ref_pos = (last + first) / 2;
    const int ref_elem = arr[ref_pos];

    int l = first, r = last;
    do {
        while (arr[l] < ref_elem) ++l;
        while (arr[r] > ref_elem) --r;
        if (l <= r)
            swap(arr + l++, arr + r--);
    } while (l <= r);

    if (last - first >= 10'000)
    {
        #pragma omp parallel shared(arr) firstprivate(first, r, l, last)
        {
            #pragma omp single nowait
            {
                #pragma omp task
                    Qsort(arr, first, r);
                #pragma omp task
                    Qsort(arr, l, last);
            }
        }
    }
    else
    {
        Qsort(arr, first, r);
        Qsort(arr, l, last);
    }
}

// #define ENABLE_RES_CHECK

int main()
{
    const int size = 100'000;

    std::mt19937_64 mers{15+0*std::random_device{}()};
    std::uniform_int_distribution<int> dist(0, 1000);

    std::vector<int> vec(size);

    for (int i = 0; i < 1000; ++i)
    {
    #ifdef ENABLE_RES_CHECK
        for (auto& val : vec)
            val = dist(mers);

        auto copy_vec = vec;

        Qsort(copy_vec.data(), 0, copy_vec.size() - 1);
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
            Qsort(vec.data(), 0, vec.size() - 1);
    #endif
    }
}
