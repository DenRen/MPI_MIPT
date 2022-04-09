#include <gtest/gtest.h>
#include "lib.hpp"

template <typename T>
void
print_2d_array (const T* begin,
                std::size_t row_size,
                std::size_t full_size)
{
    std::size_t col_size = full_size / row_size;

    for (std::size_t t = col_size - 1; t + 1 > 0; --t) {
        std::cout << begin[t * row_size];
        for (std::size_t x = 1; x < row_size; ++x) {
            std::cout << ' ' << begin[x + t * row_size];
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

TEST (ROW_COL, STATIC) {
    int x_size = 10, t_size = 15;
    std::vector <double> v (x_size * t_size);

    for (int i = 0; i < x_size; ++i) {
        v[i] = i + 1;
    }

    // print_2d_array (v.data (), x_size, v.size ());
    copy_row_2_col (v.data () + 2 * x_size + 2, v.data (), x_size, x_size);
    // print_2d_array (v.data (), x_size, v.size ());

    for (int i = 1; i <= x_size; ++i) {
        ASSERT_EQ (v[i * x_size + x_size + 2], i);
    }
}

TEST (COL_ROW, STATIC) {
    int x_size = 10, t_size = 15;
    std::vector <double> v (x_size * t_size);

    for (int i = 1; i <= x_size; ++i) {
        v[i * x_size + x_size + 2] = i;
    }
    
    // print_2d_array (v.data (), x_size, v.size ());
    copy_col_2_row (v.data (), v.data () + 2 * x_size + 2, x_size, x_size);
    // print_2d_array (v.data (), x_size, v.size ());

    for (int i = 0; i < x_size; ++i) {
        ASSERT_EQ (v[i], i + 1);
    }
}