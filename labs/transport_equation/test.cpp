#include <gtest/gtest.h>
#include "lib.hpp"

TEST (COPY_ROW_2_COL, STATIC) {
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

TEST (COPY_COL_2_ROW, STATIC) {
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

TEST (COPY_ROW_2_RECT, STATIC) {
    int u_buf_x_size = 15;
    int u_buf_t_size = 20;

    int row_buf_size = 16;

    std::vector <double> u_buf (u_buf_x_size * u_buf_t_size);
    std::vector <double> row_buf (row_buf_size);

    for (int i = 0; i < row_buf.size (); ++i) {
        row_buf[i] = i + 1;
    }

    // print_2d_array (u_buf, u_buf_x_size);
    
    double* dest = u_buf.data () + u_buf_x_size + 2;
    const int width = 8;

    copy_row_2_rect (dest, row_buf.data (), width, u_buf_x_size, row_buf.size ());
    
    // print_2d_array (u_buf, u_buf_x_size);

    for (int i_row = 0; i_row < row_buf.size () / width; ++i_row) {
        double* begin = dest + i_row * u_buf_x_size;
        for (int i_col = 0; i_col < width; ++i_col) {
            ASSERT_EQ (begin[i_col], 1 + i_col + i_row * width);
        }
    }
}

TEST (CALC_END_INDEX, STATIC) {
    ASSERT_EQ (calc_end_index ( 0, 4,  5), 4);
    ASSERT_EQ (calc_end_index (10, 4, 13), 13);
    ASSERT_EQ (calc_end_index (10, 4, 14), 14);
    ASSERT_EQ (calc_end_index (10, 4, 15), 14);
    ASSERT_EQ (calc_end_index (10, 4, 5),  5);
}

TEST (REMOVE_LEFT_DOWN_BOUND, STATIC) {
    std::vector <double> v = {
         1,  2,  3,  4,  5,
         6,  7,  8,  9, 10,
        11, 12, 13, 14, 15
    };

    remove_left_down_bound (v, 5);

    std::vector <double> v_ref = {
        7,  8,  9, 10,
        12, 13, 14, 15,
    };

    ASSERT_EQ (v, v_ref);
}