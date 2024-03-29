#pragma once

#include <cstdlib>
#include <algorithm>

template <typename T>
void
copy_row_2_col (T* dst_col,
                const T* src_row,
                std::size_t row_size,
                std::size_t num)
{
    for (std::size_t i_row = 0; i_row < num; ++i_row) {
        *dst_col = src_row[i_row];
        dst_col += row_size;
    }
}

template <typename T>
void
copy_col_2_row (T* dst_row,
                const T* src_col,
                std::size_t row_size,
                std::size_t num)
{
    for (std::size_t i_row = 0; i_row < num; ++i_row) {
        dst_row[i_row] = *src_col;
        src_col += row_size;
    }
}
/*
                #########################
                #                       #
                #                       #
                #        #####          #
                #        #   #          #
                #        #   #          #
                #        #   #          #
dst_rect -------#------->#####          #
                #        width          #
                #                       #
                #                       #
                #########################
                <------ dst_x_size ----->
*/
template <typename T>
void
copy_row_2_rect (T* dst_rect,
                 const T* src_row,
                 std::size_t width,
                 std::size_t dst_x_size,
                 std::size_t row_size)
{
    const std::size_t num_row = row_size / width;
    for (std::size_t i_row = 0; i_row < num_row; ++i_row) {
        const T* src_begin = src_row + i_row * width;
        T* dst_begin = dst_rect + i_row * dst_x_size;
        std::copy_n (src_begin, width, dst_begin);
    }
}

template <typename T>
void
remove_left_down_bound (std::vector <T>& vec,
                        std::size_t row_size)
{
    T* dst = vec.data ();
    const std::size_t num_rows = vec.size () / row_size;
    for (std::size_t i_row = 1; i_row < num_rows; ++i_row) {
        T* src = vec.data () + i_row * row_size + 1;
        std::copy_n (src, row_size - 1, dst);
        dst += row_size - 1;
    }
    vec.resize (dst - vec.data ());
}