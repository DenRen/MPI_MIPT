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
} // TEST (COPY_ROW_2_COL, STATIC)

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
} // TEST (COPY_COL_2_ROW, STATIC)

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
} // TEST (COPY_ROW_2_RECT, STATIC)

TEST (CALC_END_INDEX, STATIC) {
    ASSERT_EQ (calc_end_index ( 0, 4,  5), 4);
    ASSERT_EQ (calc_end_index (10, 4, 13), 13);
    ASSERT_EQ (calc_end_index (10, 4, 14), 14);
    ASSERT_EQ (calc_end_index (10, 4, 15), 14);
    ASSERT_EQ (calc_end_index (10, 4, 5),  5);
} // TEST (CALC_END_INDEX, STATIC)

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
} // TEST (REMOVE_LEFT_DOWN_BOUND, STATIC)

namespace spec_caise {

treq::map_manager_t mgr {31, 24,
                         0, 30, 0, 24,
                         5, 5};

void
test_zero_area (const treq::map_manager_t& mgr) {
    treq::area_params_t area = mgr.get_area_params (0);
    EXPECT_EQ (area.x_i_end, 31);
    EXPECT_EQ (area.t_i_end, 24);
    EXPECT_EQ (area.get_zero_chunk_x_size (), 7);
    EXPECT_EQ (area.get_zero_chunk_t_size (), 5);

    treq::rect_t& zero_rect = area.zero_rect;
    EXPECT_EQ (zero_rect.x_i_min, 0);
    EXPECT_EQ (zero_rect.t_i_min, 0);
    EXPECT_EQ (zero_rect.x_chunk_size, 7);
    EXPECT_EQ (zero_rect.t_chunk_size, 5);
    EXPECT_EQ (zero_rect.get_x_i_end (), 7);
    EXPECT_EQ (zero_rect.get_t_i_end (), 5);
    EXPECT_EQ (zero_rect.get_offset_rect (10), 0);
    EXPECT_EQ (zero_rect.get_offset_rect (10, zero_rect), 0);

    auto rect_right = area.get_right_rect (1);
    EXPECT_EQ (rect_right.x_i_min, 7);
    EXPECT_EQ (rect_right.t_i_min, 0);
    EXPECT_EQ (rect_right.x_chunk_size, 7);
    EXPECT_EQ (rect_right.t_chunk_size, 5);
    EXPECT_EQ (rect_right.get_x_i_end (), 14);
    EXPECT_EQ (rect_right.get_t_i_end (), 5);
    EXPECT_EQ (rect_right.get_offset_rect (10), 7);
    EXPECT_EQ (rect_right.get_offset_rect (10, zero_rect), 7);
    EXPECT_EQ (rect_right.get_offset_rect (10, rect_right), 0);

    auto rect_up = area.get_up_rect (1);
    EXPECT_EQ (rect_up.x_i_min, 0);
    EXPECT_EQ (rect_up.t_i_min, 5);
    EXPECT_EQ (rect_up.x_chunk_size, 7);
    EXPECT_EQ (rect_up.t_chunk_size, 5);
    EXPECT_EQ (rect_up.get_x_i_end (), 7);
    EXPECT_EQ (rect_up.get_t_i_end (), 10);
    EXPECT_EQ (rect_up.get_offset_rect (10), 10 * 5);
    EXPECT_EQ (rect_up.get_offset_rect (10, zero_rect), 10 * 5);
    EXPECT_EQ (rect_up.get_offset_rect (10, rect_up), 0);

    auto rect_right_last = area.get_right_rect (4);
    EXPECT_EQ (rect_right_last.x_i_min, 28);
    EXPECT_EQ (rect_right_last.t_i_min, 0);
    EXPECT_EQ (rect_right_last.x_chunk_size, 3);
    EXPECT_EQ (rect_right_last.t_chunk_size, 5);//
    EXPECT_EQ (rect_right_last.get_x_i_end (), 31);
    EXPECT_EQ (rect_right_last.get_t_i_end (), 5);//
    EXPECT_EQ (rect_right_last.get_offset_rect (10), 28);
    EXPECT_EQ (rect_right_last.get_offset_rect (10, zero_rect), 28);
    EXPECT_EQ (rect_right_last.get_offset_rect (10, rect_right_last), 0);

    auto rect_up_last = area.get_up_rect (4);
    EXPECT_EQ (rect_up_last.x_i_min, 0);
    EXPECT_EQ (rect_up_last.t_i_min, 20);
    EXPECT_EQ (rect_up_last.x_chunk_size, 7);
    EXPECT_EQ (rect_up_last.t_chunk_size, 4);
    EXPECT_EQ (rect_up_last.get_x_i_end (), 7);
    EXPECT_EQ (rect_up_last.get_t_i_end (), 24);
    EXPECT_EQ (rect_up_last.get_offset_rect (10), 10 * 20);
    EXPECT_EQ (rect_up_last.get_offset_rect (10, zero_rect), 10 * 20);
    EXPECT_EQ (rect_up_last.get_offset_rect (10, rect_up_last), 0);

    EXPECT_EQ (area.get_zero_and_right_rect_size (), 31 * 5);
    EXPECT_EQ (area.get_zero_and_right_rect_x_size (), 31);
    EXPECT_EQ (area.get_up_rect_size (), (23 - 4) * 7);
    EXPECT_EQ (area.get_up_rect_x_size (), 7);
} // void test_zero_area (const treq::map_manager_t& mgr)

void
test_first_area (const treq::map_manager_t& mgr) {
    treq::area_params_t area = mgr.get_area_params (1);
    EXPECT_EQ (area.x_i_end, 31);
    EXPECT_EQ (area.t_i_end, 24);
    EXPECT_EQ (area.get_zero_chunk_x_size (), 7);
    EXPECT_EQ (area.get_zero_chunk_t_size (), 5);

    treq::rect_t& zero_rect = area.zero_rect;
    EXPECT_EQ (zero_rect.x_i_min, 7);
    EXPECT_EQ (zero_rect.t_i_min, 5);
    EXPECT_EQ (zero_rect.x_chunk_size, 7);
    EXPECT_EQ (zero_rect.t_chunk_size, 5);
    EXPECT_EQ (zero_rect.get_x_i_end (), 14);
    EXPECT_EQ (zero_rect.get_t_i_end (), 10);
    EXPECT_EQ (zero_rect.get_offset_rect (10), 57);
    EXPECT_EQ (zero_rect.get_offset_rect (10, zero_rect), 0);

    auto rect_right = area.get_right_rect (1);
    EXPECT_EQ (rect_right.x_i_min, 14);
    EXPECT_EQ (rect_right.t_i_min, 5);
    EXPECT_EQ (rect_right.x_chunk_size, 7);
    EXPECT_EQ (rect_right.t_chunk_size, 5);
    EXPECT_EQ (rect_right.get_x_i_end (), 21);
    EXPECT_EQ (rect_right.get_t_i_end (), 10);
    EXPECT_EQ (rect_right.get_offset_rect (10), 50 + 14);
    EXPECT_EQ (rect_right.get_offset_rect (10, zero_rect), 7);
    EXPECT_EQ (rect_right.get_offset_rect (10, rect_right), 0);

    auto rect_up = area.get_up_rect (1);
    EXPECT_EQ (rect_up.x_i_min, 7);
    EXPECT_EQ (rect_up.t_i_min, 10);
    EXPECT_EQ (rect_up.x_chunk_size, 7);
    EXPECT_EQ (rect_up.t_chunk_size, 5);
    EXPECT_EQ (rect_up.get_x_i_end (), 14);
    EXPECT_EQ (rect_up.get_t_i_end (), 15);
    EXPECT_EQ (rect_up.get_offset_rect (10), 100 + 7);
    EXPECT_EQ (rect_up.get_offset_rect (10, zero_rect), 50);
    EXPECT_EQ (rect_up.get_offset_rect (10, rect_up), 0);

    auto rect_right_last = area.get_right_rect (3);
    EXPECT_EQ (rect_right_last.x_i_min, 28);
    EXPECT_EQ (rect_right_last.t_i_min, 5);
    EXPECT_EQ (rect_right_last.x_chunk_size, 3);
    EXPECT_EQ (rect_right_last.t_chunk_size, 5);
    EXPECT_EQ (rect_right_last.get_x_i_end (), 31);
    EXPECT_EQ (rect_right_last.get_t_i_end (), 10);
    EXPECT_EQ (rect_right_last.get_offset_rect (10), 50 + 28);
    EXPECT_EQ (rect_right_last.get_offset_rect (10, zero_rect), 21);
    EXPECT_EQ (rect_right_last.get_offset_rect (10, rect_right_last), 0);

    auto rect_up_last = area.get_up_rect (3);
    EXPECT_EQ (rect_up_last.x_i_min, 7);
    EXPECT_EQ (rect_up_last.t_i_min, 20);
    EXPECT_EQ (rect_up_last.x_chunk_size, 7);
    EXPECT_EQ (rect_up_last.t_chunk_size, 4);
    EXPECT_EQ (rect_up_last.get_x_i_end (), 14);
    EXPECT_EQ (rect_up_last.get_t_i_end (), 24);
    EXPECT_EQ (rect_up_last.get_offset_rect (10), 200 + 7);
    EXPECT_EQ (rect_up_last.get_offset_rect (10, zero_rect), 150);
    EXPECT_EQ (rect_up_last.get_offset_rect (10, rect_up_last), 0);

    EXPECT_EQ (area.get_zero_and_right_rect_size (), (31 - 7) * 5);
    EXPECT_EQ (area.get_zero_and_right_rect_x_size (), 31 - 7);
    EXPECT_EQ (area.get_up_rect_size (), (23 - 9) * 7);
    EXPECT_EQ (area.get_up_rect_x_size (), 7);
} // void test_first_area (const treq::map_manager_t& mgr)

void
test_last_area (const treq::map_manager_t& mgr) {
    treq::area_params_t area = mgr.get_area_params (4);
    EXPECT_EQ (area.x_i_end, 31);
    EXPECT_EQ (area.t_i_end, 24);
    EXPECT_EQ (area.get_zero_chunk_x_size (), 3);
    EXPECT_EQ (area.get_zero_chunk_t_size (), 4);

    treq::rect_t& zero_rect = area.zero_rect;
    EXPECT_EQ (zero_rect.x_i_min, 28);
    EXPECT_EQ (zero_rect.t_i_min, 20);
    EXPECT_EQ (zero_rect.x_chunk_size, 3);
    EXPECT_EQ (zero_rect.t_chunk_size, 4);
    EXPECT_EQ (zero_rect.get_x_i_end (), 31);
    EXPECT_EQ (zero_rect.get_t_i_end (), 24);
    EXPECT_EQ (zero_rect.get_offset_rect (100), 100 * 20 + 28);
    EXPECT_EQ (zero_rect.get_offset_rect (10, zero_rect), 0);

    EXPECT_EQ (area.get_up_rect_size (), 0);
    EXPECT_EQ (area.get_zero_and_right_rect_size (), 3 * 4);
    EXPECT_EQ (area.get_zero_and_right_rect_x_size (), 3);
} // void test_last_area (const treq::map_manager_t& mgr)

} // namespace namespace spec_caise

TEST (MAP_MANAGER_SPEC_CASE, STATIC) {
    using namespace spec_caise;

    test_zero_area (mgr);
    test_first_area (mgr);
    test_last_area (mgr);
}