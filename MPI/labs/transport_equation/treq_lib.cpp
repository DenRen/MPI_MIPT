#include "treq_lib.hpp"

#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <thread>
#include <mpi/mpi.h>

namespace treq_seq {

void
check_error (int err_code, int line, const char* file_name) {
    if (err_code != MPI_SUCCESS) {
        char err_str[MPI_MAX_ERROR_STRING] = "";
        int str_len = 0;
        MPI_Error_string (err_code, err_str, &str_len);

        std::string err_msg {file_name};
        err_msg += ": " + std::to_string (line);
        err_msg += ": ";
        err_msg += err_str;

        throw std::runtime_error (err_msg);
    }
}

std::vector <double>
solve_transport_eq (int num_points_coord_x,
                    int num_points_time,
                    double x_max,
                    double t_max)
{
    const double dx = x_max / (num_points_coord_x - 1);
    const double dt = t_max / (num_points_time - 1);

    using funcs::u_t_0;
    using funcs::u_0_x;
    using funcs::f;
    using funcs::u_next;
    return calc_u (num_points_coord_x, num_points_time, dx, dt,
                   u_t_0, u_0_x, f, u_next);
}

std::vector <double>
get_solve_ref (int num_points_coord_x,
               int num_points_time,
               double x_max,
               double t_max)
{
    const double dx = x_max / (num_points_coord_x - 1);
    const double dt = t_max / (num_points_time - 1);

    using funcs::u_solve;
    return calc_u (num_points_coord_x, num_points_time, dx, dt,
                   u_solve);
}

//                        h      err
std::vector <std::pair <double, double>>
calc_h_err_vec (int num_points_coord_x_min,
                int num_points_coord_dx,
                int num_steps,
                int num_points_time,
                double x_max,
                double t_max)
{
    std::vector <std::pair <double, double>> h_err (num_steps);
    for (int step = 0; step < num_steps; ++step) {
        int num_points_coord_x =
            num_points_coord_x_min + num_points_coord_dx * step;

        auto u = solve_transport_eq (num_points_coord_x,
                                     num_points_time, x_max, t_max);

        double h = x_max / (num_points_coord_x - 1);
        double err = calc_err_use_solve_func (u, num_points_coord_x,
                                              num_points_time, x_max, t_max,
                                              funcs::u_solve);

        h_err[step] = {h, err};
    }

    return h_err;
}

void
check_single_thread_solution () {
    unsigned M = 800;       // coord
    unsigned K = 1'000'000; // time

    double X_max = 6;       // X_min = 0
    double T_max = 6;       // T_min = 0

// #define ONLY_SOLVE

#ifdef ONLY_SOLVE
    auto u = solve_transport_eq (M, K, X_max, T_max);
    { volatile double x = u[1]; }
#else
    auto h_err = calc_h_err_vec (M, M / 10, 3, K, X_max, T_max);
    for (auto[h, err] : h_err) {
        std::cout << std::log (h) << ' ' << std::log (err) << '\n';
    }

    for (auto&[h, err] : h_err) {
        h = std::log (h);
        err = std::log (err);
    }
    std::cout << "O(h^" << linearize (h_err).first << ")\n";
#endif

#undef ONLY_SOLVE
} // void check_single_thread_solution ()

} // namespace treq_seq

namespace treq {

area_params_t
map_manager_t::get_area_params (int i_area) const noexcept {
    assert (i_area >= 0);

    area_params_t params = {};
    rect_t& zero_rect = params.zero_rect;

    int x_chunk_size = calc_chunk_size (x_size, num_areas);
    int t_chunk_size = calc_chunk_size (t_size, num_areas);

    zero_rect.x_i_min = i_area * x_chunk_size;
    zero_rect.t_i_min = i_area * t_chunk_size;

    if (x_chunk_size + zero_rect.x_i_min <= x_size) {
        zero_rect.x_chunk_size = x_chunk_size;
        zero_rect.t_chunk_size = t_chunk_size;
    } else {
        zero_rect.x_chunk_size = x_size - zero_rect.x_i_min;
        zero_rect.t_chunk_size = t_size - zero_rect.t_i_min;
    }

    params.x_i_end = x_size;
    params.t_i_end = t_size;

    return params;
} // area_params_t map_manager_t::get_area_params (int i_area)

void
trans_eq_solver::solve (int rank) {
    if (rank == 0) {
        std::vector <double> u_buf (map_mgr.x_size * map_mgr.t_size);

        calc_zero_rank_func_border (u_buf);
        receive_u_bufs (u_buf);

        // if (map_mgr.num_threads == 1) {
        //     calc_zero_rank_func_border (u_buf);
        //     receive_u_bufs (u_buf);
        // } else {
        //     auto compute_func = [&] {
        //         calc_zero_rank_func_border (u_buf);
        //         calc_zero_rank_grid_border (u_buf);
        //     };

        //     auto bufferize_func = [&] {
        //         receive_u_bufs (u_buf);
        //     };

        //     std::thread compute_thread {compute_func};
        //     std::thread bufferer_thread {bufferize_func};

        //     compute_thread.join ();
        //     bufferer_thread.join ();
        // }

        { volatile double tr = u_buf[3]; }

        print_2d_array (u_buf.data (), map_mgr.x_size, u_buf.size ());
    } else {
        calc_non_zero_rank_grid_border (rank);
    }
} // void trans_eq_solver::solve (int rank)

void
trans_eq_solver::calc_zero_rank_func_border (std::vector <double>& u_buf) {
    const area_params_t area_params = map_mgr.get_area_params (0);
    const auto[dx, dt] = map_mgr.calc_dx_dt ();
    std::size_t u_buf_x_size = map_mgr.x_size;
    std::vector <double> u_t_buf;

    // 1) Calc u, where need only functional border conditions
    // Fill u(t = 0, x)
    for (int x_i = 0; x_i < area_params.get_zero_chunk_x_size (); ++x_i) {
        double x = x_i * dx;
        u_buf[x_i] = funcs::u_0_x (x);
    }

    // Fill u(t, x = 0)
    for (int t_i = 1; t_i < area_params.get_zero_chunk_t_size (); ++t_i) {
        double t = t_i * dt;
        u_buf[t_i * u_buf_x_size] = funcs::u_t_0 (t);
    }

    // Calc u on (x_size - 1) * (t_size - 1)
    const rect_t& zero_rect = area_params.zero_rect;
    calc_u_rect (u_buf.data () + u_buf_x_size + 1, u_buf_x_size,
                    zero_rect.x_i_min + 1, zero_rect.t_i_min + 1,
                    zero_rect.x_chunk_size - 1, zero_rect.t_chunk_size - 1,
                    funcs::f, funcs::u_next);

    // Calc u rights and ups in zero area
    const int next_rank = 1;
    for (int i_chunk = 1; i_chunk < map_mgr.num_areas; ++i_chunk) {
        // Fill u(t = 0, x)
        const rect_t rect_right = area_params.get_right_rect (i_chunk);
        for (int x_i = rect_right.x_i_min; x_i < rect_right.get_x_i_end (); ++x_i) {
            double x = x_i * dx;
            u_buf[x_i] = funcs::u_0_x (x);
        }

        // Calc right rect u
        double* u_buf_right_begin = u_buf.data () + rect_right.x_i_min;
        calc_u_rect (u_buf_right_begin + u_buf_x_size, u_buf_x_size,
                        rect_right.x_i_min, rect_right.t_i_min + 1,
                        rect_right.x_chunk_size, rect_right.t_chunk_size - 1,
                        funcs::f, funcs::u_next);

        // Send up
        if (i_chunk == 1) { // right_rect.x_chunk_size + 1
            MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_x_size - 1,
                        rect_right.x_chunk_size + 1,
                        MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
        } else {
            MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_x_size,
                        rect_right.x_chunk_size,
                        MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
        }

        // Fill u(t = 0, x)
        const rect_t rect_up = area_params.get_up_rect (i_chunk);
        for (int t_i = rect_up.t_i_min; t_i < rect_up.get_t_i_end (); ++t_i) {
            double t = t_i * dt;
            u_buf[t_i * u_buf_x_size] = funcs::u_t_0 (t);
        }

        // Calc up rect u
        double* u_buf_up_begin = u_buf.data () + rect_up.t_i_min * u_buf_x_size;
        calc_u_rect (u_buf_up_begin + 1, u_buf_x_size,
                        rect_up.x_i_min + 1, rect_up.t_i_min,
                        rect_up.x_chunk_size - 1, rect_up.t_chunk_size,
                        funcs::f, funcs::u_next);

        // Copy to u_t_buf
        u_t_buf.resize (rect_up.t_chunk_size);
        copy_col_2_row (u_t_buf.data (),
                        u_buf_up_begin + rect_up.x_chunk_size - 1,
                        u_buf_x_size, rect_up.t_chunk_size);
        // print_2d_array (u_buf, task.x_size);
        // Send right
        MPI_Send (u_t_buf.data (), rect_up.t_chunk_size,
                    MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
    }
} // void trans_eq_solver::calc_zero_rank_func_border ()

void
trans_eq_solver::calc_zero_rank_grid_border (std::vector <double>& u_buf) {
    std::size_t u_buf_x_size = map_mgr.x_size;
    std::vector <double> u_t_buf, u_x_buf;

    const int next_rank = 1, prev_rank = map_mgr.num_threads - 1;

    // 2) Calc u, where need gridden border conditions
    // area - is set of chuks: |__
    MPI_Status status = {};
    for (int i_area = map_mgr.num_threads; i_area < map_mgr.num_areas; i_area += map_mgr.num_threads) {
        const area_params_t area_params = map_mgr.get_area_params (i_area);
        const rect_t& zero_rect = area_params.zero_rect;

        u_t_buf.resize (zero_rect.t_chunk_size + 1);
        u_x_buf.resize (zero_rect.x_chunk_size + 1);

        std::size_t zero_pos = zero_rect.t_i_min * u_buf_x_size + zero_rect.x_i_min;
        double* u_buf_begin = u_buf.data () + zero_pos;

        // Get x and t array from prev rank
        MPI_Recv (u_x_buf.data (), zero_rect.x_chunk_size + 1,
                    MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
        MPI_Recv (u_t_buf.data (), zero_rect.t_chunk_size,
                    MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

        // Calc zero area
        calc_u_rect_x_row_t_col (u_buf_begin, u_x_buf.data (), u_t_buf.data (),
                                    u_buf_x_size, zero_rect,
                                    funcs::f, funcs::u_next);

        u_t_buf[zero_rect.t_chunk_size] = u_t_buf[zero_rect.t_chunk_size - 1];
        for (int i_chunk_g = 1 + i_area; i_chunk_g < map_mgr.num_areas; ++i_chunk_g) {
            int i_chunk_l = i_chunk_g - i_area;

            const rect_t rect_right = area_params.get_right_rect (i_chunk_l);
            const rect_t rect_up = area_params.get_up_rect (i_chunk_l);

            double* u_buf_right_begin = u_buf.data () + rect_right.get_offset_rect (u_buf_x_size);
            double* u_buf_up_begin = u_buf.data () + rect_up.get_offset_rect (u_buf_x_size);

            // Calc u odd
            u_x_buf[0] = u_x_buf[zero_rect.x_chunk_size];
            MPI_Recv (u_x_buf.data () + 1, rect_right.x_chunk_size,
                        MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

            calc_u_rect_x_row (u_buf_right_begin, u_x_buf.data (), u_buf_x_size, rect_right,
                                funcs::f, funcs::u_next);

            // Send up row
            if (i_chunk_l == 1) {
                MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_x_size - 1,
                            rect_right.x_chunk_size + 1,
                            MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            } else {
                MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_x_size,
                            rect_right.x_chunk_size,
                            MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            }

            // Calc u even
            u_t_buf[0] = u_t_buf[zero_rect.t_chunk_size];
            MPI_Recv (u_t_buf.data () + 1, rect_up.t_chunk_size,
                        MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

            calc_u_rect_t_col (u_buf_up_begin, u_t_buf.data (), u_buf_x_size, rect_up,
                                funcs::f, funcs::u_next);

            // Send right col
            copy_col_2_row (u_t_buf.data (),
                            u_buf_up_begin + rect_up.x_chunk_size - 1,
                            u_buf_x_size, rect_up.t_chunk_size);
            MPI_Send (u_t_buf.data (), rect_up.t_chunk_size,
                        MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
        }
    }
} // void trans_eq_solver::calc_zero_rank_grid_border ()

void
trans_eq_solver::calc_non_zero_rank_grid_border (int rank) {
    std::size_t u_buf_x_size = map_mgr.x_size;

    const int next_rank = rank == map_mgr.num_threads - 1 ? 0 : rank + 1;
    const int prev_rank = rank - 1;

    std::vector <double> u_buf_right, u_buf_up, u_x_buf, u_t_buf;

    // area - is set of chuks: |__
    MPI_Status status = {};
    for (int i_area = rank; i_area < map_mgr.num_areas; i_area += map_mgr.num_threads) {
        const area_params_t area_params = map_mgr.get_area_params (i_area);
        const rect_t& zero_rect = area_params.zero_rect;

        // Get x and t array from prev rank
        u_t_buf.resize (zero_rect.t_chunk_size + 1);
        u_x_buf.resize (zero_rect.x_chunk_size + 1);

        u_buf_right.resize (area_params.get_zero_and_right_rect_size ());
        std::size_t u_buf_right_x_size = area_params.get_zero_and_right_rect_x_size ();
        std::size_t u_buf_up_x_size = area_params.get_up_rect_x_size ();

        MPI_Recv (u_x_buf.data (), zero_rect.x_chunk_size + 1,
                    MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
        MPI_Recv (u_t_buf.data (), zero_rect.t_chunk_size,
                    MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

        // Calc zero area
        calc_u_rect_x_row_t_col (u_buf_right.data (), u_x_buf.data (), u_t_buf.data (),
                                    u_buf_right_x_size, zero_rect,
                                    funcs::f, funcs::u_next);

        if (area_params.get_up_rect_size () != 0) {
            u_buf_up.resize (area_params.get_up_rect_size () + u_buf_up_x_size);
            std::copy_n (u_buf_right.data () + (zero_rect.t_chunk_size - 1) * u_buf_right_x_size,
                            u_buf_up_x_size, u_buf_up.data ());
            u_t_buf[zero_rect.t_chunk_size] = u_t_buf[zero_rect.t_chunk_size - 1];
        }

        double* u_buf_up_begin_ = u_buf_up.data () + u_buf_up_x_size;
        for (int i_chunk_g = 1 + i_area; i_chunk_g < map_mgr.num_areas; ++i_chunk_g) {
            int i_chunk_l = i_chunk_g - i_area;

            const rect_t rect_right = area_params.get_right_rect (i_chunk_l);
            const rect_t rect_up = area_params.get_up_rect (i_chunk_l);

            double* u_buf_right_begin = u_buf_right.data () +
                                        rect_right.get_offset_rect (u_buf_right_x_size, zero_rect);
            double* u_buf_up_begin = u_buf_up_begin_ +
                                        (i_chunk_l - 1) * u_buf_up_x_size * zero_rect.t_chunk_size;

            // Calc u odd
            u_x_buf[0] = u_x_buf[zero_rect.x_chunk_size];
            MPI_Recv (u_x_buf.data () + 1, rect_right.x_chunk_size,
                        MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

            calc_u_rect_x_row (u_buf_right_begin, u_x_buf.data (),
                                u_buf_right_x_size, rect_right,
                                funcs::f, funcs::u_next);

            // Send up row
            if (i_chunk_l == 1) {
                MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_right_x_size - 1,
                            rect_right.x_chunk_size + 1,
                            MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            } else {
                MPI_Send (u_buf_right_begin + (rect_right.t_chunk_size - 1) * u_buf_right_x_size,
                            rect_right.x_chunk_size,
                            MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            }

            // Calc u even
            u_t_buf[0] = u_t_buf[zero_rect.t_chunk_size];
            MPI_Recv (u_t_buf.data () + 1, rect_up.t_chunk_size,
                        MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

            calc_u_rect_t_col (u_buf_up_begin, u_t_buf.data (),
                                u_buf_up_x_size, rect_up,
                                funcs::f, funcs::u_next);

            // Send right col
            copy_col_2_row (u_t_buf.data (),
                            u_buf_up_begin + u_buf_up_x_size - 1,
                            u_buf_up_x_size, rect_up.t_chunk_size);

            MPI_Send (u_t_buf.data (), rect_up.t_chunk_size,
                        MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
        }

        // Send u_buf_right and u_buf_up to process with rank 0
        MPI_Send (u_buf_right.data (), area_params.get_zero_and_right_rect_size (),
                    MPI_DOUBLE, 0, TAG_SAVE_ON_HOST, MPI_COMM_WORLD);

        if (area_params.get_up_rect_size () != 0) {
            MPI_Send (u_buf_up_begin_, area_params.get_up_rect_size (),
                        MPI_DOUBLE, 0, TAG_SAVE_ON_HOST, MPI_COMM_WORLD);
        }
    }
} // void trans_eq_solver::calc_non_zero_rank_grid_border (int rank)

void
trans_eq_solver::receive_u_bufs (std::vector <double>& u_buf) {
    const std::size_t u_buf_x_size = map_mgr.x_size;
    std::vector <double> u_buf_right, u_buf_up;

    MPI_Status status = {};
    for (int i_area = 1; i_area < map_mgr.num_areas; ++i_area) {
        int source = i_area % map_mgr.num_threads;
        if (source == 0) {
            continue;
        }
        // DUMP (source);

        const area_params_t area_params = map_mgr.get_area_params (i_area);
        const rect_t& zero_rect = area_params.zero_rect;

        const std::size_t u_buf_right_x_size = area_params.get_zero_and_right_rect_x_size ();
        const std::size_t u_buf_right_size = area_params.get_zero_and_right_rect_size ();
        const std::size_t u_buf_up_x_size = area_params.get_up_rect_x_size ();
        const std::size_t u_buf_up_size = area_params.get_up_rect_size ();

        u_buf_right.resize (u_buf_right_size);  // todo delete
        u_buf_up.resize (u_buf_up_size);        // todo delete

        // Receive right buffer
        MPI_Recv (u_buf_right.data (), u_buf_right_size,
                    MPI_DOUBLE, source, TAG_SAVE_ON_HOST, MPI_COMM_WORLD, &status);

        double* u_buf_right_begin = u_buf.data () + zero_rect.get_offset_rect (u_buf_x_size);
        copy_row_2_rect (u_buf_right_begin, u_buf_right.data (),
                            u_buf_right_x_size, u_buf_x_size, u_buf_right_size);

        if (u_buf_up_size != 0) {
            // Receive up buffer
            MPI_Recv (u_buf_up.data (), u_buf_up_size,
                        MPI_DOUBLE, source, TAG_SAVE_ON_HOST, MPI_COMM_WORLD, &status);

            double* u_buf_up_begin = u_buf_right_begin + zero_rect.t_chunk_size * u_buf_x_size;
            copy_row_2_rect (u_buf_up_begin, u_buf_up.data (),
                                u_buf_up_x_size, u_buf_x_size, u_buf_up_size);
        }
    }
} // void trans_eq_solver::receive_u_bufs ()

} // namespace treq

void
solve_trans_eq_parallel (const treq::trans_eq_task_t& task,
                         int* argc_ptr,
                         char** argv_ptr[])
{
    MPI_Init (argc_ptr, argv_ptr);

    int num_threads = 0, rank = 0;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &num_threads);

    int k_zone = 1;

    treq::trans_eq_solver solver {task, num_threads, k_zone};
    solver.solve (rank);

    MPI_Finalize ();
} // void solve_trans_eq_parallel
