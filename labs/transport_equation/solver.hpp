#include "lib.hpp"

namespace treq {

struct rect_t {
    int x_i_min, t_i_min;
    int x_chunk_size, t_chunk_size;

    int
    get_x_i_end () const noexcept {
        return x_i_min + x_chunk_size;
    }

    int
    get_t_i_end () const noexcept {
        return t_i_min + t_chunk_size;
    }

    std::size_t
    get_offset_rect (std::size_t buf_x_size) const noexcept {
        return x_i_min + t_i_min * buf_x_size;
    }

    std::size_t
    get_offset_rect (std::size_t buf_x_size,
                     const rect_t& zero_rect) const noexcept {
        return (x_i_min - zero_rect.x_i_min) + (t_i_min - zero_rect.t_i_min) * buf_x_size;
    }
}; // struct rect_t

struct area_params_t {

    rect_t zero_rect;
    int x_i_end, t_i_end;

    int
    get_x_zero_chunk_size () const noexcept {
        return std::min (x_i_end - zero_rect.x_i_min, zero_rect.x_chunk_size);
    }

    int
    get_t_zero_chunk_size () const noexcept {
        return std::min (t_i_end - zero_rect.t_i_min, zero_rect.t_chunk_size);
    }

    int
    get_zero_and_right_rect_x_size () const noexcept {
        return x_i_end - zero_rect.x_i_min;
    }

    int
    get_zero_and_right_rect_size () const noexcept {
        return get_t_zero_chunk_size () * get_zero_and_right_rect_x_size ();
    }

    int
    get_up_rect_size () const noexcept {
        int area_t_size = t_i_end - zero_rect.t_i_min;
        if (area_t_size <= zero_rect.t_chunk_size) {
            return 0;
        } else {
            return (area_t_size - zero_rect.t_chunk_size) * zero_rect.x_chunk_size;
        }
    }

    rect_t
    get_right_rect (int i_chunk) const noexcept {
        rect_t right_rect = {
            .x_i_min = zero_rect.x_i_min + i_chunk * zero_rect.x_chunk_size,
            .t_i_min = zero_rect.t_i_min
        };

        if (right_rect.x_i_min + zero_rect.x_chunk_size > x_i_end) {
            right_rect.x_chunk_size = x_i_end - right_rect.x_i_min;
            right_rect.t_chunk_size = t_i_end - right_rect.t_i_min;
        } else {
            right_rect.x_chunk_size = zero_rect.x_chunk_size;
            right_rect.t_chunk_size = zero_rect.t_chunk_size;
        }

        return right_rect;
    } // rect_t get_right_rect (int i_chunk)

    rect_t
    get_up_rect (int i_chunk) const noexcept {
        rect_t up_rect = {
            .x_i_min = zero_rect.x_i_min,
            .t_i_min = zero_rect.t_i_min + i_chunk * zero_rect.t_chunk_size,
            .x_chunk_size = zero_rect.x_chunk_size
        };

        if (up_rect.t_i_min + zero_rect.t_chunk_size > t_i_end) {
            up_rect.t_chunk_size = t_i_end - up_rect.t_i_min;
        } else {
            up_rect.t_chunk_size = zero_rect.t_chunk_size;
        }

        return up_rect;
    } // rect_t get_up_rect (int i_chunk)
}; // struct area_params_t

struct map_manager_t {
    int x_size;
    int t_size;

    double x_min;
    double x_max;
    double t_min;
    double t_max;

    int num_threads;
    int num_areas;

    map_manager_t (int x_size,
                    int t_size,
                    double x_min,
                    double x_max,
                    double t_min,
                    double t_max,
                    int num_threads,
                    int num_areas) :
        x_size (x_size),
        t_size (t_size),
        x_min (x_min),
        x_max (x_max),
        t_min (t_min),
        t_max (t_max),
        num_threads (num_threads),
        num_areas (num_areas)
    {
        if (x_size <= 1 || t_size <= 1 || num_threads <= 0 ||
            num_areas < num_threads) {
            throw std::invalid_argument ("Incorrect params");
        }
    }

    map_manager_t (const trans_eq_task_t& task,
                    int num_threads,
                    int num_areas) :
        map_manager_t (task.x_size, task.t_size,
                        task.x_min, task.x_max,
                        task.t_min, task.x_max,
                        num_threads, num_areas)
    {}

    area_params_t
    get_area_params (int i_area) const noexcept {
        assert (i_area >= 0);

        area_params_t params = {};

        params.zero_rect.x_chunk_size = calc_chunk_size (x_size, num_areas);
        params.zero_rect.t_chunk_size = calc_chunk_size (t_size, num_areas);

        params.zero_rect.x_i_min = i_area * params.zero_rect.x_chunk_size;
        params.zero_rect.t_i_min = i_area * params.zero_rect.t_chunk_size;

        params.x_i_end = x_size;
        params.t_i_end = t_size;

        return params;
    }

    double
    calc_dx () const noexcept {
        return (x_max - x_min) / (x_size - 1);
    }

    double
    calc_dt () const noexcept {
        return (t_max - t_min) / (t_size - 1);
    }

    std::pair <double, double>
    calc_dx_dt () const noexcept {
        return {calc_dx (), calc_dt ()};
    }
}; // struct map_manager_t

struct trans_eq_solver {
    const map_manager_t map_mgr;
    std::vector <double> u_buf;

    trans_eq_solver (const trans_eq_task_t& task,
                     int num_threads,
                     int koef_mult_threads) :
        map_mgr (task, num_threads, koef_mult_threads * num_threads)
    {}

    void
    solve (int rank) {
        if (rank == 0) {
            auto compute_func = [&] {
                calc_zero_rank_func_border ();
                calc_zero_rank_grid_border ();
            };

            auto bufferize_func = [&] {
                receive_u_bufs ();
            };

            std::thread compute_thread {compute_func};
            std::thread bufferer_thread {bufferize_func};

            compute_thread.join ();
            bufferer_thread.join ();
        } else {
            calc_non_zero_rank_grid_border (rank);
        }
    } // void solve

    void
    calc_zero_rank_func_border () {
        const area_params_t area_params = map_mgr.get_area_params (0);
        const auto[dx, dt] = map_mgr.calc_dx_dt ();
        std::size_t u_buf_x_size = map_mgr.x_size;
        std::vector <double> u_t_buf;

        // 1) Calc u, where need only functional border conditions
        // Fill u(t = 0, x)
        for (int x_i = 0; x_i < area_params.get_x_zero_chunk_size (); ++x_i) {
            double x = x_i * dx;
            u_buf[x_i] = funcs::u_0_x (x);
        }

        // Fill u(t, x = 0)
        for (int t_i = 1; t_i < area_params.get_t_zero_chunk_size (); ++t_i) {
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
    } // void calc_zero_rank_func_border ()

    void
    calc_zero_rank_grid_border () {
        std::size_t u_buf_x_size = map_mgr.x_size;
        std::vector <double> u_t_buf, u_x_buf;

        const int next_rank = 1, prev_rank = map_mgr.num_threads - 1;

        // 2) Calc u, where need gridden border conditions
        // area - is set of chuks: |__
        MPI_Status status = {};
        for (int i_area = map_mgr.num_threads; i_area < map_mgr.num_areas; i_area += map_mgr.num_threads) {
            const area_params_t area_params = map_mgr.get_area_params (i_area);
            const rect_t& zero_rect = area_params.zero_rect;

            u_t_buf.resize (zero_rect.t_chunk_size);
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

            for (int i_chunk_g = 1 + i_area; i_chunk_g < map_mgr.num_areas; ++i_chunk_g) {
                int i_chunk_l = i_chunk_g - i_area;
                const rect_t rect_right = area_params.get_right_rect (i_chunk_g);
                const rect_t rect_up = area_params.get_up_rect (i_chunk_g);
                double* u_buf_right_begin = u_buf.data () + rect_right.get_offset_rect (u_buf_x_size);
                double* u_buf_up_begin = u_buf.data () + rect_up.get_offset_rect (u_buf_x_size);

                u_x_buf[0] = u_x_buf[zero_rect.x_chunk_size];
                u_x_buf.resize (rect_right.x_chunk_size + 1);   // todo delete
                u_t_buf.resize (rect_up.t_chunk_size + 1);      // todo delete

                // Calc u odd
                MPI_Recv (u_x_buf.data () + 1, rect_right.x_chunk_size,
                          MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

                calc_u_rect_x_row (u_buf_right_begin, u_x_buf.data (), u_buf_x_size, rect_right,
                                   funcs::f, funcs::u_next);

                // Send up
                if (i_chunk_l == 1) {
                    MPI_Send (u_buf_right_begin + (rect_up.t_chunk_size - 1) * u_buf_x_size - 1,
                              rect_up.x_chunk_size + 1,
                              MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
                } else {
                    MPI_Send (u_buf_right_begin + (rect_up.t_chunk_size - 1) * u_buf_x_size,
                              rect_up.x_chunk_size,
                              MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
                }

                // Calc u even
                u_t_buf[0] = u_t_buf[zero_rect.t_chunk_size];
                MPI_Recv (u_t_buf.data () + 1, rect_up.t_chunk_size,
                          MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

                calc_u_rect_t_col (u_buf_up_begin, u_t_buf.data (), u_buf_x_size, rect_up,
                                   funcs::f, funcs::u_next);

                // Send right
                copy_col_2_row (u_t_buf.data (),
                                u_buf_up_begin + rect_up.x_chunk_size - 1,
                                u_buf_x_size, rect_up.t_chunk_size);
                MPI_Send (u_t_buf.data (), rect_up.t_chunk_size,
                          MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            }
        }
    } // void calc_zero_rank_grid_border ()

    void
    calc_non_zero_rank_grid_border (int rank) {
        std::size_t u_buf_x_size = map_mgr.x_size;

        const int next_rank = rank == map_mgr.num_threads - 1 ? 0 : rank + 1;
        const int prev_rank = rank - 1;

        std::vector <double> u_buf_right, u_buf_up, u_x_buf, u_t_buf;

        // area - is set of chuks: |__
        MPI_Status status = {};
        for (int i_area = rank; i_area < map_mgr.num_areas; i_area += map_mgr.num_threads) {
            const area_params_t area_params = map_mgr.get_area_params (i_area);
            const rect_t& zero_rect = area_params.zero_rect;

            double* u_buf_begin = u_buf.data () + zero_rect.get_offset_rect (u_buf_x_size);

            // Get x and t array from prev rank
            u_t_buf.resize (zero_rect.t_chunk_size);
            u_x_buf.resize (zero_rect.x_chunk_size + 1);

            u_buf_right.resize (area_params.get_zero_and_right_rect_size ());
            std::size_t u_buf_right_x_size = area_params.get_zero_and_right_rect_x_size ();
            std::size_t u_buf_up_x_size = zero_rect.x_chunk_size;

            MPI_Recv (u_x_buf.data (), zero_rect.x_chunk_size + 1,
                      MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);
            MPI_Recv (u_t_buf.data (), zero_rect.t_chunk_size,
                      MPI_DOUBLE, prev_rank, TAG_BORDER_COND, MPI_COMM_WORLD, &status);

            // Calc zero area
            calc_u_rect_x_row_t_col (u_buf_begin, u_x_buf.data (), u_t_buf.data (),
                                     u_buf_x_size, zero_rect,
                                     funcs::f, funcs::u_next);

            if (area_params.get_up_rect_size () != 0) {
                u_buf_up.resize (area_params.get_up_rect_size () + zero_rect.x_chunk_size);
                std::copy_n (u_buf_right.data () + (zero_rect.t_chunk_size - 1) * u_buf_right_x_size,
                             zero_rect.x_chunk_size, u_buf_up.data ());
            }

            double* u_buf_up_begin_ = u_buf_up.data () + zero_rect.x_chunk_size;
            for (int i_chunk_g = 1 + i_area; i_chunk_g < map_mgr.num_areas; ++i_chunk_g) {
                int i_chunk_l = i_chunk_g - i_area;

                const rect_t rect_right = area_params.get_right_rect (i_chunk_g);
                const rect_t rect_up = area_params.get_up_rect (i_chunk_g);

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

                // Send up
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

                // Send right
                copy_col_2_row (u_t_buf.data (),
                                u_buf_up_begin + u_buf_up_x_size - 1,
                                u_buf_up_x_size, rect_up.t_chunk_size);
                MPI_Send (u_t_buf.data (), rect_up.t_chunk_size,
                          MPI_DOUBLE, next_rank, TAG_BORDER_COND, MPI_COMM_WORLD);
            }

            // Send u_buf_right and u_buf_up to process with rank 0
            MPI_Send (u_buf_right.data (), u_buf_right.size (),
                      MPI_DOUBLE, 0, TAG_SAVE_ON_HOST, MPI_COMM_WORLD);

            MPI_Send (u_buf_up_begin_, area_params.get_up_rect_size (),
                      MPI_DOUBLE, 0, TAG_SAVE_ON_HOST, MPI_COMM_WORLD);
        }
    } // void calc_non_zero_rank_grid_border (int rank)

    void
    receive_u_bufs () {

    } // void receive_u_bufs ()

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_x_row_t_col (T* u_buf,
                             const T* u_x_buf,
                             const T* u_t_buf,
                             std::size_t u_buf_x_size,
                             std::size_t x_i_min,
                             std::size_t t_i_min,
                             std::size_t x_size_calc,
                             std::size_t t_size_calc,
                             F f,
                             U_next u_next)
    {
        // u_k1_m1 = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh)

        auto[dx, dt] = map_mgr.calc_dx_dt ();
        double x_min = map_mgr.x_min + x_i_min * dx;
        double t_min = map_mgr.t_min + t_i_min * dt;

        // Calc left down point
        u_buf[0] = u_next (dx, dt, u_x_buf[0], u_x_buf[1], u_t_buf[0],
                           f (x_min - dx/2, t_min - dt/2));

        // Calc first row of u
        double u_k1_m = u_buf[0];
        for (std::size_t x_i_rel = 1; x_i_rel < x_size_calc; ++x_i_rel) {
            double u_k_m  = u_x_buf[x_i_rel];
            double u_k_m1 = u_x_buf[x_i_rel + 1];

            double x = x_min + x_i_rel * dx;
            double f_kh_mh = f (t_min - dt/2, x - dx/2);

            u_k1_m = u_buf[x_i_rel]
                   = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }

        // Calc first col of u
        double u_k_m1 = u_buf[0];
        for (std::size_t t_i_rel = 1; t_i_rel < t_size_calc; ++t_i_rel) {
            double u_k_m  = u_t_buf[t_i_rel - 1];
            double u_k1_m = u_t_buf[t_i_rel];

            double t = t_min + t_i_rel * dt;
            double f_kh_mh = f (t - dt/2, x_min - dx/2);

            u_k_m1 = u_buf[t_i_rel * u_buf_x_size]
                   = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }

        // Calc rect (x_size_calc - 1) x (t_size_calc - 1)
        calc_u_rect (u_buf + u_buf_x_size + 1, u_buf_x_size,
                     x_i_min + 1, t_i_min + 1, x_size_calc -1, t_size_calc -1,
                     funcs::f, funcs::u_next);
    } // void calc_u_rect_x_row_t_col (T* u_buf, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_x_row_t_col (T* u_buf,
                             const T* u_x_buf,
                             const T* u_t_buf,
                             std::size_t u_buf_x_size,
                             const rect_t& rect_calc,
                             F f,
                             U_next u_next)
    {
        calc_u_rect_x_row_t_col (u_buf, u_x_buf, u_t_buf, u_buf_x_size,
                                 rect_calc.x_i_min, rect_calc.t_i_min,
                                 rect_calc.x_chunk_size, rect_calc.t_chunk_size,
                                 f, u_next);
    } // void calc_u_rect_x_row_t_col (T* u_buf, ..., const rect_t& rect_calc, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_x_row (T* u_buf,
                       const T* u_x_buf,
                       std::size_t u_buf_x_size,
                       std::size_t x_i_min,
                       std::size_t t_i_min,
                       std::size_t x_size_calc,
                       std::size_t t_size_calc,
                       F f,
                       U_next u_next)
    {
        // u_k1_m1 = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh)

        auto[dx, dt] = map_mgr.calc_dx_dt ();
        double x_min = map_mgr.x_min + x_i_min * dx;
        double t_min = map_mgr.t_min + t_i_min * dt;

        // Calc first row of u
        for (std::size_t x_i_rel = 0; x_i_rel < x_size_calc; ++x_i_rel) {
            double u_k_m  = u_x_buf[x_i_rel];
            double u_k_m1 = u_x_buf[x_i_rel + 1];
            double u_k1_m = u_buf[x_i_rel - 1];

            double x = x_min + x_i_rel * dx;
            double f_kh_mh = f (t_min - dt/2, x - dx/2);

            u_buf[x_i_rel] = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }

        // Calc (rect x_size_calc) x (t_size_calc - 1) // todo
        calc_u_rect (u_buf + u_buf_x_size, u_buf_x_size, x_i_min, t_i_min + 1,
                     x_size_calc, t_size_calc - 1,
                     funcs::f, funcs::u_next);
    } // void calc_u_rect_x_row (T* u_buf, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_x_row (T* u_buf,
                       const T* u_x_buf,
                       std::size_t u_buf_x_size,
                       const rect_t& rect_calc,
                       F f,
                       U_next u_next)
    {
        calc_u_rect_x_row (u_buf, u_x_buf, u_buf_x_size,
                           rect_calc.x_i_min, rect_calc.t_i_min,
                           rect_calc.x_chunk_size, rect_calc.t_chunk_size,
                           f, u_next);
    } // void calc_u_rect_x_row (T* u_buf, ..., const rect_t& rect_calc, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_t_col (T* u_buf,
                       const T* u_t_buf,
                       std::size_t u_buf_x_size,
                       std::size_t x_i_min,
                       std::size_t t_i_min,
                       std::size_t x_size_calc,
                       std::size_t t_size_calc,
                       F f,
                       U_next u_next)
    {
        // u_k1_m1 = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh)

        auto[dx, dt] = map_mgr.calc_dx_dt ();
        double x_min = map_mgr.x_min + x_i_min * dx;
        double t_min = map_mgr.t_min + t_i_min * dt;

        // Calc first col of u
        for (std::size_t t_i_rel = 0; t_i_rel < t_size_calc; ++t_i_rel) {
            double* u_buf_cur = u_buf + t_i_rel * u_buf_x_size;
            double* u_buf_prev = u_buf_cur - u_buf_x_size;

            double u_k_m  = u_t_buf[t_i_rel];
            double u_k1_m = u_t_buf[t_i_rel + 1];
            double u_k_m1 = *u_buf_prev;

            double t = t_min + t_i_rel * dt;
            double f_kh_mh = f (t - dt/2, x_min - dx/2);

            *u_buf_cur = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
        }

        // Calc rect of u
        calc_u_rect (u_buf + 1, u_buf_x_size, x_i_min + 1, t_i_min,
                     x_size_calc - 1, t_size_calc,
                     funcs::f, funcs::u_next);
    } // void calc_u_rect_t_col (T* u_buf, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect_t_col (T* u_buf,
                       const T* u_t_buf,
                       std::size_t u_buf_x_size,
                       const rect_t& rect_calc,
                       F f,
                       U_next u_next)
    {
        calc_u_rect_t_col (u_buf, u_t_buf, u_buf_x_size,
                           rect_calc.x_i_min, rect_calc.t_i_min,
                           rect_calc.x_chunk_size, rect_calc.t_chunk_size,
                           f, u_next);
    } // void calc_u_rect_t_col (T* u_buf, ..., const rect_t& rect_calc, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect (T* u_buf,
                 std::size_t u_buf_x_size,
                 std::size_t x_i_min,
                 std::size_t t_i_min,
                 std::size_t x_size_calc,
                 std::size_t t_size_calc,
                 F f,
                 U_next u_next)
    {
        // u_k1_m1 = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh)

        auto[dx, dt] = map_mgr.calc_dx_dt ();
        double x_min = map_mgr.x_min + x_i_min * dx;
        double t_min = map_mgr.t_min + t_i_min * dt;

        // Calc rect (x_size_calc - 1) x (t_size_calc - 1)
        for (std::size_t t_i_rel = 0; t_i_rel < t_size_calc; ++t_i_rel) {
            double* u_row_cur = u_buf + t_i_rel * u_buf_x_size;
            double* u_row_prev = u_buf - u_buf_x_size;

            double t = t_min + t_i_rel * dt;
            for (std::size_t x_i_rel = 0; x_i_rel < x_size_calc; ++x_i_rel) {
                double u_k_m  = u_row_prev[-1];
                double u_k_m1 = u_row_prev[0];
                double u_k1_m = u_row_cur[-1];

                double x = x_min + x_i_rel * dx;
                double f_kh_mh = f (t - dt/2, x - dx/2);

                u_row_cur[x_i_rel] = u_next (dx, dt, u_k_m, u_k_m1, u_k1_m, f_kh_mh);
            }
        }
    } // void calc_u_rect (T* u_buf, ..., U_next u_next)

    template <typename T, typename F, typename U_next>
    void
    calc_u_rect (T* u_buf,
                 std::size_t u_buf_x_size,
                 const rect_t& rect_calc,
                 F f,
                 U_next u_next)
    {
        calc_u_rect (u_buf, u_buf_x_size,
                     rect_calc.x_i_min, rect_calc.t_i_min,
                     rect_calc.x_chunk_size, rect_calc.t_chunk_size,
                     f, u_next);
    } // void calc_u_rect (T* u_buf, ..., const rect_t& rect_calc, ..., U_next u_next)

}; // struct trans_eq_solver

} // namespace treq