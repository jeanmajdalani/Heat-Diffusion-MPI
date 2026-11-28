#pragma once

#ifdef HEAT_HAVE_MPI

#include <mpi.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "heat_common.hpp"

namespace heat_blocks {

inline int run(const heat_common::Options& opt,
               const heat_common::TimeSetup& ts) {
    int rank = 0, P = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    const int N = opt.N;
    const int M = N - 1;
    const int side = N + 1;

    // 1. Create a non-periodic Cartesian process grid.
    // Dimension 0 corresponds to y; dimension 1 to x.
    int dims[2] = {0, 0};
    MPI_Dims_create(P, 2, dims);

    if (dims[0] > M || dims[1] > M) {
        if (rank == 0)
            std::cerr << "Error: too many blocks for this grid\n";
        return 1;
    }

    int periods[2] = {0, 0};
    MPI_Comm cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims,
                    periods, 0, &cart);

    MPI_Comm_rank(cart, &rank);

    int coords[2];
    MPI_Cart_coords(cart, rank, 2, coords);

    // 2. Distribute interior nodes evenly.
    auto count_of = [&](int dim, int c) {
        return M / dims[dim] + (c < M % dims[dim] ? 1 : 0);
    };

    auto first_of = [&](int dim, int c) {
        return 1 + c * (M / dims[dim])
                 + std::min(c, M % dims[dim]);
    };

    const int ny = count_of(0, coords[0]);
    const int nx = count_of(1, coords[1]);
    const int y0 = first_of(0, coords[0]);
    const int x0 = first_of(1, coords[1]);
    const int pitch = nx + 2;

    // 3. Identify the four neighboring processes.
    int down, up, left, right;
    MPI_Cart_shift(cart, 0, 1, &down, &up);
    MPI_Cart_shift(cart, 1, 1, &left, &right);

    auto idx = [&](int i, int j) {
        return heat_common::idx(i, j, pitch);
    };

    // 4. Allocate owned nodes and ghost layers.
    const std::size_t local_size =
        static_cast<std::size_t>(ny + 2) * pitch;

    std::vector<double> old_grid(local_size);
    std::vector<double> new_grid(local_size);

    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i) {
            const double v = heat_common::initial_value(
                x0 + i - 1, y0 + j - 1, N, opt.test_case);

            old_grid[idx(i, j)] = v;
            new_grid[idx(i, j)] = v;
        }
    }

    // 5. Derived MPI datatype for non-contiguous columns.
    MPI_Datatype column;
    MPI_Type_vector(ny, 1, pitch, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);

    const double r = ts.effective_r;
    double t_comm = 0.0;
    double t_comp = 0.0;

    MPI_Barrier(cart);
    const double t_start = MPI_Wtime();

    // 6. Main time-stepping loop.
    for (long long step = 0; step < ts.steps; ++step) {
        const double a = MPI_Wtime();

        // Exchange horizontal rows.
        MPI_Sendrecv(
            &old_grid[idx(1, ny)], nx, MPI_DOUBLE, up, 0,
            &old_grid[idx(1, 0)], nx, MPI_DOUBLE, down, 0,
            cart, MPI_STATUS_IGNORE);

        MPI_Sendrecv(
            &old_grid[idx(1, 1)], nx, MPI_DOUBLE, down, 1,
            &old_grid[idx(1, ny + 1)], nx, MPI_DOUBLE, up, 1,
            cart, MPI_STATUS_IGNORE);

        // Exchange vertical columns without manual packing.
        MPI_Sendrecv(
            &old_grid[idx(nx, 1)], 1, column, right, 2,
            &old_grid[idx(0, 1)], 1, column, left, 2,
            cart, MPI_STATUS_IGNORE);

        MPI_Sendrecv(
            &old_grid[idx(1, 1)], 1, column, left, 3,
            &old_grid[idx(nx + 1, 1)], 1, column, right, 3,
            cart, MPI_STATUS_IGNORE);

        const double b = MPI_Wtime();

        // Same arithmetic order as the serial solver.
        for (int j = 1; j <= ny; ++j) {
            for (int i = 1; i <= nx; ++i) {
                const std::size_t p = idx(i, j);

                const double neighbors =
                    old_grid[p - 1] + old_grid[p + 1]
                  + old_grid[p - pitch] + old_grid[p + pitch];

                new_grid[p] =
                    (1.0 - 4.0 * r) * old_grid[p]
                    + r * neighbors;
            }
        }

        old_grid.swap(new_grid);

        const double c = MPI_Wtime();
        t_comm += b - a;
        t_comp += c - b;
    }

    MPI_Barrier(cart);
    const double t_total_local = MPI_Wtime() - t_start;

    // 7. Collect maximum execution times across processes.
    double t_total = 0.0;
    double t_comp_max = 0.0;
    double t_comm_max = 0.0;

    MPI_Reduce(&t_total_local, &t_total,
               1, MPI_DOUBLE, MPI_MAX, 0, cart);
    MPI_Reduce(&t_comp, &t_comp_max,
               1, MPI_DOUBLE, MPI_MAX, 0, cart);
    MPI_Reduce(&t_comm, &t_comm_max,
               1, MPI_DOUBLE, MPI_MAX, 0, cart);

    MPI_Type_free(&column);

    // 8. Reconstruct the global field when required.
    const bool need_global =
        !opt.output_path.empty() || opt.test_case == "V";

    std::vector<double> global;

    if (need_global) {
        if (1LL * M * M > std::numeric_limits<int>::max()) {
            if (rank == 0)
                std::cerr << "Error: grid too large for MPI_Gatherv\n";
            MPI_Comm_free(&cart);
            return 1;
        }

        // Pack owned interior nodes only, excluding ghosts.
        std::vector<double> packed(
            static_cast<std::size_t>(nx) * ny);

        for (int j = 1; j <= ny; ++j) {
            std::copy_n(
                &old_grid[idx(1, j)], nx,
                packed.data() + static_cast<std::size_t>(j - 1) * nx);
        }

        std::vector<int> counts, displs;
        std::vector<double> gathered;

        if (rank == 0) {
            counts.resize(P);
            displs.resize(P);

            int offset = 0;

            for (int q = 0; q < P; ++q) {
                int qc[2];
                MPI_Cart_coords(cart, q, 2, qc);

                counts[q] =
                    count_of(0, qc[0]) * count_of(1, qc[1]);

                displs[q] = offset;
                offset += counts[q];
            }

            gathered.resize(offset);
        }

        MPI_Gatherv(
            packed.data(), nx * ny, MPI_DOUBLE,
            rank == 0 ? gathered.data() : nullptr,
            counts.data(), displs.data(),
            MPI_DOUBLE, 0, cart);

        if (rank == 0) {
            global.resize(static_cast<std::size_t>(side) * side);
            heat_common::initialize(global, N, opt.test_case);

            for (int q = 0; q < P; ++q) {
                int qc[2];
                MPI_Cart_coords(cart, q, 2, qc);

                const int qny = count_of(0, qc[0]);
                const int qnx = count_of(1, qc[1]);
                const int qy0 = first_of(0, qc[0]);
                const int qx0 = first_of(1, qc[1]);

                for (int j = 0; j < qny; ++j) {
                    for (int i = 0; i < qnx; ++i) {
                        global[heat_common::idx(
                            qx0 + i, qy0 + j, side)] =
                            gathered[displs[q] + j * qnx + i];
                    }
                }
            }
        }
    }

    // 9. Output results using the shared project format.
    if (rank == 0) {
        try {
            if (!opt.output_path.empty())
                heat_common::write_binary(opt.output_path, global);

            std::cout << std::setprecision(17)
                      << "decomp=blocks P=" << P
                      << " N=" << N
                      << " steps=" << ts.steps
                      << " r=" << ts.effective_r
                      << " dt=" << ts.dt
                      << " T=" << ts.final_time
                      << " t_total=" << t_total
                      << " t_comp=" << t_comp_max
                      << " t_comm=" << t_comm_max;

            if (opt.test_case == "V") {
                double err_inf = 0.0;
                double err_l2_h = 0.0;

                heat_common::error_norms(
                    global, N, ts.final_time, err_inf, err_l2_h);

                std::cout << " err_inf=" << err_inf
                          << " err_l2_h=" << err_l2_h;
            }

            std::cout << '\n';

        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << '\n';
            MPI_Abort(cart, 1);
        }
    }

    MPI_Comm_free(&cart);
    return 0;
}

} // namespace heat_blocks

#endif // HEAT_HAVE_MPI
