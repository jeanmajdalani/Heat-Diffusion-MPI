#pragma once
// MPI 1D "strips" decomposition: the N-1 interior rows are split among P
// processes. Each process stores its rows plus one ghost row above and below.
#ifdef HEAT_HAVE_MPI

#include <mpi.h>

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "heat_common.hpp"

namespace heat_strips {

// Returns the process exit code. Must be called between MPI_Init and MPI_Finalize.
inline int run(const heat_common::Options& opt, const heat_common::TimeSetup& ts) {
    int rank = 0, P = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    const int N = opt.N;
    const int side = N + 1;
    const int M = N - 1;  // interior rows (global rows 1..N-1)

    // All ranks know P and M, so all take this branch together.
    if (P > M) {
        if (rank == 0) {
            std::cerr << "Error: strips need P <= N-1 (P=" << P << ", N=" << N << ")\n";
        }
        return 1;
    }

    // Row distribution: the first (M % P) ranks get one extra row.
    const int base = M / P;
    const int rem = M % P;
    auto rows_of = [&](int r) { return base + (r < rem ? 1 : 0); };
    auto first_row_of = [&](int r) { return 1 + r * base + (r < rem ? r : rem); };
    const int nloc = rows_of(rank);
    const int j0 = first_row_of(rank);  // global index of local row 1

    const int down = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
    const int up = (rank < P - 1) ? rank + 1 : MPI_PROC_NULL;

    // Local rows 0..nloc+1 hold global rows j0-1 .. j0+nloc.
    const std::size_t local_size = static_cast<std::size_t>(nloc + 2) * side;
    std::vector<double> old_grid(local_size), new_grid(local_size);
    for (int l = 0; l <= nloc + 1; ++l) {
        const int j = j0 + l - 1;
        for (int i = 0; i <= N; ++i) {
            const double v = heat_common::initial_value(i, j, N, opt.test_case);
            old_grid[heat_common::idx(i, l, side)] = v;
            new_grid[heat_common::idx(i, l, side)] = v;
        }
    }

    const double r = ts.effective_r;
    double t_comm = 0.0, t_comp = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    const double t_start = MPI_Wtime();
    for (long long step = 0; step < ts.steps; ++step) {
        const double a = MPI_Wtime();
        // Send my top owned row up, receive my bottom ghost row from below.
        MPI_Sendrecv(&old_grid[heat_common::idx(0, nloc, side)], side, MPI_DOUBLE, up, 0,
                     &old_grid[heat_common::idx(0, 0, side)], side, MPI_DOUBLE, down, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Send my bottom owned row down, receive my top ghost row from above.
        MPI_Sendrecv(&old_grid[heat_common::idx(0, 1, side)], side, MPI_DOUBLE, down, 1,
                     &old_grid[heat_common::idx(0, nloc + 1, side)], side, MPI_DOUBLE, up, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        const double b = MPI_Wtime();

        // Same arithmetic, in the same order, as the serial code.
        for (int l = 1; l <= nloc; ++l) {
            for (int i = 1; i < N; ++i) {
                const std::size_t p = heat_common::idx(i, l, side);
                const double neighbors = old_grid[p - 1] + old_grid[p + 1]
                                       + old_grid[p - side] + old_grid[p + side];
                new_grid[p] = (1.0 - 4.0 * r) * old_grid[p] + r * neighbors;
            }
        }
        old_grid.swap(new_grid);
        const double c = MPI_Wtime();
        t_comm += b - a;
        t_comp += c - b;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    const double t_total_local = MPI_Wtime() - t_start;

    double t_total = 0.0, t_comp_max = 0.0, t_comm_max = 0.0;
    MPI_Reduce(&t_total_local, &t_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_comp, &t_comp_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_comm, &t_comm_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Gather only when needed (output file or error norms): this is not timed.
    const bool need_global = !opt.output_path.empty() || opt.test_case == "V";
    std::vector<double> global;
    if (need_global) {
        std::vector<int> counts, displs;
        if (rank == 0) {
            global.assign(static_cast<std::size_t>(side) * side, 0.0);
            heat_common::initialize(global, N, opt.test_case);  // fixes rows 0 and N
            counts.resize(P);
            displs.resize(P);
            for (int q = 0; q < P; ++q) {
                counts[q] = rows_of(q) * side;
                displs[q] = first_row_of(q) * side;
            }
        }
        MPI_Gatherv(&old_grid[heat_common::idx(0, 1, side)], nloc * side, MPI_DOUBLE,
                    rank == 0 ? global.data() : nullptr, counts.data(), displs.data(),
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    int rc = 0;
    if (rank == 0) {
        try {
            if (!opt.output_path.empty()) heat_common::write_binary(opt.output_path, global);
            std::cout << std::setprecision(17)
                      << "decomp=strips P=" << P << " N=" << N
                      << " steps=" << ts.steps
                      << " r=" << ts.effective_r
                      << " dt=" << ts.dt
                      << " T=" << ts.final_time
                      << " t_total=" << t_total
                      << " t_comp=" << t_comp_max
                      << " t_comm=" << t_comm_max;
            if (opt.test_case == "V") {
                double err_inf = 0.0, err_l2_h = 0.0;
                heat_common::error_norms(global, N, ts.final_time, err_inf, err_l2_h);
                std::cout << " err_inf=" << err_inf << " err_l2_h=" << err_l2_h;
            }
            std::cout << '\n';
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << '\n';
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    return rc;
}

}  // namespace heat_strips

#endif  // HEAT_HAVE_MPI