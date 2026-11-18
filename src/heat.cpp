#include "heat_common.hpp"
#ifdef HEAT_HAVE_MPI
#include "heat_strips.hpp"
#endif
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using heat_common::Options;
using heat_common::parse_options;
using heat_common::usage;

using heat_common::idx;
using heat_common::initialize;
using heat_common::write_binary;
using heat_common::error_norms;

} // namespace

int main(int argc, char** argv) {
    try {
        const Options opt = parse_options(argc, argv);
        const int side = opt.N + 1;
        const auto timing = heat_common::make_time_setup(opt);
        
        if (opt.decomp != "serial") {
#ifdef HEAT_HAVE_MPI
            if (opt.decomp == "blocks") {
                throw std::invalid_argument("--decomp blocks is not implemented yet");
            }
            MPI_Init(&argc, &argv);
            const int rc = heat_strips::run(opt, timing);
            MPI_Finalize();
            return rc;
#else
            throw std::invalid_argument("this build has no MPI support: install Open MPI and rebuild");
#endif
        }
        const long long steps = timing.steps;
        const double dt = timing.dt;
        const double effective_r = timing.effective_r;
        const double final_time = timing.final_time;

        const std::size_t size = static_cast<std::size_t>(side) * static_cast<std::size_t>(side);
        std::vector<double> old_grid(size, 0.0);
        std::vector<double> new_grid(size, 0.0);
        initialize(old_grid, opt.N, opt.test_case);
        initialize(new_grid, opt.N, opt.test_case);

        const auto start = std::chrono::steady_clock::now();
        for (long long step = 0; step < steps; ++step) {
            for (int j = 1; j < opt.N; ++j) {
                for (int i = 1; i < opt.N; ++i) {
                    const std::size_t p = idx(i, j, side);
                    const double neighbors = old_grid[p - 1] + old_grid[p + 1]
                                           + old_grid[p - side] + old_grid[p + side];
                    new_grid[p] = (1.0 - 4.0 * effective_r) * old_grid[p]
                                + effective_r * neighbors;
                }
            }
            old_grid.swap(new_grid);
        }
        const auto stop = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(stop - start).count();

        if (!opt.output_path.empty()) write_binary(opt.output_path, old_grid);
        std::cout << std::setprecision(17)
                  << "decomp=serial P=1 N=" << opt.N
                  << " steps=" << steps
                  << " r=" << effective_r
                  << " dt=" << dt
                  << " T=" << final_time
                  << " t_total=" << elapsed
                  << " t_comp=" << elapsed
                  << " t_comm=0";
        if (opt.test_case == "V") {
            double err_inf = 0.0;
            double err_l2_h = 0.0;
            error_norms(old_grid, opt.N, final_time, err_inf, err_l2_h);
            std::cout << " err_inf=" << err_inf << " err_l2_h=" << err_l2_h;
        }
        std::cout << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        usage(argv[0]);
        return 1;
    }
}
