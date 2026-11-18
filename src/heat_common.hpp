#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace heat_common {

inline std::size_t idx(int i, int j, int side) {
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(side)
         + static_cast<std::size_t>(i);
}

// Initial value of node (i,j) for a given case. Shared by Serial and MPI versions
// so that all decompositions start from exactly the same field.
inline double initial_value(int i, int j, int N, const std::string& test_case) {
    if (test_case == "1") {
        return (i == 0 || i == N || j == 0 || j == N) ? 1.0 : 0.0;
    }
    if (test_case == "2") {
        // Top/bottom edges take precedence at the four corners.
        return (j == 0 || j == N) ? 1.0 : 0.0;
    }
    if (i > 0 && i < N && j > 0 && j < N) {
        const double pi = std::acos(-1.0);
        const double x = static_cast<double>(i) / N;
        const double y = static_cast<double>(j) / N;
        return std::sin(pi * x) * std::sin(pi * y);
    }
    return 0.0;
}

inline void initialize(std::vector<double>& grid, int N, const std::string& test_case) {
    const int side = N + 1;
    for (int j = 0; j <= N; ++j) {
        for (int i = 0; i <= N; ++i) {
            grid[idx(i, j, side)] = initial_value(i, j, N, test_case);
        }
    }
}

inline void write_binary(const std::string& path, const std::vector<double>& grid) {
    const std::uint16_t probe = 1;
    if (*reinterpret_cast<const unsigned char*>(&probe) != 1) {
        throw std::runtime_error("output requires a little-endian machine");
    }
    if (sizeof(double) != 8) {
        throw std::runtime_error("output requires 64-bit doubles");
    }
    std::ofstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot open output file: " + path);
    stream.write(reinterpret_cast<const char*>(grid.data()),
                 static_cast<std::streamsize>(grid.size() * sizeof(double)));
    if (!stream) throw std::runtime_error("failed to write output file: " + path);
}

inline void error_norms(const std::vector<double>& grid, int N, double t,
                 double& err_inf, double& err_l2_h) {
    const int side = N + 1;
    const double pi = std::acos(-1.0);
    const double factor = std::exp(-2.0 * pi * pi * t);
    double sum_sq = 0.0;
    err_inf = 0.0;
    for (int j = 0; j <= N; ++j) {
        for (int i = 0; i <= N; ++i) {
            double exact = 0.0;
            if (i > 0 && i < N && j > 0 && j < N) {
                exact = std::sin(pi * static_cast<double>(i) / N)
                      * std::sin(pi * static_cast<double>(j) / N) * factor;
            }
            const double error = grid[idx(i, j, side)] - exact;
            err_inf = std::max(err_inf, std::abs(error));
            sum_sq += error * error;
        }
    }
    err_l2_h = std::sqrt(sum_sq) / static_cast<double>(N);
}

struct Options {
    std::string test_case;
    int N = 0;
    double requested_r = 0.0;
    double final_time = 0.0;
    long long steps = 0;
    bool has_T = false;
    bool has_steps = false;
    std::string output_path;
    std::string decomp = "serial";
};

inline void usage(const char* name) {
    std::cout << "Usage: " << name
              << " --case {1|2|V} --N N --r r (--steps K | --T T)"
                 " [--out file.bin] [--decomp serial]\n";
}

inline int parse_int(const std::string& value, const std::string& flag) {
    std::size_t used = 0;
    int result = std::stoi(value, &used);
    if (used != value.size()) throw std::invalid_argument("invalid " + flag + ": " + value);
    return result;
}

inline long long parse_long(const std::string& value, const std::string& flag) {
    std::size_t used = 0;
    long long result = std::stoll(value, &used);
    if (used != value.size()) throw std::invalid_argument("invalid " + flag + ": " + value);
    return result;
}

inline double parse_double(const std::string& value, const std::string& flag) {
    std::size_t used = 0;
    double result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result)) {
        throw std::invalid_argument("invalid " + flag + ": " + value);
    }
    return result;
}

inline Options parse_options(int argc, char** argv) {
    Options opt;
    for (int k = 1; k < argc; k += 2) {
        const std::string flag(argv[k]);
        if (flag == "--help" || flag == "-h") {
            usage(argv[0]);
            std::exit(0);
        }
        if (k + 1 >= argc) throw std::invalid_argument("missing value for " + flag);
        const std::string value(argv[k + 1]);
        if (flag == "--case") opt.test_case = value;
        else if (flag == "--N") opt.N = parse_int(value, flag);
        else if (flag == "--r") opt.requested_r = parse_double(value, flag);
        else if (flag == "--T") {
            opt.has_T = true;
            opt.final_time = parse_double(value, flag);
        } else if (flag == "--steps") {
            opt.has_steps = true;
            opt.steps = parse_long(value, flag);
        } else if (flag == "--out") opt.output_path = value;
        else if (flag == "--decomp") opt.decomp = value;
        else throw std::invalid_argument("unknown option: " + flag);
    }

    if (opt.test_case != "1" && opt.test_case != "2" && opt.test_case != "V") {
        throw std::invalid_argument("--case must be 1, 2 or V");
    }
    if (opt.N < 2) throw std::invalid_argument("--N must be at least 2");
    if (!(opt.requested_r > 0.0 && opt.requested_r <= 0.25)) {
        throw std::invalid_argument("unstable r: require 0 < r <= 0.25 for the 2D explicit heat scheme");
    }
    if (opt.has_T == opt.has_steps) {
        throw std::invalid_argument("specify exactly one of --T or --steps");
    }
    if (opt.has_T && opt.final_time <= 0.0) {
        throw std::invalid_argument("--T must be positive");
    }
    if (opt.has_steps && opt.steps <= 0) {
        throw std::invalid_argument("--steps must be positive");
    }
    if (opt.decomp != "serial" && opt.decomp != "strips" && opt.decomp != "blocks") {
        throw std::invalid_argument("--decomp must be serial, strips or blocks");
    }
    return opt;
}

struct TimeSetup {
    long long steps;
    double dt;
    double effective_r;
    double final_time;
};

inline TimeSetup make_time_setup(const Options& opt) {
    const double h = 1.0 / static_cast<double>(opt.N);
    const double dt_max = opt.requested_r * h * h; // alpha = 1

    long long steps = opt.steps;
    double dt = dt_max;
    if (opt.has_T) {
        const double ratio = opt.final_time / dt_max;
        if (!std::isfinite(ratio) || ratio > static_cast<double>(std::numeric_limits<long long>::max())) {
            throw std::invalid_argument("--T implies too many time steps");
        }
        steps = std::max(1LL, static_cast<long long>(std::ceil(ratio)));
        dt = opt.final_time / static_cast<double>(steps);
        // Guard against a rare floating-point rounding edge case.
        if (dt > dt_max) {
            ++steps;
            dt = opt.final_time / static_cast<double>(steps);
        }
    }
    const double effective_r = dt / (h * h);
    const double final_time = opt.has_T ? opt.final_time : steps * dt;
    return {steps, dt, effective_r, final_time};
}

}  // namespace heat_common
