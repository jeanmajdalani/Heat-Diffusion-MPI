// Task 6.7 : faire casser le schema explicite exprès (r > 1/4).
// Programme separe : le solveur principal refuse r > 0.25.
// Usage : unstable_demo --r 0.2625 --N 64 --steps 200 --case 1 --snap 0,50,100,200 --outdir results/unstable
// Sorties : <outdir>/maxabs_r<r>.csv (step,max|u|,max|u_interieur|) et snap_r<r>_step<n>.bin
//           ((N+1)^2 doubles little-endian, ligne par ligne, i le plus rapide, comme le solveur).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static std::vector<int> parse_list(const std::string& s) {
    std::vector<int> v; std::stringstream ss(s); std::string t;
    while (std::getline(ss, t, ',')) if (!t.empty()) v.push_back(std::stoi(t));
    return v;
}

int main(int argc, char** argv) {
    int N = 64, steps = 200, test_case = 1;
    double r = 0.2625;                       // 1.05 * 0.25
    std::string outdir = "results/unstable", snap_arg = "0,50,100,150,200";
    for (int a = 1; a + 1 < argc; a += 2) {
        std::string f = argv[a], v = argv[a + 1];
        if (f == "--N") N = std::stoi(v); else if (f == "--r") r = std::stod(v);
        else if (f == "--steps") steps = std::stoi(v); else if (f == "--case") test_case = std::stoi(v);
        else if (f == "--snap") snap_arg = v; else if (f == "--outdir") outdir = v;
        else { std::cerr << "option inconnue : " << f << "\n"; return 1; }
    }
    if (N < 4 || steps < 1 || (test_case != 1 && test_case != 2)) { std::cerr << "arguments invalides\n"; return 1; }
    std::vector<int> snaps = parse_list(snap_arg);
    const int side = N + 1;
    auto idx = [side](int i, int j) { return static_cast<std::size_t>(j) * side + i; };

    // Conditions initiales : memes que le solveur (bords chauds = 1, interieur = 0).
    std::vector<double> u(side * side, 0.0), v(side * side, 0.0);
    for (int j = 0; j <= N; ++j) for (int i = 0; i <= N; ++i) {
        bool hot = (test_case == 1) ? (i == 0 || i == N || j == 0 || j == N) : (j == 0 || j == N);
        if (hot) u[idx(i, j)] = v[idx(i, j)] = 1.0;
    }
    std::ostringstream tag; tag << "r" << r;
    std::ofstream csv(outdir + "/maxabs_" + tag.str() + ".csv");
    if (!csv) { std::cerr << "dossier introuvable : " << outdir << " (faites mkdir -p)\n"; return 1; }
    csv << "step,max_abs,max_abs_interior\n";

    auto record = [&](int step) {
        double m = 0, mi = 0;
        for (int j = 0; j <= N; ++j) for (int i = 0; i <= N; ++i) {
            double a = std::fabs(u[idx(i, j)]); m = std::max(m, a);
            if (i > 0 && i < N && j > 0 && j < N) mi = std::max(mi, a);
        }
        csv << step << ',' << m << ',' << mi << '\n';
        if (std::find(snaps.begin(), snaps.end(), step) != snaps.end()) {
            std::ostringstream nm; nm << outdir << "/snap_" << tag.str() << "_step" << step << ".bin";
            std::ofstream out(nm.str(), std::ios::binary);
            out.write(reinterpret_cast<const char*>(u.data()), sizeof(double) * u.size());
        }
        return mi;
    };

    double last = record(0);
    for (int n = 1; n <= steps; ++n) {
        for (int j = 1; j < N; ++j) for (int i = 1; i < N; ++i) {
            std::size_t p = idx(i, j);
            double nb = u[p - 1] + u[p + 1] + u[p - side] + u[p + side];
            v[p] = (1.0 - 4.0 * r) * u[p] + r * nb;
        }
        u.swap(v);                               // les bords restent fixes (Dirichlet)
        last = record(n);
        if (!std::isfinite(last)) { std::cout << "overflow au pas " << n << "\n"; break; }
    }
    std::printf("r=%.5f N=%d steps<=%d  max|u| interieur final = %.6e\n", r, N, steps, last);
    return 0;
}
