# Heat-Diffusion-MPI

Parallel 2D heat diffusion on the unit square with MPI over Tailscale
(Parallel Programming course project). Explicit scheme, alpha = 1,
Dirichlet boundaries.

## Status

| Part | Content | State |
| --- | --- | --- |
| 0 | Cluster (Ubuntu 24.04, Tailscale, SSH, Open MPI 4.1.6) | in progress (Person 3) |
| 1 | Serial solver + convergence study | done (observed order about 2.0, `results/convergence.csv`) |
| 2 | MPI 1D strips | done (bit-for-bit identical to serial, `tests/compare_all.sh strips`) |
| 3 | MPI 2D blocks | in progress (Person 2) |
| 4 | Performance study | in progress (Person 3), runner: `scripts/benchmark.sh` |

## Build

Requires a C++17 compiler and CMake. MPI is optional for the serial
solver and will be required for the MPI versions.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Flags: `-O3 -ffp-contract=off` (needed so that results are identical
on x86-64 and ARM).

## Run

```bash
./build/heat --case 1 --N 64 --r 0.2 --steps 500 --out case1.bin
./build/heat --case 2 --N 64 --r 0.2 --steps 500 --out case2.bin
./build/heat --case V --N 32 --r 0.2 --T 0.05 --out V32.bin
python3 tests/test_serial.py
```
## Run with MPI (strips)

```bash
mpirun -np 4 ./build/heat --case 1 --N 256 --r 0.24 --steps 500 --decomp strips
tests/compare_all.sh strips
scripts/benchmark.sh -N "512 1024" -P "1 2 4 8" -R 3 -M "node slot"
```

Strips need P <= N-1. For MPI runs, `t_total`, `t_comp` and `t_comm` are the
maximum over ranks. On the cluster, use `-H cluster/hosts -I tailscale0`
(see `cluster/hosts.example`).

## Command-line options

| Option | Meaning |
| --- | --- |
| `--case 1\|2\|V` | 1: four hot edges. 2: top and bottom edges at 1 (corners included), left and right at 0. V: verification problem sin(pi x) sin(pi y), exact solution known |
| `--N` | number of intervals per direction (grid is (N+1) x (N+1)) |
| `--r` | r = alpha*dt/h^2. Values above 0.25 are refused (unstable) |
| `--steps K` or `--T T` | exactly one of: fixed number of steps, or final time. With `--T`, dt is chosen so that T/dt is an integer and r is not exceeded |
| `--out file.bin` | optional output file for the final field |
| `--decomp serial\|strips\|blocks` | decomposition: `serial` (default), `strips` (MPI, 1D rows, run with `mpirun`), `blocks` (MPI, 2D, in progress) |

## Output

- Field file: (N+1)^2 IEEE-754 doubles, little-endian, row-major, i
  (x direction) varies fastest, no header. Boundary nodes included.
- Standard output: one `key=value` line, for example
  `decomp=serial P=1 N=64 steps=500 r=... dt=... T=... t_total=... t_comp=... t_comm=0`.
  Test V adds `err_inf` and `err_l2_h`.
- Timing covers the time-stepping loop only (no initialisation, error
  computation or output).

## Repository layout

| Path | Content |
| --- | --- |
| `src/heat.cpp` | main program (serial) |
| `src/heat_common.hpp` | options, time step, initialisation, output, error norms |
| `tests/test_serial.py` | cases 1, 2, V, stability check, convergence study |
| `scripts/`, `results/` | numerical experiments and their outputs |
| `solutions/` | written solutions of the background-notes Tasks |
| `llm_log/` | AI interaction log (required by the course) |
| `AGENTS.md` | AI logging rules |
| `src/heat_strips.hpp` | MPI strips decomposition |
| `tests/compare_all.sh` | serial vs MPI comparison |
| `scripts/benchmark.sh` | benchmark runner (repetitions, node/slot mapping) |
| `scripts/heatmap.py`, `scripts/plot_unstable.py` | heat maps, instability figure |
| `tools/unstable_demo.cpp` | Task 6.7 instability experiment (r above 1/4) |
| `cluster/` | MPI hello world, ping-pong, hostfile template |

## Verification

`python3 tests/test_serial.py` checks boundary values, binary format,
rejection of r > 0.25 and the observed convergence order (about 2.0
for N = 16..256 at T = 0.05, r = 0.2).

## Team

- Jean Majdalani: serial solver, MPI strips, verification, benchmark runner, figures, integration.
- Ali Hassoun : MPI blocks (Cartesian decomposition), block-level validation.
- Yorgo : cluster (Tailscale, SSH, hostfile), inter-VM benchmarks, performance model.

## AI use

All AI interactions are recorded in `llm_log/` (see `AGENTS.md`).
