# 2D Heat Equation — Part 1, serial C++ starter

This code implements only the **serial** version of the assignment. It must later be extended to MPI strips and MPI blocks, and tested on the required VM cluster.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Run

```bash
./build/heat --case 1 --N 64 --r 0.2 --steps 500 --out case1.bin
./build/heat --case 2 --N 64 --r 0.2 --steps 500 --out case2.bin
./build/heat --case V --N 32 --r 0.2 --T 0.05 --out V32.bin
python3 tests/test_serial.py
```

The program accepts `--case 1|2|V`, `--N`, `--r`, exactly one of `--steps` or `--T`, optional `--out`, and optional `--decomp serial`. Values `r > 0.25` are refused. With `--T`, it chooses an integer number of steps such that final time is exactly `T` and the effective `r` does not exceed the requested value.

The output file contains `(N+1)^2` IEEE-754 doubles in **little-endian, row-major** order, with `i` (horizontal coordinate) varying fastest. No header. The program prints one `key=value` summary line. Test `V` adds error norms. Timing excludes initialization, error calculation, and output.

For Case 2, the **top and bottom edges including corners** are 1, and other left/right edge nodes are 0 (corner precedence as specified by the teacher).

## Submission requirements still pending

- Complete **all non-optional background Tasks** and submit a separate PDF.
- Implement MPI Strips with `MPI_Sendrecv`, `MPI_Gatherv`, identical output to serial.
- Implement MPI Blocks, network ping-pong, hostfile, performance experiments.
- Cluster according to Part 0 of the PDF and live demo on at least 3 VMs.
- Export complete AI prompt/response transcripts to `llm_log/`. The directory included here is only a **placeholder**, not a complete log.

Read and test every line of this AI-assisted starter before submitting, and preserve the full conversation logs required by the course.
