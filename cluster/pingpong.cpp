// Ping-pong between rank 0 and rank 1, message sizes 8 bytes .. 16 MiB.
// Prints CSV: bytes,reps,t_oneway_s,bandwidth_MBps. Run with exactly 2 processes.
// Fit t_msg = t_s + m * t_w on the output (latency t_s, 1/bandwidth t_w).
#include <mpi.h>
#include <cstdio>
#include <vector>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
        if (rank == 0) std::fprintf(stderr, "run with exactly 2 processes\n");
        MPI_Finalize();
        return 1;
    }
    const int peer = 1 - rank;
    std::vector<char> buf(16 * 1024 * 1024, 1);
    if (rank == 0) std::printf("bytes,reps,t_oneway_s,bandwidth_MBps\n");
    for (long long bytes = 8; bytes <= 16LL * 1024 * 1024; bytes *= 2) {
        const int reps = bytes <= 65536 ? 1000 : (bytes <= 1048576 ? 100 : 20);
        // warm-up
        for (int k = 0; k < 5; ++k) {
            if (rank == 0) {
                MPI_Send(buf.data(), (int)bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD);
                MPI_Recv(buf.data(), (int)bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(buf.data(), (int)bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(buf.data(), (int)bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        const double t0 = MPI_Wtime();
        for (int k = 0; k < reps; ++k) {
            if (rank == 0) {
                MPI_Send(buf.data(), (int)bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD);
                MPI_Recv(buf.data(), (int)bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(buf.data(), (int)bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(buf.data(), (int)bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD);
            }
        }
        const double t = MPI_Wtime() - t0;
        if (rank == 0) {
            const double one_way = t / (2.0 * reps);
            std::printf("%lld,%d,%.9g,%.3f\n", bytes, reps, one_way, bytes / one_way / 1e6);
        }
    }
    MPI_Finalize();
    return 0;
}