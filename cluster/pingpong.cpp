// Ping-pong between rank 0 and rank 1, message sizes 8 bytes .. 16 MiB.
// Prints CSV: bytes,reps,t_oneway_s,bandwidth_MBps. Run with exactly 2 processes.
// Fit t_msg = t_s + m * t_w on the output (latency t_s, 1/bandwidth t_w).
//
// Usage: pingpong [--budget SECONDS] [--max-bytes M]
//   --budget S    time budget per message size.  Default 0 = the original fixed
//                 repetition counts (1000 up to 64 KiB, 100 up to 1 MiB, 20 above).
//                 With S > 0 the number of repetitions is chosen from a measured
//                 round trip so that one size takes about S seconds, between
//                 3 repetitions and the fixed count above.  Use it on slow links
//                 (on a ~0.3 MB/s relayed link the fixed counts take ~2.5 hours).
//   --max-bytes M largest message size (default 16777216).
// Rank 0 decides the number of repetitions and broadcasts it, so both ranks
// always run the same number of round trips (a mismatch would deadlock).
#include <mpi.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static void exchange(int rank, int peer, char* buf, int bytes) {
    if (rank == 0) {
        MPI_Send(buf, bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD);
        MPI_Recv(buf, bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(buf, bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buf, bytes, MPI_CHAR, peer, 0, MPI_COMM_WORLD);
    }
}

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
    double budget = 0.0;
    long long max_bytes = 16LL * 1024 * 1024;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--budget") && i + 1 < argc) budget = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--max-bytes") && i + 1 < argc) max_bytes = std::atoll(argv[++i]);
    }
    const int peer = 1 - rank;
    std::vector<char> buf(16 * 1024 * 1024, 1);
    if (rank == 0) std::printf("bytes,reps,t_oneway_s,bandwidth_MBps\n");
    for (long long bytes = 8; bytes <= max_bytes; bytes *= 2) {
        const int fixed = bytes <= 65536 ? 1000 : (bytes <= 1048576 ? 100 : 20);
        // warm-up (also gives one timing probe)
        const double w0 = MPI_Wtime();
        for (int k = 0; k < 5; ++k) exchange(rank, peer, buf.data(), (int)bytes);
        const double probe = (MPI_Wtime() - w0) / 5.0;  // one round trip
        int reps = fixed;
        if (budget > 0.0 && rank == 0)
            reps = std::min(fixed, std::max(3, (int)(budget / probe)));
        MPI_Bcast(&reps, 1, MPI_INT, 0, MPI_COMM_WORLD);   // same reps on both ranks
        MPI_Barrier(MPI_COMM_WORLD);
        const double t0 = MPI_Wtime();
        for (int k = 0; k < reps; ++k) exchange(rank, peer, buf.data(), (int)bytes);
        const double t = MPI_Wtime() - t0;
        if (rank == 0) {
            const double one_way = t / (2.0 * reps);
            std::printf("%lld,%d,%.9g,%.3f\n", bytes, reps, one_way, bytes / one_way / 1e6);
            std::fflush(stdout);   // keep partial results if a long run is interrupted
        }
    }
    MPI_Finalize();
    return 0;
}
