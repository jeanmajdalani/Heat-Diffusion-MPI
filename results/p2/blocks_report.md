# MPI Blocks - Report Contribution (Person 2)

## Architecture and domain decomposition

The Blocks solver distributes the interior grid nodes
among P MPI processes arranged in a two-dimensional
Cartesian topology.

MPI_Dims_create determines the process grid dimensions.
MPI_Cart_create creates a non-periodic communicator,
and MPI_Cart_shift identifies the four neighboring ranks.

Both spatial dimensions are partitioned as evenly as
possible. Remainder nodes are distributed among the
first process coordinates in each dimension.

Each process stores its owned rectangular block,
surrounded by one layer of ghost cells on all four sides.

Example with P = 4:

+-------------+-------------+
|             |             |
|   Rank 0    |   Rank 1    |
|             |             |
+-------------+-------------+
|             |             |
|   Rank 2    |   Rank 3    |
|             |             |
+-------------+-------------+

## Halo exchange

At every time step, each process exchanges boundary
values with its neighboring processes using MPI_Sendrecv.

Rows are contiguous in memory and are exchanged directly.

Columns are non-contiguous and are exchanged using
a derived datatype created with MPI_Type_vector.
No manual column packing is used for halo exchanges.

At physical boundaries, MPI_PROC_NULL prevents
communication outside the plate. The prescribed
Dirichlet boundary values remain unchanged.

After exchanging halos, each process updates its
interior nodes using the same five-point stencil
and arithmetic order as the serial implementation.

Two arrays represent the old and new time levels.
They are swapped after every iteration.

## Why corner ghost cells are unnecessary

The five-point stencil uses only the central node
and its four direct neighbors.

           Top
            |
   Left -- Center -- Right
            |
          Bottom

Diagonal nodes are not involved in the update.
Therefore, exchanging corner ghost values would
introduce unnecessary communication.

A nine-point stencil involving diagonal neighbors
would require access to corner ghost values.

## Global reconstruction

When global output or verification is required,
MPI_Gatherv collects the owned interior values.

Rank 0 reconstructs the complete global field using
the Cartesian coordinates and block dimensions.

The final binary output follows the same format
as the serial and Strips implementations.

## Correctness verification

The solver was compiled and tested locally on an
ARM64 Ubuntu 24.04 VM with four virtual CPUs,
Open MPI 4.1.6 and GCC 13.3.0.

The shared automated comparison script reported:

102/102 bit-for-bit identical results.

An independent validation campaign tested Cases 1,
2 and V with N = 8, 9, 16, 33 and 64, using valid
process counts between 1 and 8.

The independent campaign checked execution failures,
removed previous output files, imposed MPI timeouts
and compared binary outputs against the serial solver.

Result: 117/117 successful comparisons, exit code 0.

These experiments cover uneven partitions,
prime process counts and rectangular process grids.

Evidence: results/p2/blocks_validation.log

## Local performance measurements

Case 1, 500 steps, r = 0.24.
Three repetitions per configuration.
Local slot mapping without oversubscription.

The following values are median times in seconds.

| N | P | Total time | Computation | Communication |
|---|---|---|---|---|
| 512 | 1 | 0.060561 | 0.060502 | 0.000033 |
| 512 | 2 | 0.034739 | 0.032090 | 0.003126 |
| 512 | 4 | 0.020353 | 0.017190 | 0.004492 |
| 2048 | 1 | 1.236094 | 1.235366 | 0.000452 |
| 2048 | 2 | 0.690292 | 0.658386 | 0.047638 |
| 2048 | 4 | 0.714414 | 0.508514 | 0.221208 |

Computation, communication and total times are
reduced independently using the maximum across ranks.
Their reported values need not add up exactly.

For N = 512, four processes achieved a speed-up
of approximately 2.98 relative to Blocks with P = 1.

For N = 2048, the speed-up was approximately 1.79
with two processes and 1.73 with four processes.

The four-process runs at N = 2048 varied between
approximately 0.534 and 0.719 seconds.

Communication overhead, memory bandwidth and VM
scheduling are possible explanations for the
observed scaling behavior and variability.

Raw measurements:
results/p2/blocks_local.csv
results/p2/blocks_local.raw.log

## Distributed validation and performance model

The reported results are local single-VM measurements.

Cross-machine correctness testing and the comparison
between Strips and Blocks for P = 4 and 8 and
N = 512, 2048 and 8192 require the distributed
measurements and network parameters from Person 3.

These results will be integrated into the final
team performance analysis.

## Person 2 contribution

Implementation of the MPI Blocks solver, Cartesian
decomposition, derived column datatype, halo exchange,
global reconstruction, correctness verification,
local benchmarks and background exercises
3.5, 4.1, 5.3, 6.1 and 6.2.

