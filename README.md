# Parallel Needleman-Wunsch

This project implements the Needleman-Wunsch global sequence alignment algorithm in C and evaluates four execution models:

- Sequential baseline
- Pthreads wavefront parallelization
- OpenMP wavefront parallelization
- MPI master-worker task distribution

The benchmark workflow is fully scriptable: it compiles all implementations, generates reproducible input data, runs the benchmark grid, and writes raw outputs plus summary CSV files.

## Repository Structure

```text
.
|-- Makefile
|-- README.md
|-- run_experiments.sh
|-- data/
|   |-- input_generator.c
|   `-- generated/                  generated benchmark inputs
|-- scripts/
|   `-- run_experiments.sh           benchmark driver
|-- src/
|   |-- sequential/main_seq.c
|   |-- pthreads/main_pth.c
|   |-- openmp/main_omp.c
|   `-- mpi/main_mpi.c
`-- results/
    |-- benchmark_results.csv        summary statistics
    |-- benchmark_raw_runs.csv       per-repetition measurements
    |-- benchmark_configurations.csv expected benchmark grid
    `-- raw_outputs/                 full program outputs and logs
```

All commands use paths relative to the project root.

## Requirements

Required:

```text
GCC with C99 support
GNU Make
Bash
AWK
POSIX Pthreads
OpenMP runtime support
OpenMPI >= 4.0 or MPICH equivalent
```

Install dependencies on Debian, Ubuntu, or Kali:

```bash
sudo apt update
sudo apt install -y build-essential openmpi-bin libopenmpi-dev
```

Optional debugging tools:

```bash
sudo apt install -y gdb valgrind
```

## Build

Build all implementations:

```bash
make clean
make all
```

Build individual targets:

```bash
make sequential
make pthreads
make openmp
make mpi
make input_generator
```

Generated binaries:

```text
bin/main_seq
bin/main_pth
bin/main_omp
bin/main_mpi
bin/input_generator
```

Equivalent manual compile commands:

```bash
mkdir -p bin
gcc -std=c99 -Wall -Wextra -pedantic -O0 -g src/sequential/main_seq.c -o bin/main_seq
gcc -std=c99 -Wall -Wextra -pedantic -O0 -g src/pthreads/main_pth.c -o bin/main_pth -pthread
gcc -std=c99 -Wall -Wextra -pedantic -O0 -g src/openmp/main_omp.c -o bin/main_omp -fopenmp
mpicc -std=c99 -Wall -Wextra -pedantic -O0 -g src/mpi/main_mpi.c -o bin/main_mpi
gcc -std=c99 -Wall -Wextra -pedantic -O0 -g data/input_generator.c -o bin/input_generator
```

## Quick Correctness Check

Run a small known alignment:

```bash
printf "ACGT\nACGT\n" | ./bin/main_seq
printf "ACGT\nACGT\n" | ./bin/main_pth --threads 4
printf "ACGT\nACGT\n" | OMP_NUM_THREADS=4 ./bin/main_omp
printf "1\nACGT\nACGT\n" | mpirun -np 4 ./bin/main_mpi
```

Expected result:

```text
Alignment score: 4
Aligned sequence A: ACGT
Aligned sequence B: ACGT
```

Reference example:

```bash
printf "GATTACA\nGCATGCG\n" | ./bin/main_seq
printf "GATTACA\nGCATGCG\n" | ./bin/main_pth --threads 4
printf "GATTACA\nGCATGCG\n" | OMP_NUM_THREADS=4 ./bin/main_omp
printf "1\nGATTACA\nGCATGCG\n" | mpirun -np 4 ./bin/main_mpi
```

Expected score:

```text
Alignment score: 0
```

## Running the Benchmark

The root-level benchmark command is:

```bash
chmod +x run_experiments.sh scripts/run_experiments.sh
./run_experiments.sh
```

The root script delegates to `scripts/run_experiments.sh`.

The benchmark driver automatically:

- Compiles all implementations.
- Generates deterministic DNA sequence inputs.
- Runs sequential, Pthreads, OpenMP, and MPI configurations.
- Repeats each configuration for statistical stability.
- Saves full program output.
- Saves raw measurements and summary statistics.

Default benchmark configuration:

```text
Input sizes:       100, 500, 1000, 5000, 10000, 50000
Worker counts:     1, 2, 4, 8
Pairs per size:    3
Random seed:       2026
Repetitions:       10
Warm-up runs:      1
MPI oversubscribe: enabled for process counts greater than 4
```

Full reproducibility run:

```bash
PAIRS_PER_SIZE=3 SEED=2026 REPETITIONS=10 WARMUP_RUNS=1 ./run_experiments.sh
```

Quick dry run:

```bash
INPUT_SIZES_OVERRIDE="100" THREAD_COUNTS_OVERRIDE="1 2" REPETITIONS=2 PAIRS_PER_SIZE=1 ./run_experiments.sh
```

Use dry runs only to validate the environment or script behavior. Use the full benchmark configuration for final reported results.

## Benchmark Options

The benchmark script accepts environment variables:

| Variable | Default | Description |
|---|---:|---|
| `PAIRS_PER_SIZE` | `3` | Number of sequence pairs generated per input size |
| `SEED` | `2026` | Random seed for input generation |
| `REPETITIONS` | `10` | Measured repetitions per configuration |
| `WARMUP_RUNS` | `1` | Warm-up runs before measured repetitions |
| `MPI_OVERSUBSCRIBE_ABOVE` | `4` | Adds `--oversubscribe` above this MPI process count |
| `INPUT_SIZES_OVERRIDE` | unset | Space-separated replacement list of input sizes |
| `THREAD_COUNTS_OVERRIDE` | unset | Space-separated replacement list of worker counts |

Example:

```bash
PAIRS_PER_SIZE=5 SEED=1234 REPETITIONS=5 ./run_experiments.sh
```

## Output Files

The benchmark creates all output under `results/`.

| Path | Description |
|---|---|
| `results/raw_outputs/compile.log` | Compilation log and benchmark settings |
| `results/raw_outputs/*.txt` | Full stdout/stderr for individual runs |
| `results/raw_outputs/*.times` | Per-pair timings for sequential, Pthreads, and OpenMP runs |
| `results/benchmark_raw_runs.csv` | One row per measured repetition |
| `results/benchmark_results.csv` | Aggregated timing statistics, speedup, and efficiency |
| `results/benchmark_configurations.csv` | Benchmark grid that was executed |

Summary CSV schema:

```csv
Algorithm,InputSize,Threads,Processes,MeanTime,StdDev,MedianTime,MinTime,MaxTime,OutliersRemoved,ValidRuns,Speedup,Efficiency,Notes
```

Raw run CSV schema:

```csv
RunID,Algorithm,InputSize,Threads,Processes,PairsPerSize,ExecutionTime,Status,Notes
```

## Validating Results

Check for failed repetitions:

```bash
awk -F, '$8 != "Status" && $8 != "OK" { print }' results/benchmark_raw_runs.csv
```

No output means all measured repetitions completed successfully.

Check for incomplete summary rows:

```bash
awk -F, '$11 == 0 || $5 == "NA" { print }' results/benchmark_results.csv
```

No output means every summary row has valid timing data.

Review the executed benchmark grid:

```bash
cat results/benchmark_configurations.csv
```

## Input Generation

The input generator creates random DNA sequences over:

```text
A C G T
```

Manual generation:

```bash
make input_generator
mkdir -p data/generated
./bin/input_generator data/generated 3 2026
```

Input file format:

```text
number_of_pairs
sequence_A_pair_0
sequence_B_pair_0
sequence_A_pair_1
sequence_B_pair_1
...
```

Generated files use this naming pattern:

```text
data/generated/input_len_<length>_pairs_<count>.txt
```

## Manual Program Usage

Sequential:

```bash
printf "ACGT\nACGT\n" | ./bin/main_seq
```

Pthreads:

```bash
printf "ACGT\nACGT\n" | ./bin/main_pth --threads 4
```

OpenMP:

```bash
printf "ACGT\nACGT\n" | OMP_NUM_THREADS=4 ./bin/main_omp
```

MPI:

```bash
printf "1\nACGT\nACGT\n" | mpirun -np 4 ./bin/main_mpi
```

MPI with generated input:

```bash
mpirun -np 4 ./bin/main_mpi < data/generated/input_len_100_pairs_3.txt
```

MPI with oversubscription:

```bash
mpirun --oversubscribe -np 8 ./bin/main_mpi < data/generated/input_len_100_pairs_3.txt
```

## Reproducibility Notes

For reproducible reporting, record:

```text
CPU model
Core and thread count
RAM size
Operating system and kernel
GCC version
MPI implementation and version
OpenMP support
Compiler flags
Input seed
Pairs per size
Repetition count
```

Suggested commands:

```bash
lscpu
free -h
uname -a
gcc --version
mpicc --version
mpirun --version
```

Performance notes:

- Small inputs may show poor parallel speedup because synchronization and startup overhead dominate.
- MPI process counts above available physical cores may require oversubscription and should be interpreted carefully.
- The `50000` input can require very large memory because full score and traceback matrices are stored.
- If a benchmark run is interrupted, rerun the script to regenerate a complete result set.

## Troubleshooting

`make: command not found`

```bash
sudo apt install -y build-essential
```

`mpicc: command not found` or `mpirun: command not found`

```bash
sudo apt install -y openmpi-bin libopenmpi-dev
```

OpenMPI reports insufficient slots for 8 processes:

```bash
mpirun --oversubscribe -np 8 ./bin/main_mpi < data/generated/input_len_100_pairs_3.txt
```

The benchmark script applies `--oversubscribe` automatically for MPI process counts greater than 4 by default.

Benchmark run is too long:

```bash
INPUT_SIZES_OVERRIDE="100 500 1000" REPETITIONS=3 ./run_experiments.sh
```

Use shortened runs for environment checks only. Final reported results should use the documented full configuration.

Very large input fails:

Needleman-Wunsch stores a full `(m + 1) x (n + 1)` dynamic programming matrix. The `50000` case may exceed available memory on small VMs or lab machines.
