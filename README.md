# Parallel Needleman-Wunsch Benchmarking Project

This project implements and benchmarks the Needleman-Wunsch global sequence alignment algorithm using four C99 implementations:

- Sequential baseline
- Pthreads wavefront parallelization
- OpenMP wavefront parallelization
- MPI master-worker task farming

The repository is designed for reproducible university evaluation. All build and experiment commands use relative paths from the project root.

## Project Layout

```text
.
|-- Makefile
|-- README.md
|-- run_experiments.sh
|-- BENCHMARKING_METHODOLOGY.md
|-- NW_Technical_Research.md
|-- data/
|   |-- input_generator.c
|   `-- generated/                 created by the experiment script
|-- scripts/
|   `-- run_experiments.sh
|-- src/
|   |-- sequential/main_seq.c
|   |-- pthreads/main_pth.c
|   |-- openmp/main_omp.c
|   `-- mpi/main_mpi.c
`-- results/
    |-- benchmark_results.csv       created by the experiment script
    |-- benchmark_raw_runs.csv      created by the experiment script
    |-- benchmark_configurations.csv created by the experiment script
    `-- raw_outputs/                raw program outputs and compile log
```

## Dependencies

Required on the university lab/cluster machine:

```text
GCC with C99 support
GNU Make
Bash
AWK
POSIX Pthreads
OpenMP runtime support
MPI implementation, preferably OpenMPI >= 4.0 or MPICH equivalent
```

On Debian, Ubuntu, or Kali:

```bash
sudo apt update
sudo apt install -y build-essential openmpi-bin libopenmpi-dev
```

Optional debugging tools:

```bash
sudo apt install -y gdb valgrind
```

## Compilation

The recommended way to compile all implementations is:

```bash
make clean
make all
```

The Makefile builds these binaries:

```text
bin/main_seq
bin/main_pth
bin/main_omp
bin/main_mpi
bin/input_generator
```

Exact compile commands used by the Makefile:

```bash
mkdir -p bin
gcc -std=c99 -Wall -Wextra -pedantic -O0 -g src/sequential/main_seq.c -o bin/main_seq
gcc -std=c99 -Wall -Wextra -pedantic -O0 -g src/pthreads/main_pth.c -o bin/main_pth -pthread
gcc -std=c99 -Wall -Wextra -pedantic -O0 -g src/openmp/main_omp.c -o bin/main_omp -fopenmp
mpicc -std=c99 -Wall -Wextra -pedantic -O0 -g src/mpi/main_mpi.c -o bin/main_mpi
gcc -std=c99 -Wall -Wextra -pedantic -O0 -g data/input_generator.c -o bin/input_generator
```

Individual targets:

```bash
make sequential
make pthreads
make openmp
make mpi
make input_generator
```

## Manual Correctness Smoke Tests

Run these from the project root after compilation:

```bash
printf "ACGT\nACGT\n" | ./bin/main_seq
printf "ACGT\nACGT\n" | ./bin/main_pth --threads 4
printf "ACGT\nACGT\n" | OMP_NUM_THREADS=4 ./bin/main_omp
printf "1\nACGT\nACGT\n" | mpirun -np 4 ./bin/main_mpi
```

Expected result for all four:

```text
Alignment score: 4
Aligned sequence A: ACGT
Aligned sequence B: ACGT
```

Known reference case:

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

## Reproducing the Benchmark Logbook

The TA should run the benchmark from the project root:

```bash
chmod +x run_experiments.sh scripts/run_experiments.sh
./run_experiments.sh
```

Equivalent command:

```bash
bash scripts/run_experiments.sh
```

The script automatically:

- Creates required output directories.
- Compiles the input generator and all four implementations.
- Generates deterministic input files.
- Runs the documented benchmark grid.
- Saves raw program output under `results/raw_outputs/`.
- Saves raw repetition data to `results/benchmark_raw_runs.csv`.
- Saves summary statistics to `results/benchmark_results.csv`.
- Saves the expected configuration grid to `results/benchmark_configurations.csv`.

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

The default run is the configuration that should be used for the final logbook and paper tables. If any parameter is changed, regenerate the logbook and report the changed parameter.

## Experiment Script Options

The script is controlled through environment variables:

```bash
PAIRS_PER_SIZE=3 SEED=2026 REPETITIONS=10 WARMUP_RUNS=1 ./run_experiments.sh
```

Available options:

| Variable | Default | Meaning |
|---|---:|---|
| `PAIRS_PER_SIZE` | `3` | Number of sequence pairs generated for each input size |
| `SEED` | `2026` | Random seed used by the input generator |
| `REPETITIONS` | `10` | Measured repetitions per configuration |
| `WARMUP_RUNS` | `1` | Unrecorded warm-up runs per configuration |
| `MPI_OVERSUBSCRIBE_ABOVE` | `4` | Adds `--oversubscribe` for MPI process counts above this value |
| `INPUT_SIZES_OVERRIDE` | unset | Space-separated replacement list of input sizes |
| `THREAD_COUNTS_OVERRIDE` | unset | Space-separated replacement list of worker counts |

Quick dry run for checking the script:

```bash
INPUT_SIZES_OVERRIDE="100" THREAD_COUNTS_OVERRIDE="1 2" REPETITIONS=2 PAIRS_PER_SIZE=1 ./run_experiments.sh
```

Full reproducibility run:

```bash
PAIRS_PER_SIZE=3 SEED=2026 REPETITIONS=10 WARMUP_RUNS=1 ./run_experiments.sh
```

## Output Files

All result paths are relative to the project root.

| File or Directory | Purpose |
|---|---|
| `results/raw_outputs/compile.log` | Build log and experiment settings |
| `results/raw_outputs/*.txt` | Raw stdout/stderr for individual program runs |
| `results/raw_outputs/*.times` | Per-pair times for sequential, Pthreads, and OpenMP runs |
| `results/benchmark_raw_runs.csv` | One row per measured repetition |
| `results/benchmark_results.csv` | Mean, standard deviation, median, min, max, speedup, efficiency |
| `results/benchmark_configurations.csv` | Expected benchmark grid |

Summary CSV columns:

```csv
Algorithm,InputSize,Threads,Processes,MeanTime,StdDev,MedianTime,MinTime,MaxTime,OutliersRemoved,ValidRuns,Speedup,Efficiency,Notes
```

Raw repetition CSV columns:

```csv
RunID,Algorithm,InputSize,Threads,Processes,PairsPerSize,ExecutionTime,Status,Notes
```

## Verifying Completion

After the benchmark finishes, check for failed runs:

```bash
awk -F, '$8 != "Status" && $8 != "OK" { print }' results/benchmark_raw_runs.csv
```

No output means every recorded repetition completed successfully.

Check for incomplete summary rows:

```bash
awk -F, '$11 == 0 || $5 == "NA" { print }' results/benchmark_results.csv
```

No output means every summary row has valid timing data.

Check the produced grid:

```bash
cat results/benchmark_configurations.csv
```

## Manual Input Generation

The input generator creates files in the format required by the MPI program:

```text
number_of_pairs
sequence_A_pair_0
sequence_B_pair_0
sequence_A_pair_1
sequence_B_pair_1
...
```

Manual generation:

```bash
make input_generator
mkdir -p data/generated
./bin/input_generator data/generated 3 2026
```

Generated files:

```text
data/generated/input_len_100_pairs_3.txt
data/generated/input_len_500_pairs_3.txt
data/generated/input_len_1000_pairs_3.txt
data/generated/input_len_5000_pairs_3.txt
data/generated/input_len_10000_pairs_3.txt
data/generated/input_len_50000_pairs_3.txt
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

MPI with oversubscription for 8 processes:

```bash
mpirun --oversubscribe -np 8 ./bin/main_mpi < data/generated/input_len_100_pairs_3.txt
```

## Reproducibility Notes for the TA

The numbers in the paper/logbook should be generated using:

```bash
PAIRS_PER_SIZE=3 SEED=2026 REPETITIONS=10 WARMUP_RUNS=1 ./run_experiments.sh
```

If the TA reruns the script on different hardware, exact times may differ, but the relative behavior should be reasonable:

- Very small inputs may show poor parallel speedup because thread/process overhead dominates.
- MPI with 8 processes may be marked `oversubscribed` on machines with fewer than 8 available slots.
- The `50000` input may require very large memory because the algorithm stores full score and traceback matrices.
- If a run is interrupted, the CSV may be incomplete. Rerun the script to regenerate the complete result set.

For reproducible reporting, include:

```text
CPU model
Core/thread count
RAM size
OS and kernel
GCC version
MPI implementation/version
OpenMP support
Compiler flags
Seed
Pairs per size
Repetitions
```

## Troubleshooting

`make: command not found`

Install GNU Make and build tools:

```bash
sudo apt install -y build-essential
```

`mpicc: command not found` or `mpirun: command not found`

Install MPI:

```bash
sudo apt install -y openmpi-bin libopenmpi-dev
```

OpenMPI refuses 8 processes because there are not enough slots:

```bash
mpirun --oversubscribe -np 8 ./bin/main_mpi < data/generated/input_len_100_pairs_3.txt
```

The experiment script already applies `--oversubscribe` for MPI process counts greater than 4 by default.

Benchmark takes too long:

```bash
INPUT_SIZES_OVERRIDE="100 500 1000" REPETITIONS=3 ./run_experiments.sh
```

Use shortened runs only for testing the script. Do not use shortened runs for the final logbook unless the paper states the changed configuration.

Very large input fails:

Needleman-Wunsch uses a full `(m + 1) x (n + 1)` dynamic programming matrix. The `50000` case may exceed VM or lab-machine memory. Report this as a memory limitation if it occurs.
