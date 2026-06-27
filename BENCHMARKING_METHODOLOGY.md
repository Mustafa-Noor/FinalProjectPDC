# Benchmarking Methodology for Parallel Needleman-Wunsch

## 1. Objective

This methodology defines a reproducible, IEEE-style experimental plan for evaluating four implementations of the Needleman-Wunsch global sequence alignment algorithm:

- Sequential baseline
- Pthreads wavefront implementation
- OpenMP wavefront implementation
- MPI master-worker implementation

The experiments measure execution time, strong scaling, weak scaling, speedup, efficiency, and implementation overheads. Correctness must be verified before any timing result is accepted.

## 2. Implementations Under Test

| Implementation | Parallel Model | Description |
|---|---|---|
| Sequential | Single thread | Baseline implementation used for correctness and speedup |
| Pthreads | Shared-memory threads | Manual wavefront parallelization using anti-diagonal barriers |
| OpenMP | Shared-memory threads | Compiler-directed wavefront parallelization using `parallel for` |
| MPI | Distributed processes | Master-worker task farming over independent sequence pairs |

All implementations must use the same scoring scheme:

```text
Match:    +1
Mismatch: -1
Gap:      -1
```

All implementations must use the same deterministic traceback tie-breaking policy. This ensures that alignment scores and traceback outputs can be compared directly.

## 3. Experimental Platform

Before collecting benchmark results, record the full execution environment:

```text
CPU model
Number of physical cores
Number of hardware threads
RAM size
Cache hierarchy if available
Operating system and kernel version
Compiler name and version
Compiler flags
OpenMP runtime version
MPI implementation and version
Virtual machine configuration if applicable
Random seed
Number of repetitions
Outlier removal method
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

If MPI runs with 8 processes on a VM with fewer than 8 available execution slots, use OpenMPI oversubscription:

```bash
mpirun --oversubscribe -np 8 ./bin/main_mpi
```

Oversubscribed runs should be labeled clearly. They are useful for functional scaling comparison, but they do not represent true physical-core scaling.

## 4. Inputs

Use random DNA sequences over this alphabet:

```text
A C G T
```

Use a fixed seed for reproducibility. The required input sizes are:

| Input Size | Matrix Dimensions | Approximate DP Cells |
|---:|---:|---:|
| 100 | 101 x 101 | 10,201 |
| 500 | 501 x 501 | 251,001 |
| 1000 | 1001 x 1001 | 1,002,001 |
| 5000 | 5001 x 5001 | 25,010,001 |
| 10000 | 10001 x 10001 | 100,020,001 |
| 50000 | 50001 x 50001 | 2,500,100,001 |

The `50000` case is extremely memory intensive because this project stores the full score and direction matrices. Run it only if the machine has sufficient RAM. If it cannot complete, report it as a memory-capacity limit rather than silently removing it.

Use the same generated sequence pairs for every implementation.

## 5. Worker Counts

Use the following worker counts:

```text
1, 2, 4, 8
```

Interpretation:

| Implementation | Worker Meaning |
|---|---|
| Sequential | Always 1 |
| Pthreads | Thread count |
| OpenMP | `OMP_NUM_THREADS` |
| MPI | Process count |

## 6. Correctness Requirement

Before benchmarking, each parallel implementation must be compared against the sequential baseline.

Record:

```text
Input sequence A
Input sequence B
Sequential score
Parallel score
Sequential aligned A
Parallel aligned A
Sequential aligned B
Parallel aligned B
Pass/fail status
```

A timing result is valid only if:

- The program exits successfully.
- The alignment score matches the sequential baseline.
- The traceback score matches the DP matrix score.
- The traceback output matches the baseline when deterministic tie-breaking is used.

## 7. Metrics

### 7.1 Execution Time

Execution time is the wall-clock time required to complete an alignment workload.

```text
T = execution time in seconds
```

Use implementation-appropriate timers:

| Implementation | Recommended Timer |
|---|---|
| Sequential | wall-clock timer or `clock()` if already used |
| Pthreads | `clock_gettime(CLOCK_MONOTONIC)` |
| OpenMP | `omp_get_wtime()` |
| MPI | `MPI_Wtime()` |

For fair reporting, define whether timing includes only computation or includes input/output. Use the same policy throughout the study.

### 7.2 Speedup

Speedup measures how much faster a parallel run is than the sequential baseline.

```text
Speedup(p) = T_seq / T_parallel(p)
```

Where:

```text
p = number of threads or processes
T_seq = sequential execution time
T_parallel(p) = parallel execution time using p workers
```

### 7.3 Efficiency

Efficiency measures how effectively the workers are used.

```text
Efficiency(p) = Speedup(p) / p
```

Ideal efficiency is `1.0`. Values below `1.0` indicate synchronization, communication, memory, scheduling, or load-balance overhead. Values above `1.0` may occur because of measurement noise, cache effects, or workload differences, and should be interpreted cautiously.

### 7.4 Parallel Fraction

The parallel fraction is the portion of runtime that can be accelerated by parallel workers.

Sequential portions include:

- Input reading
- Matrix allocation
- Matrix initialization
- Traceback
- Thread creation
- MPI process startup
- Synchronization barriers
- Communication and result collection

Estimate the parallel fraction from measured speedup using Amdahl's Law when useful.

## 8. Scaling Laws

### 8.1 Amdahl's Law

Amdahl's Law models strong scaling for a fixed problem size:

```text
Speedup(p) = 1 / ((1 - f) + f / p)
```

Where:

```text
f = parallel fraction
p = number of workers
```

Interpretation:

Even if most of the program is parallel, the remaining serial fraction limits maximum speedup. In this project, traceback, initialization, synchronization, and communication reduce the achievable speedup.

### 8.2 Gustafson's Law

Gustafson's Law models scaled workloads where the problem size grows with the number of workers:

```text
ScaledSpeedup(p) = p - alpha(p - 1)
```

Where:

```text
alpha = serial fraction
p = number of workers
```

Interpretation:

Needleman-Wunsch has `O(mn)` work. Larger matrices contain much more parallel work, so weak scaling can look better than strong scaling if the machine has enough memory bandwidth and capacity.

## 9. Experiment 1: Execution Time vs Input Size

### Purpose

Measure how runtime grows as sequence length increases.

### Configuration

Run all implementations using representative worker counts:

| Implementation | Worker Count |
|---|---:|
| Sequential | 1 |
| Pthreads | 4 or best available |
| OpenMP | 4 or best available |
| MPI | 4 or best available |

Use input sizes:

```text
100, 500, 1000, 5000, 10000, 50000
```

### Expected Outcome

Runtime should grow approximately quadratically with sequence length because the DP matrix has `O(n^2)` cells for equal-length sequences.

### Record

```text
Algorithm
Input size
Sequence lengths
Workers
Execution time
Score
Correctness status
Peak memory if available
```

### Repetitions

Run at least 5 repetitions. For final IEEE-style results, run 10 repetitions.

## 10. Experiment 2: Strong Scaling

### Purpose

Measure how runtime changes when the input size is fixed and worker count increases.

### Configuration

For each input size:

```text
100, 500, 1000, 5000, 10000, 50000
```

Run:

```text
Pthreads: 1, 2, 4, 8 threads
OpenMP:   1, 2, 4, 8 threads
MPI:      1, 2, 4, 8 processes
```

### Expected Outcome

Small inputs may slow down as workers increase because overhead dominates. Larger inputs should show better speedup until limited by synchronization, memory bandwidth, cache effects, or communication overhead.

### Record

```text
Algorithm
Input size
Workers
Mean execution time
Standard deviation
Minimum time
Maximum time
Speedup
Efficiency
Correctness status
```

### Repetitions

Use 10 repetitions for final reporting. Use a warm-up run before measured repetitions.

## 11. Experiment 3: Weak Scaling

### Purpose

Measure whether runtime remains stable when both the workload and worker count increase.

Needleman-Wunsch has `O(n^2)` work for equal-length sequences. To keep work per worker approximately constant:

```text
n_p = n_1 * sqrt(p)
```

### Preferred Weak Scaling Plan

| Workers | Suggested Size |
|---:|---:|
| 1 | 1000 |
| 2 | 1414 |
| 4 | 2000 |
| 8 | 2828 |

### Required-Size Approximation

If only the required input sizes may be used:

| Workers | Approximate Required Size |
|---:|---:|
| 1 | 1000 |
| 2 | 1000 or 5000 |
| 4 | 5000 |
| 8 | 10000 |

The preferred plan is more mathematically accurate. The required-size approximation is acceptable if the report clearly states that it is approximate.

### Expected Outcome

Ideal weak scaling keeps runtime approximately constant. Increasing runtime indicates synchronization overhead, cache effects, communication overhead, memory bandwidth limits, or load imbalance.

### Record

```text
Algorithm
Workers
Input size
Mean execution time
Standard deviation
Weak scaling efficiency
Correctness status
```

Weak scaling efficiency:

```text
E_weak(p) = T_1 / T_p
```

### Repetitions

Use 10 repetitions for final reporting.

## 12. Experiment 4: Pthreads vs OpenMP

### Purpose

Compare manual thread management against compiler-directed parallelism for the same wavefront algorithm.

### Configuration

Run Pthreads and OpenMP at:

```text
1, 2, 4, 8 threads
```

Use all input sizes.

### Expected Outcome

OpenMP may have lower development complexity and reasonable performance. Pthreads may expose more overhead depending on barrier cost, chunking, and thread management. Both may suffer from wavefront synchronization and memory bandwidth limits.

### Record

```text
Algorithm
Input size
Thread count
Mean execution time
Speedup
Efficiency
Correctness status
```

## 13. Experiment 5: MPI vs OpenMP

### Purpose

Compare distributed-process task farming with shared-memory threading.

### Configuration

Run:

```text
OpenMP: 1, 2, 4, 8 threads
MPI:    1, 2, 4, 8 processes
```

Use the same sequence-pair workload.

### Expected Outcome

OpenMP should usually be better for a single alignment on one shared-memory machine. MPI can perform well when there are multiple independent sequence pairs because the master can distribute separate alignments across workers.

### Record

```text
Algorithm
Input size
Workers
Execution time
Speedup
Efficiency
Number of sequence pairs
Correctness status
```

## 14. Repetitions, Averaging, and Outliers

### Repetitions

Use:

```text
Minimum:     5 repetitions
Recommended: 10 repetitions
```

For publication-quality results, use 10 measured repetitions after one unrecorded warm-up run.

### Averaging

Compute the arithmetic mean after removing outliers:

```text
Mean = sum(valid runtimes) / number of valid runtimes
```

Also report:

```text
Standard deviation
Minimum
Maximum
Median
```

### Outlier Detection

Use the interquartile range method:

```text
Q1 = 25th percentile
Q3 = 75th percentile
IQR = Q3 - Q1

Lower bound = Q1 - 1.5 * IQR
Upper bound = Q3 + 1.5 * IQR
```

Reject runtimes below the lower bound or above the upper bound. Report how many samples were rejected. If only 5 repetitions are available, use either IQR or a trimmed mean by removing the minimum and maximum values.

## 15. Overheads to Analyze

### 15.1 Barrier Overhead

Pthreads and OpenMP wavefront implementations require synchronization after every anti-diagonal. This is necessary because cells on diagonal `d + 1` depend on cells from diagonal `d`.

Expected effect:

- High overhead for small matrices
- Reduced benefit from additional threads
- Better scaling for larger matrices where each diagonal has more cells

### 15.2 Communication Overhead

MPI introduces communication between the master and workers.

Sources:

- Sending sequence pairs
- Receiving scores and alignments
- Message latency
- Master-side scheduling
- Result collection

Expected effect:

- MPI overhead is visible for small inputs
- MPI becomes more useful for batches of independent alignments

### 15.3 Load Imbalance

Load imbalance occurs when workers receive unequal amounts of useful work.

Sources:

- Early and late wavefront diagonals have fewer cells
- MPI workers may receive tasks with different sequence lengths
- Oversubscribed processes compete for CPU time

Expected effect:

- Some workers finish earlier and wait
- Efficiency decreases as worker count increases

### 15.4 False Sharing

False sharing occurs when multiple threads update different memory locations that reside on the same cache line.

Expected effect:

- No correctness failure
- Increased cache coherence traffic
- Lower Pthreads/OpenMP speedup

### 15.5 Cache Locality

Sequential row-major DP filling has good spatial locality. Wavefront traversal accesses cells along anti-diagonals, which may be less contiguous in row-major memory.

Expected effect:

- Sequential may perform surprisingly well
- Parallel wavefront may lose cache efficiency
- Large matrices may suffer more cache misses

### 15.6 Memory Bandwidth

Needleman-Wunsch stores and updates large matrices. For large inputs, performance may become limited by memory bandwidth rather than CPU arithmetic.

Expected effect:

- Speedup plateaus
- Additional threads provide little benefit
- Very large inputs may fail due to memory capacity

## 16. CSV Templates

### 16.1 Raw Repetition CSV

```csv
RunID,Algorithm,InputSize,SeqA_Length,SeqB_Length,Threads,Processes,ExecutionTime,Score,Correct,Notes
1,Sequential,1000,1000,1000,1,1,0.000000,0,Yes,
1,Pthreads,1000,1000,1000,4,1,0.000000,0,Yes,
1,OpenMP,1000,1000,1000,4,1,0.000000,0,Yes,
1,MPI,1000,1000,1000,1,4,0.000000,0,Yes,
```

### 16.2 Summary CSV

```csv
Algorithm,InputSize,Threads,Processes,MeanTime,StdDev,Median,MinTime,MaxTime,OutliersRemoved,Speedup,Efficiency,Correct
Sequential,1000,1,1,0.000000,0.000000,0.000000,0.000000,0.000000,0,1.000000,1.000000,Yes
Pthreads,1000,4,1,0.000000,0.000000,0.000000,0.000000,0.000000,0,0.000000,0.000000,Yes
OpenMP,1000,4,1,0.000000,0.000000,0.000000,0.000000,0.000000,0,0.000000,0.000000,Yes
MPI,1000,1,4,0.000000,0.000000,0.000000,0.000000,0.000000,0,0.000000,0.000000,Yes
```

### 16.3 Strong Scaling CSV

```csv
Algorithm,InputSize,Workers,T1,Tp,Speedup,Efficiency,StdDev,Correct
Pthreads,5000,1,0.000000,0.000000,1.000000,1.000000,0.000000,Yes
Pthreads,5000,2,0.000000,0.000000,0.000000,0.000000,0.000000,Yes
Pthreads,5000,4,0.000000,0.000000,0.000000,0.000000,0.000000,Yes
Pthreads,5000,8,0.000000,0.000000,0.000000,0.000000,0.000000,Yes
```

### 16.4 Weak Scaling CSV

```csv
Algorithm,Workers,InputSize,MeanTime,T1,WeakScalingEfficiency,StdDev,Correct
OpenMP,1,1000,0.000000,0.000000,1.000000,0.000000,Yes
OpenMP,2,1414,0.000000,0.000000,0.000000,0.000000,Yes
OpenMP,4,2000,0.000000,0.000000,0.000000,0.000000,Yes
OpenMP,8,2828,0.000000,0.000000,0.000000,0.000000,Yes
```

### 16.5 Correctness CSV

```csv
CaseID,SeqA,SeqB,ExpectedScore,SequentialScore,PthreadsScore,OpenMPScore,MPIScore,TracebackMatch,Pass
T001,ACGT,ACGT,4,4,4,4,4,Yes,Yes
```

## 17. Publication-Quality Graphs

Use clear axis labels, units, legends, and captions. Use logarithmic y-axes when runtime spans several orders of magnitude. Include error bars using standard deviation or confidence intervals.

### 17.1 Execution Time vs Input Size

Design:

```text
X-axis: Input size
Y-axis: Mean execution time in seconds
Scale: Log-log recommended
Lines: Sequential, Pthreads, OpenMP, MPI
Error bars: Standard deviation
```

Interpretation:

This graph shows how runtime grows with sequence length. A quadratic trend is expected. Deviations may indicate cache effects, memory bandwidth saturation, or overhead domination for small inputs.

### 17.2 Speedup vs Threads

Design:

```text
X-axis: Workers
Y-axis: Speedup
Lines: Pthreads, OpenMP, MPI
Reference line: Ideal speedup y = x
Separate plots or facets: Input sizes
```

Interpretation:

Curves close to the ideal line indicate good scaling. Flattening curves indicate synchronization, communication, memory bandwidth, or serial bottlenecks.

### 17.3 Efficiency vs Threads

Design:

```text
X-axis: Workers
Y-axis: Efficiency
Reference line: y = 1
Lines: Pthreads, OpenMP, MPI
```

Interpretation:

Efficiency near `1.0` means workers are used effectively. Declining efficiency means each additional worker contributes less useful work.

### 17.4 Strong Scaling

Design:

```text
X-axis: Workers
Y-axis: Mean execution time
Lines: One line per input size
Separate graph per implementation
```

Interpretation:

For fixed input size, runtime should decrease as workers increase. If runtime increases, overhead is greater than the benefit of parallelism. Small inputs are expected to scale poorly.

### 17.5 Weak Scaling

Design:

```text
X-axis: Workers
Y-axis: Mean execution time
Lines: Pthreads, OpenMP, MPI
```

Interpretation:

Ideal weak scaling is flat. Increasing runtime indicates growing overhead, memory bandwidth saturation, or communication cost.

### 17.6 Execution Time Comparison

Design:

```text
Graph type: Grouped bar chart
X-axis: Input size
Y-axis: Mean execution time
Bars: Sequential, Pthreads, OpenMP, MPI
Worker count: Fixed, for example 4
```

Interpretation:

This graph shows which implementation is fastest at each problem size. It is useful for identifying the crossover point where parallelism becomes beneficial.

### 17.7 Pthreads vs OpenMP

Design:

```text
X-axis: Thread count
Y-axis: Runtime or speedup
Lines: Pthreads and OpenMP
Separate plots: Selected input sizes
```

Interpretation:

This graph compares manual threading with directive-based threading. Differences may come from barrier implementation, scheduling overhead, runtime behavior, and cache locality.

### 17.8 MPI vs OpenMP

Design:

```text
X-axis: Workers
Y-axis: Runtime or speedup
Lines: MPI and OpenMP
```

Interpretation:

This graph compares process-level parallelism with shared-memory threading. MPI is expected to have more communication overhead but may perform well on batches of independent alignments.

### 17.9 Hybrid Comparison

The current project does not implement MPI + OpenMP hybrid parallelism. Include this graph only as a future-work design unless a hybrid implementation is added.

Design:

```text
X-axis: Total workers
Y-axis: Runtime or speedup
Lines: MPI, OpenMP, MPI+OpenMP
```

Interpretation:

Hybrid parallelism can reduce the number of MPI processes per node while using OpenMP threads inside each process. It is expected to be useful on clusters with multiple cores per node.

## 18. Recommended Reporting Structure

Use the following structure in the final report:

```text
Experimental Setup
Correctness Validation
Execution Time Analysis
Strong Scaling Analysis
Weak Scaling Analysis
Speedup and Efficiency
Overhead Analysis
Threats to Validity
Conclusion
```

Mention threats to validity:

- VM scheduling noise
- Oversubscription for 8 MPI processes
- Limited number of repetitions
- Timing differences between implementations
- Memory pressure for very large matrices
- Random input variability

## 19. Final Benchmarking Checklist

Before collecting final results:

```text
[ ] Clean and rebuild all binaries.
[ ] Record system information.
[ ] Record compiler and MPI versions.
[ ] Generate inputs with a fixed seed.
[ ] Verify correctness against sequential baseline.
[ ] Run one warm-up per configuration.
[ ] Run at least 5 measured repetitions.
[ ] Prefer 10 measured repetitions for final reporting.
[ ] Save raw outputs.
[ ] Save raw repetition CSV.
[ ] Remove outliers using documented method.
[ ] Compute mean, median, standard deviation, min, and max.
[ ] Compute speedup and efficiency.
[ ] Mark oversubscribed MPI runs clearly.
[ ] Plot execution time vs input size.
[ ] Plot speedup vs workers.
[ ] Plot efficiency vs workers.
[ ] Plot strong scaling.
[ ] Plot weak scaling.
[ ] Plot Pthreads vs OpenMP.
[ ] Plot MPI vs OpenMP.
[ ] Discuss barrier overhead.
[ ] Discuss communication overhead.
[ ] Discuss load imbalance.
[ ] Discuss false sharing.
[ ] Discuss cache locality.
[ ] Discuss memory bandwidth.
[ ] Do not include failed or incorrect runs in final averages.
```

## 20. Expected Overall Interpretation

Small inputs are expected to show poor parallel performance because synchronization and startup overhead dominate the dynamic programming computation.

Medium inputs should begin to benefit from parallel execution, especially OpenMP and MPI task farming.

Large inputs should expose the true scalability limits of the implementations. Performance may eventually plateau because Needleman-Wunsch is memory intensive and the wavefront approach requires frequent synchronization.

MPI should be interpreted carefully. It is most meaningful for batches of independent sequence pairs. On a single VM, especially with oversubscription, MPI process counts above the available cores should be treated as functional experiments rather than true hardware scaling results.
