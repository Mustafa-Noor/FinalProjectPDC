# Benchmarking Methodology for Parallel Needleman-Wunsch

## 1. Objective

This document defines a reproducible benchmarking methodology for evaluating four implementations of the Needleman-Wunsch global sequence alignment algorithm:

- Sequential baseline
- Pthreads wavefront implementation
- OpenMP wavefront implementation
- MPI master-worker implementation

The goal is to measure execution time, strong scaling, weak scaling, speedup, and efficiency in a way suitable for a university report or IEEE-style experimental evaluation section.

## 2. Implementations Under Test

| Implementation | Parallel Model | Purpose |
|---|---|---|
| Sequential | Single-threaded | Baseline for correctness and speedup |
| Pthreads | Shared-memory threads | Manual wavefront parallelization |
| OpenMP | Shared-memory compiler directives | Directive-based wavefront parallelization |
| MPI | Distributed processes | Master-worker task farming |

All implementations must use the same scoring scheme:

```text
Match:    +1
Mismatch: -1
Gap:      -1
```

All parallel results must be compared against the sequential baseline before performance numbers are accepted.

## 3. Input Sizes

Use the following sequence lengths:

```text
100
500
1000
5000
10000
50000
```

For each input size, generate random DNA sequences using the alphabet:

```text
A C G T
```

Use a fixed random seed for reproducibility.

## 4. Thread and Process Counts

Use the following worker counts:

```text
1
2
4
8
```

For Pthreads and OpenMP, these represent thread counts.

For MPI, these represent process counts.

## 5. Metrics

### 5.1 Execution Time

Execution time is the wall-clock time required to complete an alignment workload.

```text
T = execution time in seconds
```

Use wall-clock timers:

- Sequential: standard timing function used in the implementation
- Pthreads: `clock_gettime()` or equivalent wall-clock timing
- OpenMP: `omp_get_wtime()`
- MPI: `MPI_Wtime()`

### 5.2 Speedup

Speedup measures how much faster the parallel implementation is compared with the sequential baseline.

```text
Speedup(p) = T_sequential / T_parallel(p)
```

Where `p` is the number of threads or processes.

### 5.3 Efficiency

Efficiency measures how effectively the available workers are used.

```text
Efficiency(p) = Speedup(p) / p
```

Ideal efficiency is `1.0`. In practice, efficiency decreases as worker count increases because of synchronization, communication, memory bandwidth limits, and load imbalance.

### 5.4 Parallel Fraction

The parallel fraction is the portion of the program that can benefit from parallel execution.

A high parallel fraction usually leads to better scalability.

Sequential portions include:

- Input reading
- Matrix initialization
- Traceback
- Thread/process setup
- Synchronization overhead
- Communication overhead

## 6. Scaling Laws

### 6.1 Amdahl's Law

Amdahl's Law describes the maximum possible speedup for a fixed problem size.

```text
Speedup(p) = 1 / ((1 - f) + f / p)
```

Where:

```text
f = parallel fraction
p = number of workers
```

Amdahl's Law explains why speedup eventually stops improving even when more threads or processes are added.

### 6.2 Gustafson's Law

Gustafson's Law describes scaling when the problem size increases with the number of workers.

```text
ScaledSpeedup(p) = p - alpha(p - 1)
```

Where:

```text
alpha = serial fraction
p = number of workers
```

This is important for Needleman-Wunsch because larger matrices provide more parallel work.

## 7. Strong Scaling Experiment

### Purpose

Measure how runtime changes when the input size is fixed and the number of workers increases.

### Configuration

For each input size:

```text
100, 500, 1000, 5000, 10000, 50000
```

Run each implementation with:

```text
1, 2, 4, 8 workers
```

### Expected Outcome

Runtime should decrease as worker count increases, but speedup will usually be sublinear.

Small inputs may show poor scaling because overhead dominates computation.

Large inputs should scale better because the dynamic programming matrix contains more work.

### Record

```text
Algorithm
Input size
Threads
Processes
Execution time
Speedup
Efficiency
Correctness status
```

## 8. Weak Scaling Experiment

### Purpose

Measure whether runtime remains stable as both workload and worker count increase.

Needleman-Wunsch has `O(n^2)` work. To keep work per worker approximately constant:

```text
n_p = n_1 * sqrt(p)
```

Example:

| Workers | Suggested Input Size |
|---:|---:|
| 1 | 1000 |
| 2 | 1414 |
| 4 | 2000 |
| 8 | 2828 |

If the project must use only the required input sizes, use an approximate mapping:

| Workers | Required Input Size Approximation |
|---:|---:|
| 1 | 1000 |
| 2 | 1000 or 5000 |
| 4 | 5000 |
| 8 | 10000 |

### Expected Outcome

Ideal weak scaling keeps runtime nearly constant.

Increasing runtime indicates synchronization overhead, memory bandwidth limits, cache effects, or MPI communication overhead.

### Record

```text
Algorithm
Workers
Input size
Execution time
Weak scaling efficiency
```

Weak scaling efficiency:

```text
E_weak(p) = T_1 / T_p
```

## 9. Repetitions and Averaging

Run each configuration at least:

```text
5 repetitions minimum
10 repetitions recommended
```

For IEEE-style reporting, use 10 repetitions and report:

```text
Mean runtime
Standard deviation
Minimum runtime
Maximum runtime
```

### Averaging Method

Use the arithmetic mean after outlier removal.

```text
Mean = sum(valid runtimes) / number of valid runs
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

Reject runtimes outside this range.

If only 5 repetitions are used, an alternative is to remove the minimum and maximum values and average the remaining 3.

## 10. Overheads to Analyze

### 10.1 Barrier Overhead

Wavefront parallelization requires synchronization between anti-diagonals.

The next anti-diagonal cannot begin until the previous anti-diagonal is complete.

This creates barrier overhead, especially for small matrices.

### 10.2 Communication Overhead

MPI requires communication between the master and worker processes.

Communication overhead includes:

- Sending sequence pairs
- Receiving scores and alignments
- Worker startup cost
- Message latency
- Serialization of task distribution at rank 0

MPI overhead is most visible for small inputs.

### 10.3 Load Imbalance

Load imbalance can occur because:

- Early and late anti-diagonals contain fewer cells
- Different sequence pairs may have different lengths
- MPI workers may finish at different times

Dynamic task farming helps reduce MPI load imbalance.

### 10.4 False Sharing

False sharing occurs when multiple threads write different data that lies on the same cache line.

This can slow down Pthreads and OpenMP even when there is no correctness race.

### 10.5 Cache Locality

Sequential row-major filling has strong spatial locality.

Wavefront traversal may reduce locality because cells on the same anti-diagonal are not always contiguous in memory.

### 10.6 Memory Bandwidth

Large Needleman-Wunsch matrices can become memory-bandwidth-bound.

Adding more threads may stop improving performance when memory bandwidth is saturated.

## 11. CSV Templates

### 11.1 Raw Run CSV

```csv
RunID,Algorithm,InputSize,SeqA_Length,SeqB_Length,Threads,Processes,ExecutionTime,Score,Correct
1,Sequential,1000,1000,1000,1,1,0.000000,0,Yes
1,Pthreads,1000,1000,1000,4,1,0.000000,0,Yes
1,OpenMP,1000,1000,1000,4,1,0.000000,0,Yes
1,MPI,1000,1000,1000,1,4,0.000000,0,Yes
```

### 11.2 Summary CSV

```csv
Algorithm,InputSize,Threads,Processes,MeanTime,StdDev,MinTime,MaxTime,Speedup,Efficiency
Sequential,1000,1,1,0.000000,0.000000,0.000000,0.000000,1.000000,1.000000
Pthreads,1000,4,1,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
OpenMP,1000,4,1,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
MPI,1000,1,4,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
```

### 11.3 Strong Scaling CSV

```csv
Algorithm,InputSize,Workers,T1,Tp,Speedup,Efficiency
Pthreads,5000,1,0.000000,0.000000,1.000000,1.000000
Pthreads,5000,2,0.000000,0.000000,0.000000,0.000000
Pthreads,5000,4,0.000000,0.000000,0.000000,0.000000
Pthreads,5000,8,0.000000,0.000000,0.000000,0.000000
```

### 11.4 Weak Scaling CSV

```csv
Algorithm,Workers,InputSize,ExecutionTime,WeakScalingEfficiency
OpenMP,1,1000,0.000000,1.000000
OpenMP,2,1414,0.000000,0.000000
OpenMP,4,2000,0.000000,0.000000
OpenMP,8,2828,0.000000,0.000000
```

## 12. Graph Design

### 12.1 Execution Time vs Input Size

```text
X-axis: Input size
Y-axis: Execution time in seconds
Scale: Log scale recommended
Lines: Sequential, Pthreads, OpenMP, MPI
```

Interpretation:

This graph shows how runtime grows with problem size. Since Needleman-Wunsch is `O(n^2)`, execution time should grow rapidly as input size increases.

### 12.2 Speedup vs Threads

```text
X-axis: Threads or processes
Y-axis: Speedup
Include: Ideal speedup line y = x
```

Interpretation:

If the measured curve is close to the ideal line, scaling is good. If it flattens, overhead or sequential work is limiting performance.

### 12.3 Efficiency vs Threads

```text
X-axis: Threads or processes
Y-axis: Efficiency
```

Interpretation:

Efficiency near `1.0` is ideal. Decreasing efficiency means additional workers are contributing less useful work.

### 12.4 Strong Scaling

```text
X-axis: Workers
Y-axis: Execution time
Separate lines: Different input sizes
```

Interpretation:

A good strong scaling curve decreases as workers increase. If the line flattens, the implementation has reached a scalability limit.

### 12.5 Weak Scaling

```text
X-axis: Workers
Y-axis: Execution time
```

Interpretation:

An ideal weak scaling curve is flat. If runtime increases, overhead grows with worker count.

### 12.6 Execution Time Comparison

```text
Graph type: Grouped bar chart
X-axis: Input size
Y-axis: Execution time
Bars: Sequential, Pthreads, OpenMP, MPI
```

Interpretation:

This graph shows which implementation is fastest for each input size.

### 12.7 Pthreads vs OpenMP

```text
X-axis: Thread count
Y-axis: Runtime or speedup
Lines: Pthreads, OpenMP
```

Interpretation:

This graph compares manual threading against compiler-directed threading.

### 12.8 MPI vs OpenMP

```text
X-axis: Worker count
Y-axis: Runtime or speedup
Lines: MPI, OpenMP
```

Interpretation:

This graph compares distributed-process overhead against shared-memory overhead.

### 12.9 Hybrid Comparison

The current project does not yet implement MPI + OpenMP hybrid parallelism.

For future work, use:

```text
X-axis: Total workers
Y-axis: Runtime or speedup
Lines: MPI, OpenMP, MPI+OpenMP
```

Interpretation:

Hybrid parallelism should be useful on clusters where each node has multiple CPU cores.

## 13. IEEE-Style Reporting Requirements

Record the following system information:

```text
CPU model
Number of cores
Number of hardware threads
RAM size
Operating system
Compiler name and version
Compiler flags
MPI implementation and version
OpenMP version
Random seed
Number of repetitions
Outlier removal method
```

Also report:

```text
All implementations use the same scoring scheme.
Correctness was verified before benchmarking.
Speedup is computed relative to the sequential implementation.
Raw outputs and CSV files are preserved for reproducibility.
```

## 14. Final Benchmarking Checklist

Before collecting final results:

```text
[ ] Clean and rebuild all binaries.
[ ] Record compiler and system information.
[ ] Generate inputs using a fixed seed.
[ ] Verify correctness against the sequential baseline.
[ ] Run a warm-up execution before timed runs.
[ ] Run each configuration at least 5 times.
[ ] Prefer 10 repetitions for final results.
[ ] Save all raw outputs.
[ ] Save CSV summaries.
[ ] Remove outliers using a documented method.
[ ] Compute mean runtime and standard deviation.
[ ] Compute speedup and efficiency.
[ ] Plot execution time vs input size.
[ ] Plot speedup vs workers.
[ ] Plot efficiency vs workers.
[ ] Plot strong scaling.
[ ] Plot weak scaling.
[ ] Compare Pthreads vs OpenMP.
[ ] Compare MPI vs OpenMP.
[ ] Discuss overheads and scalability limits.
[ ] Do not report failed or incorrect runs as benchmark data.
```

## 15. Expected Interpretation

Small inputs are expected to show limited or negative speedup because overhead dominates computation.

Medium and large inputs should show better speedup because the dynamic programming matrix contains more work.

Pthreads and OpenMP should generally outperform MPI on a single shared-memory machine for one alignment.

MPI should become more useful when many independent sequence pairs are processed across multiple processes or nodes.

Very large inputs may become limited by memory capacity and memory bandwidth rather than CPU compute.
