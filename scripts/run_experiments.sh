#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR" || exit 1

THREAD_COUNTS=(1 2 4 8)
INPUT_SIZES=(100 500 1000 5000 10000 50000)

if [[ -n "${THREAD_COUNTS_OVERRIDE:-}" ]]; then
    read -r -a THREAD_COUNTS <<< "$THREAD_COUNTS_OVERRIDE"
fi

if [[ -n "${INPUT_SIZES_OVERRIDE:-}" ]]; then
    read -r -a INPUT_SIZES <<< "$INPUT_SIZES_OVERRIDE"
fi

PAIRS_PER_SIZE="${PAIRS_PER_SIZE:-3}"
SEED="${SEED:-2026}"
REPETITIONS="${REPETITIONS:-10}"
WARMUP_RUNS="${WARMUP_RUNS:-1}"
MPI_OVERSUBSCRIBE_ABOVE="${MPI_OVERSUBSCRIBE_ABOVE:-4}"

INPUT_DIR="data/generated"
RESULT_DIR="results"
RAW_DIR="$RESULT_DIR/raw_outputs"
SUMMARY_CSV="$RESULT_DIR/benchmark_results.csv"
RAW_CSV="$RESULT_DIR/benchmark_raw_runs.csv"
CONFIG_CSV="$RESULT_DIR/benchmark_configurations.csv"
COMPILE_LOG="$RAW_DIR/compile.log"

mkdir -p bin "$INPUT_DIR" "$RAW_DIR"

echo "Benchmark started at $(date)" > "$COMPILE_LOG"
echo "Pairs per size: $PAIRS_PER_SIZE" >> "$COMPILE_LOG"
echo "Seed: $SEED" >> "$COMPILE_LOG"
echo "Repetitions: $REPETITIONS" >> "$COMPILE_LOG"
echo "Warm-up runs: $WARMUP_RUNS" >> "$COMPILE_LOG"
echo "MPI oversubscribe above: $MPI_OVERSUBSCRIBE_ABOVE" >> "$COMPILE_LOG"

compile_target() {
    local target="$1"

    echo "Compiling $target..."
    if make "$target" >> "$COMPILE_LOG" 2>&1; then
        echo "  ok"
        return 0
    fi

    echo "  skipped: see $COMPILE_LOG"
    return 1
}

extract_time() {
    local file="$1"

    awk '
        /Execution time:/ { value=$3 }
        /Total execution time:/ { value=$4 }
        END {
            if (value == "") {
                print "NA"
            } else {
                print value
            }
        }
    ' "$file"
}

sum_times() {
    awk '
        BEGIN { total=0; ok=1 }
        $1 == "NA" { ok=0 }
        $1 != "NA" { total += $1 }
        END {
            if (ok == 0) {
                print "NA"
            } else {
                printf "%.6f\n", total
            }
        }
    '
}

csv_ratio() {
    local numerator="$1"
    local denominator="$2"

    awk -v numerator="$numerator" -v denominator="$denominator" 'BEGIN {
        if (numerator == "NA" || denominator == "NA" || denominator <= 0) {
            print "NA";
        } else {
            printf "%.6f", numerator / denominator;
        }
    }'
}

stats_csv() {
    local values="$1"

    printf "%s\n" "$values" | awk '
        function sort_numbers(a, n, i, j, tmp) {
            for (i = 1; i <= n; i++) {
                for (j = i + 1; j <= n; j++) {
                    if (a[j] < a[i]) {
                        tmp = a[i]; a[i] = a[j]; a[j] = tmp;
                    }
                }
            }
        }
        function median(a, lo, hi, count, mid) {
            count = hi - lo + 1;
            if (count <= 0) {
                return "NA";
            }
            mid = int(count / 2);
            if (count % 2 == 1) {
                return a[lo + mid];
            }
            return (a[lo + mid - 1] + a[lo + mid]) / 2.0;
        }
        $1 != "NA" && $1 != "" {
            n++;
            x[n] = $1 + 0.0;
        }
        END {
            if (n == 0) {
                print "NA,NA,NA,NA,NA,0,0";
                exit;
            }

            sort_numbers(x, n);
            q1 = median(x, 1, int((n + 1) / 2));
            q3_start = int((n + 2) / 2);
            q3 = median(x, q3_start, n);
            iqr = q3 - q1;
            low = q1 - 1.5 * iqr;
            high = q3 + 1.5 * iqr;

            for (i = 1; i <= n; i++) {
                if (x[i] >= low && x[i] <= high) {
                    m++;
                    y[m] = x[i];
                    sum += x[i];
                } else {
                    outliers++;
                }
            }

            if (m == 0) {
                print "NA,NA,NA,NA,NA," outliers ",0";
                exit;
            }

            mean = sum / m;
            for (i = 1; i <= m; i++) {
                diff = y[i] - mean;
                sq += diff * diff;
            }
            stddev = (m > 1) ? sqrt(sq / (m - 1)) : 0.0;
            med = median(y, 1, m);
            min = y[1];
            max = y[m];

            printf "%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d\n",
                   mean, stddev, med, min, max, outliers, m;
        }
    '
}

read_pairs_into_arrays() {
    local input_file="$1"
    local expected_pairs="$2"
    local index

    mapfile -t INPUT_LINES < "$input_file"
    SEQ_A=()
    SEQ_B=()

    for ((index = 0; index < expected_pairs; index++)); do
        SEQ_A[index]="${INPUT_LINES[$((1 + index * 2))]}"
        SEQ_B[index]="${INPUT_LINES[$((2 + index * 2))]}"
    done
}

run_pair_program_batch() {
    local algorithm="$1"
    local binary="$2"
    local size="$3"
    local workers="$4"
    local extra_mode="$5"
    local pair_count="$6"
    local run_id="$7"
    local raw_prefix="$RAW_DIR/${algorithm}_len_${size}_workers_${workers}_run_${run_id}"
    local times_file="$raw_prefix.times"
    local pair
    local raw_file
    local time_value

    : > "$times_file"

    for ((pair = 0; pair < pair_count; pair++)); do
        raw_file="${raw_prefix}_pair_${pair}.txt"

        if [[ "$extra_mode" == "pthreads" ]]; then
            printf "%s\n%s\n" "${SEQ_A[$pair]}" "${SEQ_B[$pair]}" \
                | "$binary" --threads "$workers" > "$raw_file" 2>&1
        elif [[ "$extra_mode" == "openmp" ]]; then
            printf "%s\n%s\n" "${SEQ_A[$pair]}" "${SEQ_B[$pair]}" \
                | OMP_NUM_THREADS="$workers" "$binary" > "$raw_file" 2>&1
        else
            printf "%s\n%s\n" "${SEQ_A[$pair]}" "${SEQ_B[$pair]}" \
                | "$binary" > "$raw_file" 2>&1
        fi

        if [[ $? -ne 0 ]]; then
            echo "NA" >> "$times_file"
        else
            time_value="$(extract_time "$raw_file")"
            echo "$time_value" >> "$times_file"
        fi
    done

    sum_times < "$times_file"
}

run_mpi_batch() {
    local size="$1"
    local processes="$2"
    local input_file="$3"
    local run_id="$4"
    local raw_file="$RAW_DIR/mpi_len_${size}_processes_${processes}_run_${run_id}.txt"

    if command -v mpirun >/dev/null 2>&1; then
        if (( processes > MPI_OVERSUBSCRIBE_ABOVE )); then
            mpirun --oversubscribe -np "$processes" ./bin/main_mpi < "$input_file" > "$raw_file" 2>&1
        else
            mpirun -np "$processes" ./bin/main_mpi < "$input_file" > "$raw_file" 2>&1
        fi
    elif command -v mpiexec >/dev/null 2>&1; then
        if (( processes > MPI_OVERSUBSCRIBE_ABOVE )); then
            mpiexec --oversubscribe -n "$processes" ./bin/main_mpi < "$input_file" > "$raw_file" 2>&1
        else
            mpiexec -n "$processes" ./bin/main_mpi < "$input_file" > "$raw_file" 2>&1
        fi
    else
        echo "MPI launcher not found." > "$raw_file"
        echo "NA"
        return
    fi

    if [[ $? -ne 0 ]]; then
        echo "NA"
    else
        extract_time "$raw_file"
    fi
}

append_raw_run() {
    local run_id="$1"
    local algorithm="$2"
    local size="$3"
    local threads="$4"
    local processes="$5"
    local time_value="$6"
    local status="$7"
    local notes="$8"

    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$run_id" "$algorithm" "$size" "$threads" "$processes" \
        "$PAIRS_PER_SIZE" "$time_value" "$status" "$notes" >> "$RAW_CSV"
}

run_configuration() {
    local algorithm="$1"
    local size="$2"
    local threads="$3"
    local processes="$4"
    local mode="$5"
    local binary="$6"
    local input_file="$7"
    local notes="$8"
    local run
    local warmup
    local time_value
    local status
    local values=""

    echo "Running $algorithm size=$size threads=$threads processes=$processes" >&2

    for ((warmup = 1; warmup <= WARMUP_RUNS; warmup++)); do
        if [[ "$mode" == "mpi" ]]; then
            run_mpi_batch "$size" "$processes" "$input_file" "warmup_${warmup}" >/dev/null
        else
            run_pair_program_batch "$algorithm" "$binary" "$size" "$threads" "$mode" "$PAIRS_PER_SIZE" "warmup_${warmup}" >/dev/null
        fi
    done

    for ((run = 1; run <= REPETITIONS; run++)); do
        if [[ "$mode" == "mpi" ]]; then
            time_value="$(run_mpi_batch "$size" "$processes" "$input_file" "$run")"
        else
            time_value="$(run_pair_program_batch "$algorithm" "$binary" "$size" "$threads" "$mode" "$PAIRS_PER_SIZE" "$run")"
        fi

        if [[ "$time_value" == "NA" ]]; then
            status="FAIL"
        else
            status="OK"
        fi

        append_raw_run "$run" "$algorithm" "$size" "$threads" "$processes" "$time_value" "$status" "$notes"
        values="${values}${time_value}"$'\n'
    done

    echo "$values"
}

summarize_configuration() {
    local algorithm="$1"
    local size="$2"
    local threads="$3"
    local processes="$4"
    local values="$5"
    local baseline_mean="$6"
    local notes="$7"
    local stats
    local mean_time
    local stddev
    local median_time
    local min_time
    local max_time
    local outliers
    local valid_runs
    local speedup
    local efficiency
    local workers

    stats="$(stats_csv "$values")"
    IFS=',' read -r mean_time stddev median_time min_time max_time outliers valid_runs <<< "$stats"

    if [[ "$algorithm" == "Sequential" ]]; then
        speedup="1.000000"
        efficiency="1.000000"
    else
        speedup="$(csv_ratio "$baseline_mean" "$mean_time")"
        if [[ "$processes" -gt 1 ]]; then
            workers="$processes"
        else
            workers="$threads"
        fi
        efficiency="$(csv_ratio "$speedup" "$workers")"
    fi

    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$algorithm" "$size" "$threads" "$processes" "$mean_time" "$stddev" \
        "$median_time" "$min_time" "$max_time" "$outliers" "$valid_runs" \
        "$speedup" "$efficiency" "$notes" >> "$SUMMARY_CSV"

    printf "%s\n" "$mean_time"
}

SEQ_READY=0
PTH_READY=0
OMP_READY=0
MPI_READY=0
GEN_READY=0

compile_target input_generator && GEN_READY=1
compile_target sequential && SEQ_READY=1
compile_target pthreads && PTH_READY=1
compile_target openmp && OMP_READY=1
compile_target mpi && MPI_READY=1

if [[ "$GEN_READY" -ne 1 ]]; then
    echo "Input generator failed to compile. Cannot continue."
    exit 1
fi

./bin/input_generator "$INPUT_DIR" "$PAIRS_PER_SIZE" "$SEED" >> "$COMPILE_LOG" 2>&1
if [[ $? -ne 0 ]]; then
    echo "Input generation failed. See $COMPILE_LOG"
    exit 1
fi

echo "RunID,Algorithm,InputSize,Threads,Processes,PairsPerSize,ExecutionTime,Status,Notes" > "$RAW_CSV"
echo "Algorithm,InputSize,Threads,Processes,MeanTime,StdDev,MedianTime,MinTime,MaxTime,OutliersRemoved,ValidRuns,Speedup,Efficiency,Notes" > "$SUMMARY_CSV"
echo "Algorithm,InputSize,Threads,Processes,ExpectedRuns" > "$CONFIG_CSV"

declare -A BASELINE_MEAN
expected_configs=0
completed_configs=0

for size in "${INPUT_SIZES[@]}"; do
    input_file="$INPUT_DIR/input_len_${size}_pairs_${PAIRS_PER_SIZE}.txt"
    read_pairs_into_arrays "$input_file" "$PAIRS_PER_SIZE"

    if [[ "$SEQ_READY" -eq 1 ]]; then
        expected_configs=$((expected_configs + 1))
        echo "Sequential,$size,1,1,$REPETITIONS" >> "$CONFIG_CSV"
        seq_values="$(run_configuration "Sequential" "$size" 1 1 "sequential" "./bin/main_seq" "$input_file" "")"
        BASELINE_MEAN[$size]="$(summarize_configuration "Sequential" "$size" 1 1 "$seq_values" "NA" "")"
        completed_configs=$((completed_configs + 1))
    else
        BASELINE_MEAN[$size]="NA"
        echo "Sequential,$size,1,1,NA,NA,NA,NA,NA,NA,0,NA,NA,not_compiled" >> "$SUMMARY_CSV"
    fi

    for workers in "${THREAD_COUNTS[@]}"; do
        if [[ "$PTH_READY" -eq 1 ]]; then
            expected_configs=$((expected_configs + 1))
            echo "Pthreads,$size,$workers,1,$REPETITIONS" >> "$CONFIG_CSV"
            pth_values="$(run_configuration "Pthreads" "$size" "$workers" 1 "pthreads" "./bin/main_pth" "$input_file" "")"
            summarize_configuration "Pthreads" "$size" "$workers" 1 "$pth_values" "${BASELINE_MEAN[$size]}" "" >/dev/null
            completed_configs=$((completed_configs + 1))
        fi

        if [[ "$OMP_READY" -eq 1 ]]; then
            expected_configs=$((expected_configs + 1))
            echo "OpenMP,$size,$workers,1,$REPETITIONS" >> "$CONFIG_CSV"
            omp_values="$(run_configuration "OpenMP" "$size" "$workers" 1 "openmp" "./bin/main_omp" "$input_file" "")"
            summarize_configuration "OpenMP" "$size" "$workers" 1 "$omp_values" "${BASELINE_MEAN[$size]}" "" >/dev/null
            completed_configs=$((completed_configs + 1))
        fi

        if [[ "$MPI_READY" -eq 1 ]]; then
            mpi_notes=""
            if (( workers > MPI_OVERSUBSCRIBE_ABOVE )); then
                mpi_notes="oversubscribed"
            fi

            expected_configs=$((expected_configs + 1))
            echo "MPI,$size,1,$workers,$REPETITIONS" >> "$CONFIG_CSV"
            mpi_values="$(run_configuration "MPI" "$size" 1 "$workers" "mpi" "./bin/main_mpi" "$input_file" "$mpi_notes")"
            summarize_configuration "MPI" "$size" 1 "$workers" "$mpi_values" "${BASELINE_MEAN[$size]}" "$mpi_notes" >/dev/null
            completed_configs=$((completed_configs + 1))
        fi
    done
done

echo "Benchmark finished."
echo "Expected configurations: $expected_configs"
echo "Completed configurations: $completed_configs"
echo "Summary CSV: $SUMMARY_CSV"
echo "Raw repetition CSV: $RAW_CSV"
echo "Configuration grid CSV: $CONFIG_CSV"
echo "Raw outputs: $RAW_DIR"
echo
echo "Check for failed runs with:"
echo "  awk -F, '\$8 != \"Status\" && \$8 != \"OK\" { print }' $RAW_CSV"
echo
echo "Check for incomplete summary rows with:"
echo "  awk -F, '\$11 == 0 || \$5 == \"NA\" { print }' $SUMMARY_CSV"
