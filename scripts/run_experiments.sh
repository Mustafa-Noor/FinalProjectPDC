#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR" || exit 1

THREAD_COUNTS=(1 2 4 8)
INPUT_SIZES=(100 500 1000 5000 10000 50000)
PAIRS_PER_SIZE="${PAIRS_PER_SIZE:-3}"
SEED="${SEED:-2026}"
INPUT_DIR="data/generated"
RESULT_DIR="results"
RAW_DIR="$RESULT_DIR/raw_outputs"
CSV_FILE="$RESULT_DIR/benchmark_results.csv"
COMPILE_LOG="$RAW_DIR/compile.log"

mkdir -p bin "$INPUT_DIR" "$RAW_DIR"

echo "Benchmark started at $(date)" > "$COMPILE_LOG"
echo "Pairs per size: $PAIRS_PER_SIZE" >> "$COMPILE_LOG"
echo "Seed: $SEED" >> "$COMPILE_LOG"

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
    local raw_prefix="$RAW_DIR/${algorithm}_len_${size}_workers_${workers}"
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
    local raw_file="$RAW_DIR/mpi_len_${size}_processes_${processes}.txt"

    if command -v mpirun >/dev/null 2>&1; then
        if (( processes > 4 )); then
            mpirun --oversubscribe -np "$processes" ./bin/main_mpi < "$input_file" > "$raw_file" 2>&1
        else
            mpirun -np "$processes" ./bin/main_mpi < "$input_file" > "$raw_file" 2>&1
        fi
    elif command -v mpiexec >/dev/null 2>&1; then
        if (( processes > 4 )); then
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

echo "Algorithm,Input size,Threads,Processes,Execution time,Speedup,Efficiency" > "$CSV_FILE"

declare -A BASELINE_TIME

for size in "${INPUT_SIZES[@]}"; do
    input_file="$INPUT_DIR/input_len_${size}_pairs_${PAIRS_PER_SIZE}.txt"
    read_pairs_into_arrays "$input_file" "$PAIRS_PER_SIZE"

    if [[ "$SEQ_READY" -eq 1 ]]; then
        seq_time="$(run_pair_program_batch "sequential" "./bin/main_seq" "$size" 1 "sequential" "$PAIRS_PER_SIZE")"
        BASELINE_TIME[$size]="$seq_time"
        echo "Sequential,$size,1,1,$seq_time,1.000000,1.000000" >> "$CSV_FILE"
    else
        BASELINE_TIME[$size]="NA"
        echo "Sequential,$size,1,1,NA,NA,NA" >> "$CSV_FILE"
    fi

    for threads in "${THREAD_COUNTS[@]}"; do
        if [[ "$PTH_READY" -eq 1 ]]; then
            pth_time="$(run_pair_program_batch "pthreads" "./bin/main_pth" "$size" "$threads" "pthreads" "$PAIRS_PER_SIZE")"
            if [[ "${BASELINE_TIME[$size]}" != "NA" && "$pth_time" != "NA" ]]; then
                speedup="$(csv_ratio "${BASELINE_TIME[$size]}" "$pth_time")"
                efficiency="$(csv_ratio "$speedup" "$threads")"
            else
                speedup="NA"
                efficiency="NA"
            fi
            echo "Pthreads,$size,$threads,1,$pth_time,$speedup,$efficiency" >> "$CSV_FILE"
        fi

        if [[ "$OMP_READY" -eq 1 ]]; then
            omp_time="$(run_pair_program_batch "OpenMP" "./bin/main_omp" "$size" "$threads" "openmp" "$PAIRS_PER_SIZE")"
            if [[ "${BASELINE_TIME[$size]}" != "NA" && "$omp_time" != "NA" ]]; then
                speedup="$(csv_ratio "${BASELINE_TIME[$size]}" "$omp_time")"
                efficiency="$(csv_ratio "$speedup" "$threads")"
            else
                speedup="NA"
                efficiency="NA"
            fi
            echo "OpenMP,$size,$threads,1,$omp_time,$speedup,$efficiency" >> "$CSV_FILE"
        fi

        if [[ "$MPI_READY" -eq 1 ]]; then
            mpi_time="$(run_mpi_batch "$size" "$threads" "$input_file")"
            if [[ "${BASELINE_TIME[$size]}" != "NA" && "$mpi_time" != "NA" ]]; then
                speedup="$(csv_ratio "${BASELINE_TIME[$size]}" "$mpi_time")"
                efficiency="$(csv_ratio "$speedup" "$threads")"
            else
                speedup="NA"
                efficiency="NA"
            fi
            echo "MPI,$size,1,$threads,$mpi_time,$speedup,$efficiency" >> "$CSV_FILE"
        fi
    done
done

echo "Benchmark finished."
echo "CSV: $CSV_FILE"
echo "Raw outputs: $RAW_DIR"
