#!/bin/bash

# ./benchmark_threads.sh <input_frames_dir> [Options...]

set -e 
if [ $# -lt 1 ]; then
    echo "usage: $0 <input_frames_dir> [Options...]"
    exit 1
fi

INPUT_DIR="$1"
shift  

if [ ! -d "$INPUT_DIR" ]; then
    echo "error: input directory '$INPUT_DIR' does not exist"
    exit 1
fi

EXECUTABLE="./build/visual_microphone"
MPI_EXEC=${MPI_EXEC:-mpirun}
if [ ! -f "$EXECUTABLE" ]; then
    echo "error: executable file '$EXECUTABLE' does not exist"
    exit 1
fi


THREADS=(1 2 4 8 16)
RESULTS_FILE="benchmark_results_$(date +%Y%m%d_%H%M%S).txt"
TEMP_OUTPUT="temp_output_$$.wav"
TEMP_OUTPUT2="temp_output_$$_specsub.wav"

echo "Input directory: $INPUT_DIR"
echo "Results file: $RESULTS_FILE"
echo ""


cat > "$RESULTS_FILE" << EOF
Thread performance test results
==================
Test time: $(date)
Input directory: $INPUT_DIR
Options: $@

Processes | Execution Time (s) | Speedup | Efficiency (%)
---------|---------------------|---------|---------
EOF

echo "testing single-process MPI performance..."
echo -n "Processes 1: "
START_TIME=$(date +%s.%N)
$MPI_EXEC -np 1 "$EXECUTABLE" "$INPUT_DIR" "$@" -o "$TEMP_OUTPUT" > /dev/null 2>&1
END_TIME=$(date +%s.%N)
BASELINE_TIME=$(echo "$END_TIME - $START_TIME" | bc -l)
echo "${BASELINE_TIME}s"
echo "1        | $BASELINE_TIME        | 1.00x  | 100.0%" >> "$RESULTS_FILE"

for procs in "${THREADS[@]:1}"; do
    echo -n "Processes $procs: "
    START_TIME=$(date +%s.%N)
    

    MPI_CMD="$MPI_EXEC -np $procs"
    if ! $MPI_CMD "$EXECUTABLE" "$INPUT_DIR" "$@" -o "$TEMP_OUTPUT" > /dev/null 2>&1; then

        MPI_CMD="$MPI_EXEC --oversubscribe -np $procs"
        if ! $MPI_CMD "$EXECUTABLE" "$INPUT_DIR" "$@" -o "$TEMP_OUTPUT" > /tmp/mpi_error_$$.log 2>&1; then
            END_TIME=$(date +%s.%N)
            echo "FAILED"
            echo "Error: Failed to run with $procs processes." >&2
            echo "MPI error output:" >&2
            cat /tmp/mpi_error_$$.log >&2
            rm -f /tmp/mpi_error_$$.log
            printf "%-8d | FAILED     | N/A   | N/A%%\n" "$procs" >> "$RESULTS_FILE"
            continue
        fi
    fi
    
    END_TIME=$(date +%s.%N)
    EXECUTION_TIME=$(echo "$END_TIME - $START_TIME" | bc -l)

    SPEEDUP=$(echo "scale=2; $BASELINE_TIME / $EXECUTION_TIME" | bc -l)

    EFFICIENCY=$(echo "scale=1; $SPEEDUP / $procs * 100" | bc -l)

    echo "${EXECUTION_TIME}s (Speedup: ${SPEEDUP}x, Efficiency: ${EFFICIENCY}%)"
    
    printf "%-8d | %-11s | %-6s | %-6s%%\n" "$procs" "$EXECUTION_TIME" "${SPEEDUP}x" "$EFFICIENCY" >> "$RESULTS_FILE"
done

rm -f "$TEMP_OUTPUT"
rm -f "$TEMP_OUTPUT2"

cat "$RESULTS_FILE"

