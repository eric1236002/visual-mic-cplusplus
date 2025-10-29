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

Threads | Execution Time (s) | Speedup | Efficiency (%)
--------|-------------------|---------|---------
EOF

echo "testing single thread performance..."
export OMP_NUM_THREADS=1
echo -n "Threads 1: "
START_TIME=$(date +%s.%N)
"$EXECUTABLE" "$INPUT_DIR" "$@" -o "$TEMP_OUTPUT" > /dev/null 2>&1
END_TIME=$(date +%s.%N)
BASELINE_TIME=$(echo "$END_TIME - $START_TIME" | bc -l)
echo "${BASELINE_TIME}s"
echo "1        | $BASELINE_TIME    | 1.00x  | 100.0%" >> "$RESULTS_FILE"

for threads in "${THREADS[@]:1}"; do
    echo -n "Threads $threads: "
    export OMP_NUM_THREADS=$threads
    

    START_TIME=$(date +%s.%N)
    "$EXECUTABLE" "$INPUT_DIR" "$@" -o "$TEMP_OUTPUT" > /dev/null 2>&1
    END_TIME=$(date +%s.%N)
    

    EXECUTION_TIME=$(echo "$END_TIME - $START_TIME" | bc -l)
    
    SPEEDUP=$(echo "scale=2; $BASELINE_TIME / $EXECUTION_TIME" | bc -l)
    
    EFFICIENCY=$(echo "scale=1; $SPEEDUP / $threads * 100" | bc -l)
    
    echo "${EXECUTION_TIME}s (Speedup: ${SPEEDUP}x, Efficiency: ${EFFICIENCY}%)"
    
    printf "%-8d | %-11s | %-6s | %-6s%%\n" "$threads" "$EXECUTION_TIME" "${SPEEDUP}x" "$EFFICIENCY" >> "$RESULTS_FILE"
done

rm -f "$TEMP_OUTPUT"
rm -f "$TEMP_OUTPUT2"

cat "$RESULTS_FILE"

