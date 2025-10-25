#!/bin/bash

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}Visual Microphone C++ - Demo Script${NC}"
echo "======================================"
echo

# Check if executable exists
if [ ! -f "./build/visual_microphone" ]; then
    echo -e "${RED}Error: visual_microphone executable not found!${NC}"
    echo "Please build the project first:"
    echo "  cd .."
    echo "  mkdir build && cd build"
    echo "  cmake .. && make"
    exit 1
fi

EXECUTABLE="time ./build/visual_microphone"

INPUT_VIDEO="./data/Chips2-2200Hz-Mary_MIDI-input"
FPS=2200
DOWNSAMPLE=0.1
NSCALE=1
NORIENT=2
# Extract base name for output files
BASENAME=video

# Default parameters
OUTPUT_DIR="./output"
mkdir -p "$OUTPUT_DIR"

# Run with default parameters
echo -e "${YELLOW}Running with default parameters...${NC}"
echo "Output directory: $OUTPUT_DIR"
echo

$EXECUTABLE "$INPUT_VIDEO" \
    -o "$OUTPUT_DIR/${BASENAME}_recovered.wav" \
    -s $FPS \
    -d $DOWNSAMPLE \
    -n $NSCALE \
    -r $NORIENT

if [ $? -eq 0 ]; then
    echo
    echo -e "${GREEN}Success!${NC}"
    echo
    echo "Generated files:"
    echo -e "${YELLOW}  - $OUTPUT_DIR/${BASENAME}_recovered.wav (original recovered sound)${NC}"
    echo -e "${YELLOW}  - $OUTPUT_DIR/${BASENAME}_recovered_specsub.wav (enhanced with spectral subtraction)${NC}"
    echo
    echo "You can play the audio files with:"
    echo "  ffplay $OUTPUT_DIR/${BASENAME}_recovered.wav"
    echo "  ffplay $OUTPUT_DIR/${BASENAME}_recovered_specsub.wav"
else
    echo
    echo -e "${RED}Error occurred during processing!${NC}"
    exit 1
fi
