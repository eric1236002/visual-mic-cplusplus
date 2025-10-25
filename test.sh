# Usage: ./test.sh <input_video> <output_dir>

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

INPUT_VIDEO="vid_2025-10-22_13-47-33.mp4"
OUTPUT_DIR="data/${INPUT_VIDEO%.*}"
FPS=2200

cpp_INPUT_VIDEO="./data/${INPUT_VIDEO%.*}"
DOWNSAMPLE=0.1
NSCALE=1
NORIENT=2
BASENAME=${cpp_INPUT_VIDEO##*/}

sound_OUTPUT_DIR="output/${INPUT_VIDEO%.*}"


if [ ! -f "$INPUT_VIDEO" ]; then
    echo "Error: Cannot find input video file: $INPUT_VIDEO"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Analyzing video..."
FRAME_COUNT=$(ffprobe -v error -select_streams v:0 -count_frames -show_entries stream=nb_read_frames -of csv=p=0 "$INPUT_VIDEO")
WIDTH=$(ffprobe -v error -select_streams v:0 -show_entries stream=width -of csv=p=0 "$INPUT_VIDEO")
HEIGHT=$(ffprobe -v error -select_streams v:0 -show_entries stream=height -of csv=p=0 "$INPUT_VIDEO")

cat > "$OUTPUT_DIR/video_info.txt" << EOF
fps=$FPS
frame_count=$FRAME_COUNT
width=$WIDTH
height=$HEIGHT
EOF

echo ""
echo -e "${GREEN}Extracting grayscale frames...${NC}"

ffmpeg -i "$INPUT_VIDEO" \
    -vf "scale=iw*0.2:ih*0.2,format=gray" \
    -pix_fmt gray \
    "$OUTPUT_DIR/frame_%06d.pgm" \
    -loglevel warning -stats

if [ $? -eq 0 ]; then
    EXTRACTED_COUNT=$(ls "$OUTPUT_DIR"/frame_*.pgm 2>/dev/null | wc -l)
    echo ""
    echo "Done! Extracted $EXTRACTED_COUNT frames to: $OUTPUT_DIR"
    echo "Video info saved to: $OUTPUT_DIR/video_info.txt"
else
    echo -e "${RED}Error: Frame extraction failed${NC}"
    exit 1
fi



echo -e "${GREEN}Visual Microphone C++${NC}"
echo "======================================"
echo

if [ ! -f "./build/visual_microphone" ]; then
    echo -e "${RED}Error: visual_microphone executable not found!${NC}"
    echo "Please build the project first:"
    echo "  cd .."
    echo "  mkdir build && cd build"
    echo "  cmake .. && make"
    exit 1
fi

EXECUTABLE="time ./build/visual_microphone"

mkdir -p "$sound_OUTPUT_DIR"


echo -e "${YELLOW}Running with default parameters...${NC}"
echo -e "${YELLOW}Output directory: $sound_OUTPUT_DIR${NC}"
echo

$EXECUTABLE "$cpp_INPUT_VIDEO" \
    -o "$sound_OUTPUT_DIR/${BASENAME}_recovered.wav" \
    -s $FPS \
    -d $DOWNSAMPLE \
    -n $NSCALE \
    -r $NORIENT

if [ $? -eq 0 ]; then

    echo "Generated files:"
    echo -e "${YELLOW}  - $sound_OUTPUT_DIR/${BASENAME}_recovered.wav (original recovered sound)${NC}"
    echo -e "${YELLOW}  - $sound_OUTPUT_DIR/${BASENAME}_recovered_specsub.wav (enhanced with spectral subtraction)${NC}"
    echo "Play"
    echo "  ffplay $sound_OUTPUT_DIR/${BASENAME}_recovered.wav"
    echo "  ffplay $sound_OUTPUT_DIR/${BASENAME}_recovered_specsub.wav"
else
    echo
    echo -e "${RED}Error occurred during processing!${NC}"
    exit 1
fi
