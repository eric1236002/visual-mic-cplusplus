#!/bin/bash

# Video Frame Extraction Script
# Usage: ./extract_frames.sh <input_video> <output_dir>

if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_video> <output_dir> <fps>"
    exit 1
fi

INPUT_VIDEO="$1"
OUTPUT_DIR="$2"
FPS="$3"
# Check if input file exists
if [ ! -f "$INPUT_VIDEO" ]; then
    echo "Error: Cannot find input video file: $INPUT_VIDEO"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Get video information
echo "Analyzing video..."
FRAME_COUNT=$(ffprobe -v error -select_streams v:0 -count_frames -show_entries stream=nb_read_frames -of csv=p=0 "$INPUT_VIDEO")
WIDTH=$(ffprobe -v error -select_streams v:0 -show_entries stream=width -of csv=p=0 "$INPUT_VIDEO")
HEIGHT=$(ffprobe -v error -select_streams v:0 -show_entries stream=height -of csv=p=0 "$INPUT_VIDEO")

echo "Video information:"
echo "  FPS: $FPS"
echo "  Frame count: $FRAME_COUNT"
echo "  Resolution: ${WIDTH}x${HEIGHT}"

# Save video info to file
cat > "$OUTPUT_DIR/video_info.txt" << EOF
fps=$FPS
frame_count=$FRAME_COUNT
width=$WIDTH
height=$HEIGHT
EOF

echo ""
echo "Extracting grayscale frames..."

# Use FFmpeg to extract grayscale frames as PGM format
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
    echo "Error: Frame extraction failed"
    exit 1
fi
