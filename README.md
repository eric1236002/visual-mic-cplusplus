# The Visual Microphone: Passive Recovery of Sound from Video (C++ Implementation)

A C++ implementation of MIT's Visual Microphone [1], which passively recovers sound from video by analyzing subtle vibrations visible in the video frames.

This implementation is based on the original [MATLAB version](http://people.csail.mit.edu/mrub/VisualMic/#data) and the [Python version](https://github.com/antoniomuso/visual-mic). It uses FFmpeg for video frame extraction and custom implementations of steerable pyramids and signal processing algorithms.

## Key Features

- **Memory Efficient**: Streaming processing for large videos (no memory limitations)
- **Parallel Processing**: Full OpenMP parallelization for multi-core CPUs
- **No External Dependencies**: Pure C++ implementation with minimal dependencies
- **High Performance**: Optimized algorithms for real-time processing
- **Cross Platform**: Works on macOS, Linux, and Windows



## Requirements

- CMake 3.15 or higher
- C++17 compatible compiler with OpenMP support (GCC 7+, Clang 5+, MSVC 2017+)
- FFmpeg (for video frame extraction)
- OpenMP library (usually included with compiler)

### Installing Dependencies

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install cmake g++ ffmpeg libomp-dev
```

#### macOS (with Homebrew)
```bash
brew install cmake ffmpeg libomp
```

#### Windows
- Install FFmpeg from [ffmpeg.org](https://ffmpeg.org/download.html)

## Building

```bash
mkdir build
cd build
cmake ..
make -j$(nproc)
```

This will create the `visual_microphone` executable.
## data

[Chips2-2200Hz-Mary_MIDI-input.avi](https://data.csail.mit.edu/vidmag/VisualMic/Results/Chips2-2200Hz-Mary_MIDI-input.avi)

## Usage

### Step 1: Extract Video Frames

First, extract video frames using the provided script:

```bash
./scripts/extract_frames.sh <input_video> <output_directory>
```

**Example:**
```bash
./scripts/extract_frames.sh Chips2-2200Hz-Mary_MIDI-input.avi ./data/Chips2-2200Hz-Mary_MIDI-input 2200
```

This will:
- Extract all frames as grayscale PGM images
- Save video metadata (FPS, resolution, frame count) to `video_info.txt`
- Create frames like `frame_000001.pgm`, `frame_000002.pgm`, etc.

### Step 2: Process Frames

```bash
./build/visual_microphone <frames_directory> <output_file>
```
#Options:
- `-o <output>` : Output WAV file path (default: `recoveredsound.wav`)
- `-s <sampling_rate>` : Video frame rate (default: auto from video_info.txt)
- `-d <downsample>` : Downsample factor for processing (default: 0.1)
- `-n <nscale>` : Number of pyramid scales (default: 1)
- `-r <norient>` : Number of orientations (default: 2)
- `-h` : Show help message

**Example:**
```bash
./build/visual_microphone ./data/Chips2-2200Hz-Mary_MIDI-input -o ./output/Chips2-2200Hz-Mary_MIDI-input_recovered.wav -s 2200 -d 0.1 -n 1 -r 2
```

This will generate:
- `output.wav` - The recovered sound from video vibrations
- `output_specsub.wav` - Enhanced version with spectral subtraction

### Controlling Thread Count

The program automatically uses all available CPU cores. You can control the number of threads:

```bash
# Use 4 threads
export OMP_NUM_THREADS=4
./build/visual_microphone <frames_directory> -o <output>

# Use all available cores (default)
unset OMP_NUM_THREADS
./build/visual_microphone <frames_directory> -o <output>
```


## How It Works

The Visual Microphone works by:

1. **Frame Extraction**: Video frames are extracted as grayscale images using FFmpeg.

2. **Streaming Processing**: Frames are processed one by one to minimize memory usage (supports videos of any length).

3. **Steerable Pyramid Decomposition**: Each frame is decomposed using a complex steerable pyramid in the frequency domain, creating multiple scales and orientations.

4. **Phase Analysis**: For each pyramid band, the phase difference between the current frame and the first frame is computed (Formula 2 in the paper).

5. **Motion Signal Extraction**: A single motion signal is computed by multiplying the phase difference with the squared amplitude (Formula 3).

6. **Signal Alignment**: Signals from different bands are aligned using cross-correlation (Formula 4).

7. **Signal Combination**: All aligned signals are summed to produce the final sound signal (Formula 5).

8. **Post-processing**: 
   - Highpass filtering to remove low-frequency noise
   - Normalization to [-1, 1] range

## Implementation Details

### Custom Implementations

This implementation includes from-scratch versions of:

- **Steerable Pyramid**: Multi-scale and multi-orientation decomposition in frequency domain
- **FFT/IFFT**: Fast Fourier Transform and its inverse for frequency domain processing
- **Butterworth Filter**: Second-order sections (SOS) implementation for highpass filtering
- **Cross-Correlation**: For temporal alignment of signals
- **Frame Loading**: Efficient PGM image format reader

### Design Constraints

- **No External Libraries**: Pure C++ implementation with minimal dependencies
- **Memory Efficient**: Streaming processing prevents memory overflow on large videos
- **Cross Platform**: No platform-specific dependencies
- **Self-Contained**: All algorithms implemented from scratch

## Project Structure

```
visual-mic-cpp/
├── CMakeLists.txt
├── README.md
├── scripts/
│   └── extract_frames.sh          # Video frame extraction script
├── include/
│   ├── frame_loader.h             # Frame loading utilities
│   ├── math_ops.h                 # Mathematical operations and FFT
│   ├── signal_processing.h        # Butterworth filter implementation
│   ├── sound_from_video.h         # Main sound extraction interface
│   ├── steerable_pyramid.h        # Steerable pyramid decomposition
│   └── utils.h                    # Utility functions
└── src/
    ├── frame_loader.cpp           # PGM frame loading implementation
    ├── main.cpp                   # Main application entry point
    ├── math_ops.cpp               # FFT and matrix operations
    ├── signal_processing.cpp      # Filter implementation
    ├── sound_from_video.cpp       # Core sound extraction algorithm
    ├── steerable_pyramid.cpp      # Pyramid decomposition
    └── utils.cpp                  # Utility functions
```

## References

[1] DAVIS, Abe, et al. The visual microphone: Passive recovery of sound from video. ACM Transactions on Graphics (TOG), 2014.

[2] PORTILLA, Javier; SIMONCELLI, Eero P. A parametric texture model based on joint statistics of complex wavelet coefficients. International Journal of Computer Vision, 2000.

## License

This implementation follows the same licensing as the original work. Please refer to the [original project](http://people.csail.mit.edu/mrub/VisualMic/) for licensing details.

