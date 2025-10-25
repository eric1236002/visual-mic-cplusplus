#ifndef SOUND_FROM_VIDEO_H
#define SOUND_FROM_VIDEO_H

#include <vector>
#include <string>
#include "math_ops.h"

namespace visualmic {

struct VideoFrame {
    std::vector<unsigned char> data;  
    int width;
    int height;
    int channels;  
};

struct VideoInfo {
    int frame_count;
    double fps;
    int width;
    int height;
};


std::vector<double> soundFromVideoStreaming(const std::string& frames_dir,
                                            int nscale, 
                                            int norientation, 
                                            double downsample_factor = 1.0);

} // namespace visualmic

#endif // SOUND_FROM_VIDEO_H

