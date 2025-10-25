#ifndef FRAME_LOADER_H
#define FRAME_LOADER_H

#include <vector>
#include <string>
#include <functional>
#include "math_ops.h"

namespace visualmic {

struct VideoMetadata {
    double fps;
    int frame_count;
    int width;
    int height;
};

VideoMetadata loadVideoMetadata(const std::string& frames_dir);

Matrix2D<double> loadPGMFrame(const std::string& filename);


std::vector<std::string> getFrameFilesList(const std::string& frames_dir);


#endif // FRAME_LOADER_H

}