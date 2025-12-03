#ifndef SOUND_FROM_VIDEO_H
#define SOUND_FROM_VIDEO_H

#include <vector>
#include <string>
#include <map>
#include <queue>
#include <pthread.h>
#include "math_ops.h"
#include "steerable_pyramid.h"

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

struct ThreadData {
    const std::vector<std::string>* frame_files;
    int nscale;
    int norientation;
    double downsample_factor;
    const SteerablePyramidFreq* first_pyramid;
    std::map<BandKey, std::vector<double>> signals;
    pthread_t thread_id;
    
    // Timing data
    double total_load_time;
    double total_resize_time;
    double total_normalize_time;
    double total_pyramid_time;
    double total_bandproc_time;
    int frames_processed;

    // Shared data for producer-consumer
    std::queue<Matrix2D<double>>* frame_queue;
    pthread_mutex_t* queue_mutex;
    pthread_cond_t* queue_cond;
    bool* loading_complete;
    int start_frame; // For loader thread
    int end_frame;   // For loader thread
};


std::vector<double> soundFromVideoStreaming(const std::string& frames_dir,
                                            int nscale, 
                                            int norientation, 
                                            double downsample_factor = 1.0,
                                            int num_threads = 0,
                                            const std::string& parallel_mode = "all");

} // namespace visualmic

#endif // SOUND_FROM_VIDEO_H

