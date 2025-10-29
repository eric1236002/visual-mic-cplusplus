#include <iostream>
#include <fstream>
#include <string>
#include <cstdint>
#include <chrono>
#include "sound_from_video.h"
#include "frame_loader.h"
#include "signal_processing.h"
#include <omp.h>
void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " <input_frames_dir> [options]\n"
              << "\nNote: This version uses pre-extracted frames from FFmpeg.\n"
              << "      Use scripts/extract_frames.sh to extract frames first:\n"
              << "        ./scripts/extract_frames.sh video.mp4 frames_output\n"
              << "\nOptions:\n"
              << "  -o <output>         Output WAV file path (default: recoveredsound.wav)\n"
              << "  -s <sampling_rate>  Video frame rate (default: auto from video_info.txt)\n"
              << "  -d <downsample>     Downsample factor (default: 0.1)\n"
              << "  -n <nscale>         Number of pyramid scales (default: 1)\n"
              << "  -r <norient>        Number of orientations (default: 2)\n"
              << "  -h                  Show this help message\n"
              << std::endl;
}

void writeWavFile(const std::string& filename, const std::vector<double>& data, int sample_rate) {
    std::vector<int16_t> int_data(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        double val = std::max(-1.0, std::min(1.0, data[i]));
        int_data[i] = static_cast<int16_t>(val * 32767.0);
    }
    
    std::ofstream file(filename, std::ios::binary);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open output file: " + filename);
    }

    uint32_t chunk_size = 36 + int_data.size() * sizeof(int16_t);
    uint32_t subchunk2_size = int_data.size() * sizeof(int16_t);
    uint16_t audio_format = 1; 
    uint16_t num_channels = 1; 
    uint32_t byte_rate = sample_rate * num_channels * sizeof(int16_t);
    uint16_t block_align = num_channels * sizeof(int16_t);
    uint16_t bits_per_sample = 16;
    
    file.write("RIFF", 4);
    file.write(reinterpret_cast<char*>(&chunk_size), 4);
    file.write("WAVE", 4);
    file.write("fmt ", 4);
    
    uint32_t subchunk1_size = 16;
    file.write(reinterpret_cast<char*>(&subchunk1_size), 4);
    file.write(reinterpret_cast<char*>(&audio_format), 2);
    file.write(reinterpret_cast<char*>(&num_channels), 2);
    file.write(reinterpret_cast<char*>(&sample_rate), 4);
    file.write(reinterpret_cast<char*>(&byte_rate), 4);
    file.write(reinterpret_cast<char*>(&block_align), 2);
    file.write(reinterpret_cast<char*>(&bits_per_sample), 2);
    
    file.write("data", 4);
    file.write(reinterpret_cast<char*>(&subchunk2_size), 4);
    file.write(reinterpret_cast<char*>(int_data.data()), subchunk2_size);
    
    file.close();
    
    std::cout << "Saved WAV file: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }
    int num_threads = omp_get_max_threads();
    std::cout << "Using " << num_threads << " threads" << std::endl;
    std::string frames_dir = argv[1];
    std::string output_file = "recoveredsound.wav";
    int sampling_rate = -1;  
    double downsample_factor = 0.1;
    int nscale = 1;
    int norientation = 2;
    

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-o" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "-s" && i + 1 < argc) {
            sampling_rate = std::stoi(argv[++i]);
        } else if (arg == "-d" && i + 1 < argc) {
            downsample_factor = std::stod(argv[++i]);
        } else if (arg == "-n" && i + 1 < argc) {
            nscale = std::stoi(argv[++i]);
        } else if (arg == "-r" && i + 1 < argc) {
            norientation = std::stoi(argv[++i]);
        }
    }
    
    try {
        auto start_time = std::chrono::high_resolution_clock::now();
        std::cout << "Visual Microphone - Sound Recovery from Video Frames\n";
        std::cout << "====================================================\n" << std::endl;

        try {
            auto metadata = visualmic::loadVideoMetadata(frames_dir);
            std::cout << "Video metadata loaded:" << std::endl;
            std::cout << "  FPS: " << metadata.fps << std::endl;
            std::cout << "  Frame count: " << metadata.frame_count << std::endl;
            std::cout << "  Resolution: " << metadata.width << "x" << metadata.height << std::endl;
            
            if (sampling_rate < 0) {
                sampling_rate = static_cast<int>(std::round(metadata.fps));
                std::cout << "  Using FPS as sampling rate: " << sampling_rate << " Hz" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cout << "Note: Could not load video_info.txt: " << e.what() << std::endl;
            if (sampling_rate < 0) {
                sampling_rate = 30;
                std::cout << "Using default sampling rate: " << sampling_rate << " Hz" << std::endl;
            }
        }
        
        std::cout << "\nProcessing parameters:" << std::endl;
        std::cout << "  Sampling rate: " << sampling_rate << " Hz" << std::endl;
        std::cout << "  Downsample factor: " << downsample_factor << std::endl;
        std::cout << "  Pyramid scales: " << nscale << std::endl;
        std::cout << "  Orientations: " << norientation << std::endl;
        

        auto extract_start = std::chrono::high_resolution_clock::now();
        auto sound = visualmic::soundFromVideoStreaming(frames_dir, nscale, norientation, downsample_factor);
        auto extract_end = std::chrono::high_resolution_clock::now();
        auto extract_time = std::chrono::duration_cast<std::chrono::duration<double>>(extract_end - extract_start);
        
        std::cout << "\nExtraction time: " << extract_time.count() << " seconds" << std::endl;

        std::cout << "\nSaving recovered sound to: " << output_file << std::endl;
        writeWavFile(output_file, sound, sampling_rate);
        
        auto save_end = std::chrono::high_resolution_clock::now();
        auto save_time = std::chrono::duration_cast<std::chrono::milliseconds>(save_end - extract_end);
        std::cout << "Save time: " << save_time.count() << " ms" << std::endl;
        std::cout << "\nApplying spectral subtraction..." << std::endl;
        auto sound_specsub = visualmic::getSoundSpectralSubtraction(sound);

        size_t dot_pos = output_file.find_last_of('.');
        std::string base_name = output_file.substr(0, dot_pos);
        std::string extension = output_file.substr(dot_pos);
        std::string specsub_file = base_name + "_specsub" + extension;
        
        std::cout << "Saving spectral subtraction enhanced sound..." << std::endl;
        writeWavFile(specsub_file, sound_specsub, sampling_rate);
        
        std::cout << "\n" << std::string(50, '=') << std::endl;
        std::cout << "Processing complete!" << std::endl;
        std::cout << std::string(50, '=') << std::endl;
        
        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::duration<double>>(total_end - start_time);
        std::cout << "\nTotal time: " << total_time.count() << " seconds" << std::endl;
        std::cout << "Output file: " << output_file << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

