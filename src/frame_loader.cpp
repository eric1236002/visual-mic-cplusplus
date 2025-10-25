#include "frame_loader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <stdexcept>
#include <functional>

namespace fs = std::filesystem;

namespace visualmic {

VideoMetadata loadVideoMetadata(const std::string& frames_dir) {
    VideoMetadata metadata;
    
    std::string info_file = frames_dir + "/video_info.txt";
    std::ifstream file(info_file);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open video info file: " + info_file);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.find("fps=") == 0) {
            metadata.fps = std::stod(line.substr(4));
        } else if (line.find("frame_count=") == 0) {
            metadata.frame_count = std::stoi(line.substr(12));
        } else if (line.find("width=") == 0) {
            metadata.width = std::stoi(line.substr(6));
        } else if (line.find("height=") == 0) {
            metadata.height = std::stoi(line.substr(7));
        }
    }
    
    return metadata;
}

Matrix2D<double> loadPGMFrame(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open frame file: " + filename);
    }
    
    std::string magic;
    file >> magic;
    
    if (magic != "P5" && magic != "P2") {
        throw std::runtime_error("Invalid PGM format in: " + filename);
    }
    
    char c;
    file >> std::ws;
    while (file.peek() == '#') {
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        file >> std::ws;
    }
    
    int width, height, maxval;
    file >> width >> height >> maxval;
    file.ignore(1); 
    
    Matrix2D<double> mat(height, width);
    
    if (magic == "P5") {  
        std::vector<uint8_t> buffer(width * height);
        file.read(reinterpret_cast<char*>(buffer.data()), width * height);
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                mat.at(y, x) = static_cast<double>(buffer[y * width + x]);
            }
        }
    } else {  
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int pixel;
                file >> pixel;
                mat.at(y, x) = static_cast<double>(pixel);
            }
        }
    }
    
    return mat;
}

std::vector<std::string> getFrameFilesList(const std::string& frames_dir) {
    std::vector<std::string> frame_files;
    
    if (!fs::exists(frames_dir) || !fs::is_directory(frames_dir)) {
        throw std::runtime_error("Frames directory does not exist: " + frames_dir);
    }
    
    for (const auto& entry : fs::directory_iterator(frames_dir)) {
        if (entry.path().extension() == ".pgm") {
            frame_files.push_back(entry.path().string());
        }
    }
    
    if (frame_files.empty()) {
        throw std::runtime_error("No PGM frame files found in: " + frames_dir);
    }
    
    std::sort(frame_files.begin(), frame_files.end());
    
    return frame_files;
}

} // namespace visualmic

