#include "sound_from_video.h"
#include "steerable_pyramid.h"
#include "signal_processing.h"
#include "utils.h"
#include "math_ops.h"
#include "frame_loader.h"
#include <iostream>
#include <map>
#include <cmath>
#include <chrono>

namespace visualmic {

Matrix2D<double> resizeImage(const Matrix2D<double>& img, double scale_factor) {
    int new_rows = static_cast<int>(img.rows * scale_factor);
    int new_cols = static_cast<int>(img.cols * scale_factor);
    return resizeMatrix(img, new_rows, new_cols);
}


std::vector<double> soundFromVideoStreaming(const std::string& frames_dir,
                                            int nscale, 
                                            int norientation, 
                                            double downsample_factor) {
    
    std::vector<std::string> frame_files = getFrameFilesList(frames_dir);
    
    if (frame_files.empty()) {
        throw std::runtime_error("No frames found in directory");
    }
    
    int nframes = frame_files.size();
    std::cout << "Found " << nframes << " frames to process (streaming mode)" << std::endl;
    
    auto init_start = std::chrono::high_resolution_clock::now();
    std::cout << "Loading first frame..." << std::endl;
    Matrix2D<double> gray_frame = loadPGMFrame(frame_files[0]);
    
    if (downsample_factor < 1.0) {
        gray_frame = resizeImage(gray_frame, downsample_factor);
    }
    
    Matrix2D<double> norm_frame = normalizeMatrix(gray_frame);
    
    SteerablePyramidFreq first_pyramid(norm_frame, nscale, norientation - 1);
    auto first_pyramid_coeffs = first_pyramid.getPyrCoeffs();
    auto init_end = std::chrono::high_resolution_clock::now();
    auto init_time = std::chrono::duration_cast<std::chrono::duration<double>>(init_end - init_start);
    
    std::map<BandKey, std::vector<double>> signals;
    for (const auto& pair : first_pyramid_coeffs) {
        signals[pair.first] = std::vector<double>();
        signals[pair.first].reserve(nframes); 
    }
    
    std::cout << "Processing frames (streaming mode - low memory usage)..." << std::endl;
    
    int frame_count = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = std::chrono::duration<double>::zero();
    double acc_load_s = 0.0;
    double acc_resize_s = 0.0;
    double acc_normalize_s = 0.0;
    double acc_pyramid_s = 0.0;
    double acc_bandproc_s = 0.0;
    
    for (const auto& frame_file : frame_files) {
        auto t_load_start = std::chrono::high_resolution_clock::now();
        gray_frame = loadPGMFrame(frame_file);
        auto t_load_end = std::chrono::high_resolution_clock::now();
        acc_load_s += std::chrono::duration_cast<std::chrono::duration<double>>(t_load_end - t_load_start).count();
        
        if (downsample_factor < 1.0) {
            auto t_resize_start = std::chrono::high_resolution_clock::now();
            gray_frame = resizeImage(gray_frame, downsample_factor);
            auto t_resize_end = std::chrono::high_resolution_clock::now();
            acc_resize_s += std::chrono::duration_cast<std::chrono::duration<double>>(t_resize_end - t_resize_start).count();
        }
        
        auto t_norm_start = std::chrono::high_resolution_clock::now();
        norm_frame = normalizeMatrix(gray_frame);
        auto t_norm_end = std::chrono::high_resolution_clock::now();
        acc_normalize_s += std::chrono::duration_cast<std::chrono::duration<double>>(t_norm_end - t_norm_start).count();
        
        auto t_pyr_start = std::chrono::high_resolution_clock::now();
        SteerablePyramidFreq pyramid(norm_frame, nscale, norientation - 1);
        auto pyramid_coeffs = pyramid.getPyrCoeffs();
        auto t_pyr_end = std::chrono::high_resolution_clock::now();
        acc_pyramid_s += std::chrono::duration_cast<std::chrono::duration<double>>(t_pyr_end - t_pyr_start).count();
        
        auto t_bandproc_start = std::chrono::high_resolution_clock::now();
        for (const auto& band_pair : pyramid_coeffs) {
            BandKey band = band_pair.first;
            Matrix2D<Complex> coeffs = band_pair.second;
            Matrix2D<Complex> first_coeffs = first_pyramid_coeffs[band];
            
            Matrix2D<double> amp = magnitude(coeffs);
            Matrix2D<double> angle_curr = phase(coeffs);
            Matrix2D<double> angle_first = phase(first_coeffs);
            
            Matrix2D<double> dphase(angle_curr.rows, angle_curr.cols);
            for (int i = 0; i < angle_curr.rows; ++i) {
                for (int j = 0; j < angle_curr.cols; ++j) {
                    double diff = angle_curr.at(i, j) - angle_first.at(i, j);
                    dphase.at(i, j) = moduloPi(diff);
                }
            }
            
            Matrix2D<double> amp_squared = elementwiseMultiply(amp, amp);
            Matrix2D<double> sms = elementwiseMultiply(dphase, amp_squared);
            
            double total_amp_squared = matrixSum(amp_squared);
            double sum_sms = matrixSum(sms);
            
            if (total_amp_squared > 1e-10) {
                signals[band].push_back(sum_sms / total_amp_squared);
            } else {
                signals[band].push_back(0.0);
            }
        }
        auto t_bandproc_end = std::chrono::high_resolution_clock::now();
        acc_bandproc_s += std::chrono::duration_cast<std::chrono::duration<double>>(t_bandproc_end - t_bandproc_start).count();
        
        if (frame_count == 100) {
            auto end_time = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
            double total_time = duration.count() * (nframes / 100 - 1) / 60;
            std::cout << "\nCost time: " << duration.count() * 10 << " seconds \nRemaining time: " 
                     << total_time << " minutes = " << total_time * 60 << " seconds" << std::endl;
        }
        frame_count++;
        
        if (frame_count % 100 == 0) {
            std::cout << "\rProcessed " << frame_count << "/" << nframes << " frames" << std::flush;
        }
    }
    
    std::cout << "\nTotal frames processed: " << frame_count << std::endl;
        
    auto align_start = std::chrono::high_resolution_clock::now();
    std::vector<double> sound(frame_count, 0.0);
    
    BandKey reference_band(0, 0);
    if (signals.find(reference_band) == signals.end()) {
        reference_band = BandKey(-1, 0);
    }
    
    std::vector<double> reference_signal = signals[reference_band];
    
    for (auto& sig_pair : signals) {
        std::vector<double> sig = sig_pair.second;
        
            std::vector<double> sig_aligned = alignVectors(sig, reference_signal);
        
        for (size_t i = 0; i < sound.size() && i < sig_aligned.size(); ++i) {
            sound[i] += sig_aligned[i];
        }
    }
    auto align_end = std::chrono::high_resolution_clock::now();
    auto align_time = std::chrono::duration_cast<std::chrono::duration<double>>(align_end - align_start);
    
    auto filter_start = std::chrono::high_resolution_clock::now();
    auto sos = ButterworthFilter::butter(3, 0.02, "highpass");
    std::vector<double> filtered_sound = ButterworthFilter::sosfilt(sos, sound);
    auto filter_end = std::chrono::high_resolution_clock::now();
    auto filter_time = std::chrono::duration_cast<std::chrono::duration<double>>(filter_end - filter_start);
    
    auto scale_start = std::chrono::high_resolution_clock::now();
    filtered_sound = scaleSound(filtered_sound);
    auto scale_end = std::chrono::high_resolution_clock::now();
    auto scale_time = std::chrono::duration_cast<std::chrono::duration<double>>(scale_end - scale_start);

    std::cout << "\n\n=== Timing Report (soundFromVideoStreaming) ===" << std::endl;
    std::cout << "Init (first frame + pyramid): " << init_time.count() << " s" << std::endl;
    if (frame_count > 0) {
        std::cout << "Per-frame average (over " << frame_count << ")" << std::endl;
        std::cout << "  Load:       " << (acc_load_s / frame_count) << " s" << std::endl;
        std::cout << "  Resize:     " << (acc_resize_s / frame_count) << " s" << std::endl;
        std::cout << "  Normalize:  " << (acc_normalize_s / frame_count) << " s" << std::endl;
        std::cout << "  Pyramid:    " << (acc_pyramid_s / frame_count) << " s" << std::endl;
        std::cout << "  Band proc:  " << (acc_bandproc_s / frame_count) << " s" << std::endl;
    }
    std::cout << "Align + sum:  " << align_time.count() << " s" << std::endl;
    std::cout << "Filter:       " << filter_time.count() << " s" << std::endl;
    std::cout << "Scale:        " << scale_time.count() << " s" << std::endl;
    
    return filtered_sound;
}

} // namespace visualmic

