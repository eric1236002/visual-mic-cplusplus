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
#include <omp.h>

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
    
    /*------------------------------------
    Here we initialize the signals map.
    The structure of the map is:
    signals[band][frame_idx] = signal for the band in each frame
    We also create a vector of pairs of band key and a pair of complex matrices.
    The first matrix is the coefficients of the current frame, the second matrix is the coefficients of the first frame.
    This is the preprocessing step for the parallelization of the band processing.
    ------------------------------------*/
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
        signals[pair.first] = std::vector<double>(nframes, 0.0);
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
    
    /*------------------------------------
    Here we process the frames in parallel.
    First we load the frame, then we resize it if needed, 
    then we normalize it, then we create the steerable pyramid, 
    then we process the bands in parallel. 
    Finally we align the signals and sum them.
    ------------------------------------*/
    #pragma omp parallel for schedule(dynamic) reduction(+:acc_load_s,acc_resize_s,acc_normalize_s,acc_pyramid_s,acc_bandproc_s)
    for (int frame_idx = 0; frame_idx < nframes; ++frame_idx) {
        const std::string& frame_file = frame_files[frame_idx];

        auto t_load_start = std::chrono::high_resolution_clock::now();
        Matrix2D<double> gray_frame_local = loadPGMFrame(frame_file);
        auto t_load_end = std::chrono::high_resolution_clock::now();
        acc_load_s += std::chrono::duration_cast<std::chrono::duration<double>>(t_load_end - t_load_start).count();
        
        if (downsample_factor < 1.0) {
            auto t_resize_start = std::chrono::high_resolution_clock::now();
            gray_frame_local = resizeImage(gray_frame_local, downsample_factor);
            auto t_resize_end = std::chrono::high_resolution_clock::now();
            acc_resize_s += std::chrono::duration_cast<std::chrono::duration<double>>(t_resize_end - t_resize_start).count();
        }
        
        auto t_norm_start = std::chrono::high_resolution_clock::now();
        Matrix2D<double> norm_frame_local = normalizeMatrix(gray_frame_local);
        auto t_norm_end = std::chrono::high_resolution_clock::now();
        acc_normalize_s += std::chrono::duration_cast<std::chrono::duration<double>>(t_norm_end - t_norm_start).count();
        
        auto t_pyr_start = std::chrono::high_resolution_clock::now();
        SteerablePyramidFreq pyramid(norm_frame_local, nscale, norientation - 1);
        auto pyramid_coeffs = pyramid.getPyrCoeffs();
        auto t_pyr_end = std::chrono::high_resolution_clock::now();
        acc_pyramid_s += std::chrono::duration_cast<std::chrono::duration<double>>(t_pyr_end - t_pyr_start).count();
        
        auto t_bandproc_start = std::chrono::high_resolution_clock::now();


        /*------------------------------------
        Here we create a vector of pairs of band key and a pair of complex matrices.
        The first matrix is the coefficients of the current frame, the second matrix is the coefficients of the first frame.
        This is the preprocessing step for the parallelization of the band processing.
        ------------------------------------*/
        std::vector<std::pair<BandKey, std::pair<Matrix2D<Complex>, Matrix2D<Complex>>>> band_tasks;
        band_tasks.reserve(pyramid_coeffs.size());
        for (const auto& band_pair : pyramid_coeffs) {
            BandKey band = band_pair.first;
            Matrix2D<Complex> coeffs = band_pair.second;
            Matrix2D<Complex> first_coeffs = first_pyramid_coeffs.at(band);
            band_tasks.emplace_back(band, std::make_pair(coeffs, first_coeffs));
        }

        #pragma omp parallel for schedule(dynamic) if(!omp_in_parallel())
        for (size_t task_idx = 0; task_idx < band_tasks.size(); ++task_idx) {
            BandKey band = band_tasks[task_idx].first;
            Matrix2D<Complex> coeffs = band_tasks[task_idx].second.first;
            Matrix2D<Complex> first_coeffs = band_tasks[task_idx].second.second;

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

            double signal_value = (total_amp_squared > 1e-10) ? (sum_sms / total_amp_squared) : 0.0;

            signals[band][frame_idx] = signal_value;
        }

        auto t_bandproc_end = std::chrono::high_resolution_clock::now();
        acc_bandproc_s += std::chrono::duration_cast<std::chrono::duration<double>>(t_bandproc_end - t_bandproc_start).count();
    }

    frame_count = nframes;
    
    std::cout << "\nTotal frames processed: " << frame_count << std::endl;
        
    auto align_start = std::chrono::high_resolution_clock::now();
    std::vector<double> sound(frame_count, 0.0);
    
    BandKey reference_band(0, 0);
    if (signals.find(reference_band) == signals.end()) {
        reference_band = BandKey(-1, 0);
    }
    
    std::vector<double> reference_signal = signals[reference_band];
    
    /*
        for (auto& sig_pair : signals) {
        std::vector<double> sig = sig_pair.second;
        
            std::vector<double> sig_aligned = alignVectors(sig, reference_signal);
        
        for (size_t i = 0; i < sound.size() && i < sig_aligned.size(); ++i) {
            sound[i] += sig_aligned[i];
        }
    }*/
    
    std::vector<std::pair<BandKey, std::vector<double>>> signal_pairs(signals.begin(), signals.end());

    for (size_t idx = 0; idx < signal_pairs.size(); ++idx) {
        signal_pairs[idx].second = alignVectors(signal_pairs[idx].second, reference_signal);
    }
    
    for (const auto& sig_pair : signal_pairs) {
        const std::vector<double>& sig_aligned = sig_pair.second;
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
        std::cout << "  Load:       " << acc_load_s << " s" << std::endl;
        std::cout << "  Resize:     " << acc_resize_s << " s" << std::endl;
        std::cout << "  Normalize:  " << acc_normalize_s << " s" << std::endl;
        std::cout << "  Pyramid:    " << acc_pyramid_s << " s" << std::endl;
        std::cout << "  Band proc:  " << acc_bandproc_s << " s" << std::endl;
    }
    std::cout << "Align + sum:  " << align_time.count() << " s" << std::endl;
    std::cout << "Filter:       " << filter_time.count() << " s" << std::endl;
    std::cout << "Scale:        " << scale_time.count() << " s" << std::endl;
    
    return filtered_sound;
}

} // namespace visualmic

