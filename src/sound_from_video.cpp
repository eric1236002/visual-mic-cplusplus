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
#include <unistd.h>
#include <vector>

namespace visualmic {

Matrix2D<double> resizeImage(const Matrix2D<double>& img, double scale_factor) {
    int new_rows = static_cast<int>(img.rows * scale_factor);
    int new_cols = static_cast<int>(img.cols * scale_factor);
    return resizeMatrix(img, new_rows, new_cols);
}

void* process_frames_thread(void* thread_arg) {
    ThreadData* data = static_cast<ThreadData*>(thread_arg);
    auto first_pyramid_coeffs = data->first_pyramid->getPyrCoeffs();
    
    // Initialize timing accumulators
    data->total_load_time = 0.0;
    data->total_resize_time = 0.0;
    data->total_normalize_time = 0.0;
    data->total_pyramid_time = 0.0;
    data->total_bandproc_time = 0.0;
    data->frames_processed = 0;

    for (int i = data->start_frame; i < data->end_frame; ++i) {
        const std::string& frame_file = (*data->frame_files)[i];
        
        auto load_start = std::chrono::high_resolution_clock::now();
        Matrix2D<double> gray_frame = loadPGMFrame(frame_file);
        auto load_end = std::chrono::high_resolution_clock::now();
        data->total_load_time += std::chrono::duration_cast<std::chrono::duration<double>>(load_end - load_start).count();

        auto resize_start = std::chrono::high_resolution_clock::now();
        if (data->downsample_factor < 1.0) {
            gray_frame = resizeImage(gray_frame, data->downsample_factor);
        }
        auto resize_end = std::chrono::high_resolution_clock::now();
        data->total_resize_time += std::chrono::duration_cast<std::chrono::duration<double>>(resize_end - resize_start).count();

        auto normalize_start = std::chrono::high_resolution_clock::now();
        Matrix2D<double> norm_frame = normalizeMatrix(gray_frame);
        auto normalize_end = std::chrono::high_resolution_clock::now();
        data->total_normalize_time += std::chrono::duration_cast<std::chrono::duration<double>>(normalize_end - normalize_start).count();
        
        auto pyramid_start = std::chrono::high_resolution_clock::now();
        SteerablePyramidFreq pyramid(norm_frame, data->nscale, data->norientation - 1);
        auto pyramid_coeffs = pyramid.getPyrCoeffs();
        auto pyramid_end = std::chrono::high_resolution_clock::now();
        data->total_pyramid_time += std::chrono::duration_cast<std::chrono::duration<double>>(pyramid_end - pyramid_start).count();

        auto bandproc_start = std::chrono::high_resolution_clock::now();
        for (const auto& band_pair : pyramid_coeffs) {
            BandKey band = band_pair.first;
            Matrix2D<Complex> coeffs = band_pair.second;
            Matrix2D<Complex> first_coeffs = first_pyramid_coeffs.at(band);

            Matrix2D<double> amp = magnitude(coeffs);
            Matrix2D<double> angle_curr = phase(coeffs);
            Matrix2D<double> angle_first = phase(first_coeffs);

            Matrix2D<double> dphase(angle_curr.rows, angle_curr.cols);
            for (int r = 0; r < angle_curr.rows; ++r) {
                for (int c = 0; c < angle_curr.cols; ++c) {
                    double diff = angle_curr.at(r, c) - angle_first.at(r, c);
                    dphase.at(r, c) = moduloPi(diff);
                }
            }

            Matrix2D<double> amp_squared = elementwiseMultiply(amp, amp);
            Matrix2D<double> sms = elementwiseMultiply(dphase, amp_squared);

            double total_amp_squared = matrixSum(amp_squared);
            double sum_sms = matrixSum(sms);

            if (total_amp_squared > 1e-10) {
                data->signals[band].push_back(sum_sms / total_amp_squared);
            } else {
                data->signals[band].push_back(0.0);
            }
        }
        auto bandproc_end = std::chrono::high_resolution_clock::now();
        data->total_bandproc_time += std::chrono::duration_cast<std::chrono::duration<double>>(bandproc_end - bandproc_start).count();
        
        data->frames_processed++;
    }
    return nullptr;
}


std::vector<double> soundFromVideoStreaming(const std::string& frames_dir,
                                            int nscale, 
                                            int norientation, 
                                            double downsample_factor,
                                            int num_threads) {
    
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
    auto init_end = std::chrono::high_resolution_clock::now();
    auto init_time = std::chrono::duration_cast<std::chrono::duration<double>>(init_end - init_start);
    
    std::cout << "Processing frames (streaming mode - low memory usage)..." << std::endl;

    // Auto-detect CPU cores if num_threads is 0 or negative
    if (num_threads <= 0) {
        num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    }
    
    std::cout << "Using " << num_threads << " threads for processing" << std::endl;
    std::vector<ThreadData> thread_data(num_threads);
    int frames_per_thread = nframes / num_threads;

    for (int i = 0; i < num_threads; ++i) {
        thread_data[i].frame_files = &frame_files;
        thread_data[i].start_frame = i * frames_per_thread;
        thread_data[i].end_frame = (i == num_threads - 1) ? nframes : (i + 1) * frames_per_thread;
        thread_data[i].nscale = nscale;
        thread_data[i].norientation = norientation;
        thread_data[i].downsample_factor = downsample_factor;
        thread_data[i].first_pyramid = &first_pyramid;
        
        for (const auto& pair : first_pyramid.getPyrCoeffs()) {
            thread_data[i].signals[pair.first] = std::vector<double>();
        }

        int rc = pthread_create(&thread_data[i].thread_id, nullptr, process_frames_thread, &thread_data[i]);
        if (rc) {
            std::cerr << "Error:unable to create thread," << rc << std::endl;
            exit(-1);
        }
    }

    std::map<BandKey, std::vector<double>> signals;
    double acc_load_s = 0.0;
    double acc_resize_s = 0.0;
    double acc_normalize_s = 0.0;
    double acc_pyramid_s = 0.0;
    double acc_bandproc_s = 0.0;
    int frame_count = 0;
    
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(thread_data[i].thread_id, nullptr);
        
        // Aggregate timing data
        acc_load_s += thread_data[i].total_load_time;
        acc_resize_s += thread_data[i].total_resize_time;
        acc_normalize_s += thread_data[i].total_normalize_time;
        acc_pyramid_s += thread_data[i].total_pyramid_time;
        acc_bandproc_s += thread_data[i].total_bandproc_time;
        frame_count += thread_data[i].frames_processed;
        
        // Aggregate signal data
        for (const auto& pair : thread_data[i].signals) {
            signals[pair.first].insert(signals[pair.first].end(), pair.second.begin(), pair.second.end());
        }
    }
    

    
    std::cout << "\nTotal frames processed: " << nframes << std::endl;
        
    auto align_start = std::chrono::high_resolution_clock::now();
    std::vector<double> sound(nframes, 0.0);
    
    BandKey reference_band(0, 0);
    if (signals.find(reference_band) == signals.end()) {
        reference_band = BandKey(-1, 0);
    }
    
    std::vector<double> reference_signal = signals[reference_band];
    
    for (auto& sig_pair : signals) {
        std::vector<double> sig = sig_pair.second;
        
        std::vector<double> sig_aligned = alignVectors(sig, reference_signal, num_threads);
        
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
        std::cout << "Total processing time (over " << frame_count << " frames)" << std::endl;
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

