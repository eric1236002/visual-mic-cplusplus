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
#include <iomanip>
#include <mpi.h>

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
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    std::vector<std::string> frame_files = getFrameFilesList(frames_dir);
    
    if (frame_files.empty()) {
        throw std::runtime_error("No frames found in directory");
    }
    
    int nframes = frame_files.size();
    if (rank == 0)
        std::cout << "Found " << nframes << " frames to process (streaming mode)" << std::endl;
    
    auto init_start = std::chrono::high_resolution_clock::now();
    if (rank == 0)
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
    if (rank == 0)
        std::cout << "Processing frames (streaming mode - low memory usage)..." << std::endl;
    
    int frame_count = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = std::chrono::duration<double>::zero();
    double acc_load_s = 0.0;
    double acc_resize_s = 0.0;
    double acc_normalize_s = 0.0;
    double acc_pyramid_s = 0.0;
    double acc_bandproc_s = 0.0;
    
    // Divide frames among processes
    int frames_per_process = frame_files.size() / size;
    int remainder = frame_files.size() % size;
    int start_idx = rank * frames_per_process + std::min(rank, remainder);
    int end_idx = start_idx + frames_per_process + (rank < remainder ? 1 : 0);
    
    auto frames_processing_start = std::chrono::high_resolution_clock::now();
    for (int idx = start_idx; idx < end_idx; ++idx) {
        auto t_load_start = std::chrono::high_resolution_clock::now();
        const auto& frame_file = frame_files[idx];
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
        
        if (rank == 0 && frame_count == 100) {
            auto end_time = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
            double total_time = duration.count() * (nframes / 100 - 1) / 60;
            std::cout << "\nCost time: " << duration.count() * 10 << " seconds \nRemaining time: " 
                     << total_time << " minutes = " << total_time * 60 << " seconds" << std::endl;
        }
        frame_count++;
        
        if (rank == 0 && frame_count % 100 == 0) {
            std::cout << "\rProcessed " << frame_count << "/" << end_idx - start_idx << " frames" << std::flush;
        }
    }
    auto frames_processing_end = std::chrono::high_resolution_clock::now();
    auto frames_processing_time = std::chrono::duration_cast<std::chrono::duration<double>>(frames_processing_end - frames_processing_start);
    
    int global_frame_count = 0;
    MPI_Allreduce(&frame_count, &global_frame_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    double total_acc_load_s = 0.0;
    double total_acc_resize_s = 0.0;
    double total_acc_normalize_s = 0.0;
    double total_acc_pyramid_s = 0.0;
    double total_acc_bandproc_s = 0.0;
    double max_frames_processing_time = 0.0;
    
    double local_frames_processing_time = frames_processing_time.count();
    
    MPI_Reduce(&acc_load_s, &total_acc_load_s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&acc_resize_s, &total_acc_resize_s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&acc_normalize_s, &total_acc_normalize_s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&acc_pyramid_s, &total_acc_pyramid_s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&acc_bandproc_s, &total_acc_bandproc_s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_frames_processing_time, &max_frames_processing_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        
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
    
    /*filter the sound*/
    auto filter_start = std::chrono::high_resolution_clock::now();
    auto sos = ButterworthFilter::butter(3, 0.02, "highpass");
    std::vector<double> filtered_sound = ButterworthFilter::sosfilt(sos, sound);
    auto filter_end = std::chrono::high_resolution_clock::now();
    auto filter_time = std::chrono::duration_cast<std::chrono::duration<double>>(filter_end - filter_start);
    
    /*scale the sound*/
    auto scale_start = std::chrono::high_resolution_clock::now();
    filtered_sound = scaleSound(filtered_sound);
    auto scale_end = std::chrono::high_resolution_clock::now();
    auto scale_time = std::chrono::duration_cast<std::chrono::duration<double>>(scale_end - scale_start);
    
    double local_align_s = align_time.count();
    double local_filter_s = filter_time.count();
    double local_scale_s = scale_time.count();
    
    double total_align_s = 0.0;
    double total_filter_s = 0.0;
    double total_scale_s = 0.0;
    
    MPI_Reduce(&local_align_s, &total_align_s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_filter_s, &total_filter_s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_scale_s, &total_scale_s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "\n\n=== Timing Report (soundFromVideoStreaming) ===" << std::endl;
        std::cout << "Init (first frame + pyramid): " << std::fixed << std::setprecision(5) << init_time.count() << " s" << std::endl;
        std::cout << "Total frames processed: " << global_frame_count << " using " << size << " MPI ranks" << std::endl;
        
        if (global_frame_count > 0) {
            std::cout << "\nPre-processing (parallel region)" << std::endl;
            std::cout << "  Wall-clock (max rank): " << std::fixed << std::setprecision(5) << max_frames_processing_time << " s" << std::endl;
            std::cout << "  Avg per frame (wall-clock): " << std::fixed << std::setprecision(5) << (max_frames_processing_time / global_frame_count) << " s" << std::endl;
            
            auto printModule = [&](const std::string& name, double total_cpu_s) {
                double avg_ms = (total_cpu_s / global_frame_count) * 1000.0;
                double total_ms = total_cpu_s * 1000.0;
                std::cout << "  " << name << ": avg " << std::setw(9) << std::fixed << std::setprecision(4)
                          << avg_ms << " ms  | total " << std::setw(10) << total_ms << " ms" << std::endl;
            };
            
            printModule("Load     ", total_acc_load_s);
            if (downsample_factor < 1.0) {
                printModule("Resize   ", total_acc_resize_s);
            }
            printModule("Normalize", total_acc_normalize_s);
            printModule("Pyramid  ", total_acc_pyramid_s);
            printModule("Band proc", total_acc_bandproc_s);
        }
        
        std::cout << "\nPost-processing (aggregated CPU time)" << std::endl;
        double total_postproc_s = total_align_s + total_filter_s + total_scale_s;
        std::cout << "  Total: " << std::fixed << std::setprecision(5) << total_postproc_s << " s" << std::endl;
        if (global_frame_count > 0) {
            std::cout << "  Align + sum:  " << std::fixed << std::setprecision(5) << (total_align_s / global_frame_count) * 1000 << " ms/frame" << std::endl;
            std::cout << "  Filter:       " << std::fixed << std::setprecision(5) << (total_filter_s / global_frame_count) * 1000 << " ms/frame" << std::endl;
            std::cout << "  Scale:        " << std::fixed << std::setprecision(5) << (total_scale_s / global_frame_count) * 1000 << " ms/frame" << std::endl;
        }
    }
    
    return filtered_sound;
}

} // namespace visualmic

