#include "utils.h"
#include "math_ops.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <signal_processing.h>
#include <mpi.h>


namespace visualmic {

std::vector<double> scaleSound(const std::vector<double>& sound) {
    if (sound.empty()) return sound;
    
    auto result = sound;
    double maxv = *std::max_element(result.begin(), result.end());
    double minv = *std::min_element(result.begin(), result.end());
    
    if (maxv != 1.0 || minv != -1.0) {
        double rangev = maxv - minv;
        if (rangev > 1e-10) {
            for (auto& val : result) {
                val = 2.0 * val / rangev;
            }
            maxv = *std::max_element(result.begin(), result.end());
            double offset = maxv - 1.0;
            for (auto& val : result) {
                val -= offset;
            }
        }
    }
    
    return result;
}

std::vector<double> flipVector(const std::vector<double>& vec) {
    std::vector<double> flipped(vec.rbegin(), vec.rend());
    return flipped;
}

std::vector<double> convolve(const std::vector<double>& signal, 
                              const std::vector<double>& kernel) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int signal_len = signal.size();
    int kernel_len = kernel.size();
    int result_len = signal_len + kernel_len - 1;
    // Distribute work across processes
    int chunk_size = (signal_len + size - 1) / size;
    int start_idx = rank * chunk_size;
    int end_idx = std::min(start_idx + chunk_size, signal_len);
    int local_len = std::max(0, end_idx - start_idx);
    
    // Local convolution computation
    std::vector<double> local_result(result_len, 0.0);
    
    for (int i = start_idx; i < end_idx; ++i) {
        for (int j = 0; j < kernel_len; ++j) {
            local_result[i + j] += signal[i] * kernel[j];
        }
    }
    
    // Reduce all local results to get final result
    std::vector<double> result(result_len, 0.0);
    MPI_Allreduce(local_result.data(), result.data(), result_len, 
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return result;
}

// MPI-parallel rollVector
std::vector<double> rollVector(const std::vector<double>& vec, int shift) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int n = vec.size();
    if (n == 0) return vec;
    
    shift = shift % n;
    if (shift < 0) shift += n;
    
    // Distribute indices across processes
    int chunk_size = (n + size - 1) / size;
    int start_idx = rank * chunk_size;
    int end_idx = std::min(start_idx + chunk_size, n);
    
    // Each process computes its portion of the result
    std::vector<double> local_result(n, 0.0);
    
    for (int i = start_idx; i < end_idx; ++i) {
        local_result[i] = vec[(i - shift + n) % n];
    }
    
    // Gather results from all processes
    // Note: We can't use Allreduce with SUM here because each index is written by only one process
    // Instead, we use Allreduce with MAX (since non-owning processes have 0.0)
    // Better approach: use a custom operation or recognize that only one process writes each index
    std::vector<double> result(n, 0.0);
    MPI_Allreduce(local_result.data(), result.data(), n, 
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return result;
}

std::vector<double> alignVectors(const std::vector<double>& v1, 
                                  const std::vector<double>& v2) {
    auto v2_flipped = flipVector(v2);
    auto acorb = convolve(v1, v2_flipped);
    
    auto max_it = std::max_element(acorb.begin(), acorb.end());
    int maxind = std::distance(acorb.begin(), max_it);
    
    int shift = v2.size() - 1 - maxind;
    
    return rollVector(v1, shift);
}

double moduloPi(double angle) {
    double result = std::fmod(angle + M_PI, 2.0 * M_PI);
    if (result < 0) result += 2.0 * M_PI;
    return result - M_PI;
}

Matrix2D<double> absMatrix(const Matrix2D<Complex>& mat) {
    return magnitude(mat);
}

Matrix2D<double> angleMatrix(const Matrix2D<Complex>& mat) {
    return phase(mat);
}

} // namespace visualmic
