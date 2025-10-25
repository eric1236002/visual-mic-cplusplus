#include "utils.h"
#include "math_ops.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <signal_processing.h>


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
    int signal_len = signal.size();
    int kernel_len = kernel.size();
    int result_len = signal_len + kernel_len - 1;
    
    std::vector<double> result(result_len, 0.0);
    
    for (int i = 0; i < signal_len; ++i) {
        for (int j = 0; j < kernel_len; ++j) {
            result[i + j] += signal[i] * kernel[j];
        }
    }
    
    return result;
}

std::vector<double> rollVector(const std::vector<double>& vec, int shift) {
    int n = vec.size();
    if (n == 0) return vec;
    
    shift = shift % n;
    if (shift < 0) shift += n;
    
    std::vector<double> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = vec[(i - shift + n) % n];
    }
    
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
