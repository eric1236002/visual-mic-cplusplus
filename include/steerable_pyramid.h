#ifndef STEERABLE_PYRAMID_H
#define STEERABLE_PYRAMID_H

#include <vector>
#include <complex>
#include <map>
#include <tuple>
#include <pthread.h>
#include "math_ops.h"

namespace visualmic {

using BandKey = std::tuple<int, int>;

class SteerablePyramidFreq {
public:
    SteerablePyramidFreq(const Matrix2D<double>& image, int nscale, int norient);
    ~SteerablePyramidFreq();
    
    const std::map<BandKey, Matrix2D<Complex>>& getPyrCoeffs() const { return pyr_coeffs; }
    
private:
    std::map<BandKey, Matrix2D<Complex>> pyr_coeffs;
    int num_scales;
    int num_orientations;
    pthread_mutex_t pyr_mutex;

    struct ThreadData {
        SteerablePyramidFreq* pyramid;
        const Matrix2D<Complex>* fft_image;
        Matrix2D<double>* highpass_mask;
        const Matrix2D<double>* radial_mask;
        int scale;
        int orient;
    };
    
    Matrix2D<double> buildSteerableFilters(int rows, int cols, int orientation, int norient);
    Matrix2D<double> buildRadialMask(int rows, int cols, int level, int nscales);
    Matrix2D<double> buildAngularMask(int rows, int cols, int orientation, int norient);
    
    void createFrequencyGrid(int rows, int cols, Matrix2D<double>& fx, Matrix2D<double>& fy);
    
    Matrix2D<Complex> applyFrequencyFilter(const Matrix2D<Complex>& fft_image, 
                                            const Matrix2D<double>& filter);

    static void* process_orientation_thread(void* arg);
};

} // namespace visualmic

#endif // STEERABLE_PYRAMID_H
