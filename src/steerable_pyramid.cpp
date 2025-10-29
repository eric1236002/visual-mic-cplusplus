#include "steerable_pyramid.h"
#include "math_ops.h"
#include <cmath>
#include <iostream>
#include <signal_processing.h>
#include <omp.h>

namespace visualmic {


SteerablePyramidFreq::SteerablePyramidFreq(const Matrix2D<double>& image, int nscale, int norient)
    : num_scales(nscale), num_orientations(norient + 1) {
    
    Matrix2D<double> img_double = image;
    
    int rows = img_double.rows;
    int cols = img_double.cols;
    
    Matrix2D<Complex> img_complex(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            img_complex.at(i, j) = Complex(img_double.at(i, j), 0.0);
        }
    }
    
    Matrix2D<Complex> fft_image = img_complex;
    fft2D(fft_image);
    
    Matrix2D<double> fx, fy;
    createFrequencyGrid(rows, cols, fx, fy);
    
    Matrix2D<double> highpass_mask(rows, cols, 1.0);
    
    for (int scale = 0; scale < num_scales; ++scale) {
        Matrix2D<double> radial_mask = buildRadialMask(rows, cols, scale, num_scales);
        
        for (int orient = 0; orient < num_orientations; ++orient) {
            Matrix2D<double> angular_mask = buildAngularMask(rows, cols, orient, num_orientations);
            Matrix2D<double> combined_mask = elementwiseMultiply(radial_mask, angular_mask);
            Matrix2D<Complex> band_coeffs = applyFrequencyFilter(fft_image, combined_mask);
            
            pyr_coeffs[BandKey(scale, orient)] = band_coeffs;
            highpass_mask = matrixSubtract(highpass_mask, combined_mask);
        }
    }
    
    Matrix2D<double> lowpass_mask = buildRadialMask(rows, cols, num_scales, num_scales);
    Matrix2D<Complex> lowpass_coeffs = applyFrequencyFilter(fft_image, lowpass_mask);
    pyr_coeffs[BandKey(-1, 0)] = lowpass_coeffs;
    
    Matrix2D<Complex> highpass_coeffs = applyFrequencyFilter(fft_image, highpass_mask);
    pyr_coeffs[BandKey(0, -1)] = highpass_coeffs;
}

void SteerablePyramidFreq::createFrequencyGrid(int rows, int cols, 
                                                Matrix2D<double>& fx, 
                                                Matrix2D<double>& fy) {
    fx = Matrix2D<double>(rows, cols);
    fy = Matrix2D<double>(rows, cols);
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double u = (j < cols / 2) ? j : j - cols;
            double v = (i < rows / 2) ? i : i - rows;
            
            fx.at(i, j) = u / cols;
            fy.at(i, j) = v / rows;
        }
    }
}

Matrix2D<double> SteerablePyramidFreq::buildRadialMask(int rows, int cols, 
                                                        int level, int nscales) {
    Matrix2D<double> fx, fy;
    createFrequencyGrid(rows, cols, fx, fy);
    
    Matrix2D<double> radius(rows, cols);
    #pragma omp parallel for collapse(2)schedule(static)
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            radius.at(i, j) = std::sqrt(fx.at(i, j) * fx.at(i, j) + 
                                        fy.at(i, j) * fy.at(i, j));
        }
    }
    
    Matrix2D<double> mask(rows, cols, 0.0);
    
    double scale_factor = std::pow(2.0, level);
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double r = radius.at(i, j) * scale_factor * 2.0;
            
            if (r < 0.5) {
                mask.at(i, j) = 1.0;
            } else if (r < 1.0) {
                mask.at(i, j) = std::cos((r - 0.5) * M_PI) * 0.5 + 0.5;
            }
        }
    }
    
    return mask;
}

Matrix2D<double> SteerablePyramidFreq::buildAngularMask(int rows, int cols, 
                                                         int orientation, int norient) {
    Matrix2D<double> fx, fy;
    createFrequencyGrid(rows, cols, fx, fy);
    
    Matrix2D<double> angle(rows, cols);
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            angle.at(i, j) = std::atan2(fy.at(i, j), fx.at(i, j));
        }
    }
    
    Matrix2D<double> mask(rows, cols, 0.0);
    
    double target_angle = M_PI * orientation / norient;
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double theta = angle.at(i, j);
            double diff = std::abs(theta - target_angle);
            
            if (diff > M_PI) {
                diff = 2 * M_PI - diff;
            }
            
            double width = M_PI / norient;
            if (diff < width) {
                mask.at(i, j) = std::pow(std::cos(diff * norient / 2.0), 2.0);
            }
        }
    }
    
    return mask;
}

Matrix2D<Complex> SteerablePyramidFreq::applyFrequencyFilter(
    const Matrix2D<Complex>& fft_image, 
    const Matrix2D<double>& filter) {
    
    Matrix2D<Complex> filtered_fft = elementwiseMultiply(fft_image, filter);
    ifft2D(filtered_fft);
    
    return filtered_fft;
}

} // namespace visualmic
