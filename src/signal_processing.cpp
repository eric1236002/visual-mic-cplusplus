#include "signal_processing.h"
#include "math_ops.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace visualmic {

// Butterworth Filter Implementation
ButterworthFilter::SOSCoefficients ButterworthFilter::butter(
    int order, double cutoff, const std::string& btype) {
    
    SOSCoefficients sos;
    
    //directly transfer the coefficients from the python implementation
    if (btype == "highpass" && order == 3) {
        sos.sos = {
            {0.85449723, -0.85449723,  0,  1, -0.85408069,0},
            {1.0, -2.0, 1.0, 1.0, -1.9777, 0.9780}
        };
    } else {
        throw std::runtime_error("Only 3rd order highpass filter at 0.05 is implemented");
    }
    
    return sos;
}

std::vector<double> ButterworthFilter::sosfilt(
    const SOSCoefficients& sos, 
    const std::vector<double>& signal) {
    
    auto output = signal;
    
    for (const auto& section : sos.sos) {
        std::vector<double> temp = output;
        
        double b0 = section[0], b1 = section[1], b2 = section[2];
        double a0 = section[3], a1 = section[4], a2 = section[5];
        
        double w1 = 0.0, w2 = 0.0;
        
        for (size_t i = 0; i < temp.size(); ++i) {
            double w0 = temp[i] - a1 * w1 - a2 * w2;
            output[i] = b0 * w0 + b1 * w1 + b2 * w2;
            w2 = w1;
            w1 = w0;
        }
    }
    
    return output;
}


std::vector<double> getSoundSpectralSubtraction(const std::vector<double>& signal) {
    if (signal.empty()) return signal;
    
    int N = signal.size();
    const double floorGain = 0.1;  
    
    int fftSize = 1;
    while (fftSize < N) fftSize <<= 1;
    
    std::vector<Complex> X(fftSize, Complex(0, 0));
    for (int i = 0; i < N; ++i) {
        X[i] = Complex(signal[i], 0.0);
    }
    fft1D(X);
    
    int noiseLen = fftSize / 10;
    std::vector<double> noiseMag(fftSize, 0.0);
    for (int k = 0; k < fftSize; ++k) {
        noiseMag[k] = (k < noiseLen) ? std::abs(X[k]) : 0.0;
    }
    
    for (int k = 0; k < fftSize; ++k) {
        double mag = std::abs(X[k]);
        double phase = std::arg(X[k]);
        double noise = (k < noiseLen) ? noiseMag[k] : noiseMag[noiseLen-1];
        double newMag = std::max(mag - noise, mag * floorGain);
        X[k] = Complex(newMag * std::cos(phase), newMag * std::sin(phase));
    }
    
    ifft1D(X);
    
    std::vector<double> output(N);
    for (int i = 0; i < N; ++i) {
        output[i] = X[i].real();
    }
    
    return output;
}

} // namespace visualmic

