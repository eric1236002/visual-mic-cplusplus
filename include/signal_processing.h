#ifndef SIGNAL_PROCESSING_H
#define SIGNAL_PROCESSING_H

#include <vector>
#include <complex>

namespace visualmic {

class ButterworthFilter {
public:
    struct SOSCoefficients {
        std::vector<std::vector<double>> sos; 
    };
    
    static SOSCoefficients butter(int order, double cutoff, 
                                  const std::string& btype = "highpass");
    
    static std::vector<double> sosfilt(const SOSCoefficients& sos, 
                                       const std::vector<double>& signal);
};

std::vector<double> getSoundSpectralSubtraction(const std::vector<double>& signal);

} // namespace visualmic

#endif // SIGNAL_PROCESSING_H

