#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <complex>
#include "math_ops.h"

namespace visualmic {

using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;
using ComplexMatrix = std::vector<std::vector<Complex>>;

std::vector<double> scaleSound(const std::vector<double>& sound);

std::vector<double> alignVectors(const std::vector<double>& v1, 
                                  const std::vector<double>& v2);

std::vector<double> convolve(const std::vector<double>& signal, 
                              const std::vector<double>& kernel);

std::vector<double> flipVector(const std::vector<double>& vec);

std::vector<double> rollVector(const std::vector<double>& vec, int shift);

double moduloPi(double angle);

Matrix2D<double> absMatrix(const Matrix2D<Complex>& mat);
Matrix2D<double> angleMatrix(const Matrix2D<Complex>& mat);

} // namespace visualmic

#endif // UTILS_H
