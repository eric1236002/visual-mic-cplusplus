#ifndef MATH_OPS_H
#define MATH_OPS_H

#include <vector>
#include <complex>
#include <cmath>

namespace visualmic {

using Complex = std::complex<double>;

template<typename T>
class Matrix2D {
public:
    std::vector<std::vector<T>> data;
    int rows;
    int cols;
    
    Matrix2D() : rows(0), cols(0) {}
    Matrix2D(int r, int c) : rows(r), cols(c) {
        data.resize(r, std::vector<T>(c, T()));
    }
    
    Matrix2D(int r, int c, T val) : rows(r), cols(c) {
        data.resize(r, std::vector<T>(c, val));
    }
    
    T& at(int i, int j) { return data[i][j]; }
    const T& at(int i, int j) const { return data[i][j]; }
};

void fft1D(std::vector<Complex>& x);
void ifft1D(std::vector<Complex>& x);
void fft2D(Matrix2D<Complex>& mat);
void ifft2D(Matrix2D<Complex>& mat);

Matrix2D<double> matrixMultiply(const Matrix2D<double>& a, const Matrix2D<double>& b);
Matrix2D<Complex> matrixMultiply(const Matrix2D<Complex>& a, const Matrix2D<Complex>& b);
Matrix2D<double> matrixAdd(const Matrix2D<double>& a, const Matrix2D<double>& b);
Matrix2D<double> matrixSubtract(const Matrix2D<double>& a, const Matrix2D<double>& b);
Matrix2D<double> elementwiseMultiply(const Matrix2D<double>& a, const Matrix2D<double>& b);
Matrix2D<Complex> elementwiseMultiply(const Matrix2D<Complex>& a, const Matrix2D<double>& b);

Matrix2D<double> magnitude(const Matrix2D<Complex>& mat);
Matrix2D<double> phase(const Matrix2D<Complex>& mat);

Matrix2D<double> normalizeMatrix(const Matrix2D<double>& mat);
Matrix2D<double> resizeMatrix(const Matrix2D<double>& mat, int new_rows, int new_cols);
Matrix2D<double> rgbToGray(const Matrix2D<std::vector<double>>& rgb_mat);

double matrixSum(const Matrix2D<double>& mat);
double matrixMean(const Matrix2D<double>& mat);

} // namespace visualmic

#endif // MATH_OPS_H
