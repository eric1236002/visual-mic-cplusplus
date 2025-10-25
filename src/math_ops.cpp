#include "math_ops.h"
#include <algorithm>
#include <stdexcept>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace visualmic {

void fft1D(std::vector<Complex>& x) {
    int N = x.size();
    if (N <= 1) return;
    
    if ((N & (N - 1)) != 0) {
        int N_padded = 1;
        while (N_padded < N) N_padded <<= 1;
        x.resize(N_padded, Complex(0, 0));
        N = N_padded;
    }
    
    for (int i = 0, j = 0; i < N; ++i) {
        if (j > i) std::swap(x[i], x[j]);
        int m = N >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    
    for (int s = 2; s <= N; s *= 2) {
        double angle = -2.0 * M_PI / s;
        Complex w(std::cos(angle), std::sin(angle));
        
        for (int k = 0; k < N; k += s) {
            Complex wn(1.0, 0.0);
            for (int j = 0; j < s / 2; ++j) {
                Complex t = wn * x[k + j + s / 2];
                Complex u = x[k + j];
                x[k + j] = u + t;
                x[k + j + s / 2] = u - t;
                wn *= w;
            }
        }
    }
}

void ifft1D(std::vector<Complex>& x) {
    int N = x.size();
    
    for (auto& val : x) {
        val = std::conj(val);
    }
    
    fft1D(x);
    
    for (auto& val : x) {
        val = std::conj(val) / static_cast<double>(N);
    }
}

void fft2D(Matrix2D<Complex>& mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    
    for (int i = 0; i < rows; ++i) {
        std::vector<Complex> row(cols);
        for (int j = 0; j < cols; ++j) {
            row[j] = mat.at(i, j);
        }
        fft1D(row);
        for (int j = 0; j < cols; ++j) {
            mat.at(i, j) = row[j];
        }
    }
    
    for (int j = 0; j < cols; ++j) {
        std::vector<Complex> col(rows);
        for (int i = 0; i < rows; ++i) {
            col[i] = mat.at(i, j);
        }
        fft1D(col);
        for (int i = 0; i < rows; ++i) {
            mat.at(i, j) = col[i];
        }
    }
}

void ifft2D(Matrix2D<Complex>& mat) {
    int rows = mat.rows;
    int cols = mat.cols;
    
    for (int i = 0; i < rows; ++i) {
        std::vector<Complex> row(cols);
        for (int j = 0; j < cols; ++j) {
            row[j] = mat.at(i, j);
        }
        ifft1D(row);
        for (int j = 0; j < cols; ++j) {
            mat.at(i, j) = row[j];
        }
    }
    
    for (int j = 0; j < cols; ++j) {
        std::vector<Complex> col(rows);
        for (int i = 0; i < rows; ++i) {
            col[i] = mat.at(i, j);
        }
        ifft1D(col);
        for (int i = 0; i < rows; ++i) {
            mat.at(i, j) = col[i];
        }
    }
}

Matrix2D<double> elementwiseMultiply(const Matrix2D<double>& a, const Matrix2D<double>& b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        throw std::runtime_error("Matrix dimensions must match for elementwise multiplication");
    }
    
    Matrix2D<double> result(a.rows, a.cols);
    for (int i = 0; i < a.rows; ++i) {
        for (int j = 0; j < a.cols; ++j) {
            result.at(i, j) = a.at(i, j) * b.at(i, j);
        }
    }
    return result;
}

Matrix2D<Complex> elementwiseMultiply(const Matrix2D<Complex>& a, const Matrix2D<double>& b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        throw std::runtime_error("Matrix dimensions must match for elementwise multiplication");
    }
    
    Matrix2D<Complex> result(a.rows, a.cols);
    for (int i = 0; i < a.rows; ++i) {
        for (int j = 0; j < a.cols; ++j) {
            result.at(i, j) = a.at(i, j) * b.at(i, j);
        }
    }
    return result;
}

Matrix2D<double> matrixAdd(const Matrix2D<double>& a, const Matrix2D<double>& b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        throw std::runtime_error("Matrix dimensions must match for addition");
    }
    
    Matrix2D<double> result(a.rows, a.cols);
    for (int i = 0; i < a.rows; ++i) {
        for (int j = 0; j < a.cols; ++j) {
            result.at(i, j) = a.at(i, j) + b.at(i, j);
        }
    }
    return result;
}

Matrix2D<double> matrixSubtract(const Matrix2D<double>& a, const Matrix2D<double>& b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        throw std::runtime_error("Matrix dimensions must match for subtraction");
    }
    
    Matrix2D<double> result(a.rows, a.cols);
    for (int i = 0; i < a.rows; ++i) {
        for (int j = 0; j < a.cols; ++j) {
            result.at(i, j) = a.at(i, j) - b.at(i, j);
        }
    }
    return result;
}

Matrix2D<double> magnitude(const Matrix2D<Complex>& mat) {
    Matrix2D<double> result(mat.rows, mat.cols);
    for (int i = 0; i < mat.rows; ++i) {
        for (int j = 0; j < mat.cols; ++j) {
            result.at(i, j) = std::abs(mat.at(i, j));
        }
    }
    return result;
}

Matrix2D<double> phase(const Matrix2D<Complex>& mat) {
    Matrix2D<double> result(mat.rows, mat.cols);
    for (int i = 0; i < mat.rows; ++i) {
        for (int j = 0; j < mat.cols; ++j) {
            result.at(i, j) = std::arg(mat.at(i, j));
        }
    }
    return result;
}

Matrix2D<double> normalizeMatrix(const Matrix2D<double>& mat) {
    if (mat.rows == 0 || mat.cols == 0) {
        return mat;
    }
    
    double min_val = mat.data[0][0];
    double max_val = mat.data[0][0];
    
    for (int i = 0; i < mat.rows; ++i) {
        for (int j = 0; j < mat.cols; ++j) {
            double val = mat.at(i, j);
            if (val < min_val) min_val = val;
            if (val > max_val) max_val = val;
        }
    }
    
    Matrix2D<double> result(mat.rows, mat.cols);
    double range = max_val - min_val;
    
    if (range < 1e-10) {
        return result;
    }
    
    for (int i = 0; i < mat.rows; ++i) {
        for (int j = 0; j < mat.cols; ++j) {
            result.at(i, j) = (mat.at(i, j) - min_val) / range;
        }
    }
    
    return result;
}

Matrix2D<double> resizeMatrix(const Matrix2D<double>& mat, int new_rows, int new_cols) {
    Matrix2D<double> result(new_rows, new_cols);
    
    double row_ratio = static_cast<double>(mat.rows) / new_rows;
    double col_ratio = static_cast<double>(mat.cols) / new_cols;
    
    for (int i = 0; i < new_rows; ++i) {
        for (int j = 0; j < new_cols; ++j) {
            double src_row = i * row_ratio;
            double src_col = j * col_ratio;
            
            int r0 = static_cast<int>(src_row);
            int c0 = static_cast<int>(src_col);
            int r1 = std::min(r0 + 1, mat.rows - 1);
            int c1 = std::min(c0 + 1, mat.cols - 1);
            
            double dr = src_row - r0;
            double dc = src_col - c0;
            
            double v00 = mat.data[r0][c0];
            double v01 = mat.data[r0][c1];
            double v10 = mat.data[r1][c0];
            double v11 = mat.data[r1][c1];
            
            double v0 = v00 * (1 - dc) + v01 * dc;
            double v1 = v10 * (1 - dc) + v11 * dc;
            
            result.at(i, j) = v0 * (1 - dr) + v1 * dr;
        }
    }
    
    return result;
}

Matrix2D<double> rgbToGray(const Matrix2D<std::vector<double>>& rgb_mat) {
    Matrix2D<double> result(rgb_mat.rows, rgb_mat.cols);
    
    for (int i = 0; i < rgb_mat.rows; ++i) {
        for (int j = 0; j < rgb_mat.cols; ++j) {
            const auto& pixel = rgb_mat.data[i][j];
            if (pixel.size() >= 3) {
                result.at(i, j) = 0.299 * pixel[2] + 0.587 * pixel[1] + 0.114 * pixel[0];
            }
        }
    }
    
    return result;
}

double matrixSum(const Matrix2D<double>& mat) {
    double sum = 0.0;
    for (int i = 0; i < mat.rows; ++i) {
        for (int j = 0; j < mat.cols; ++j) {
            sum += mat.at(i, j);
        }
    }
    return sum;
}

double matrixMean(const Matrix2D<double>& mat) {
    if (mat.rows == 0 || mat.cols == 0) return 0.0;
    return matrixSum(mat) / (mat.rows * mat.cols);
}

} // namespace visualmic
