#include "hw1.h"

namespace algebra {
    Matrix zeros(size_t n, size_t m) {
        Matrix matrix(n, std::vector<double>(m, 0.0));
        return std::move(matrix);
    }

    Matrix ones(size_t n, size_t m) {
        Matrix matrix(n, std::vector<double>(m, 1.0));
        return std::move(matrix);
    }

    Matrix random(size_t n, size_t m, double min, double max) {
        if (max - min < 0) {
            throw std::domain_error("min number should less than max number!");
        }
        const double interval = max - min;
        std::default_random_engine e;
        Matrix matrix{};
        for (size_t i = 0; i < n; i++) {
            std::vector<double> row{};
            for (size_t j = 0; j < m; j++) {
                double random_number = static_cast<double>(e() % 100) / 100.0;
                row.emplace_back(random_number * interval + min);
            }
            matrix.emplace_back(row);
        }
        return std::move(matrix);
    }

    void show(const Matrix& matrix) {
        for (const auto& row : matrix) {
            for (const double number : row) {
                std::cout << number << std::setprecision(3) << " ";
            }
            std::cout << std::endl;
        }
    }

    Matrix sum(const Matrix& matrix, double c) {
        if (matrix.empty()) {
            return Matrix{};
        }

        Matrix new_mat = zeros(matrix.size(), matrix[0].size());
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[i].size(); j++) {
                new_mat[i][j] = matrix[i][j] + c;
            }
        }
        return new_mat;
    }

    Matrix sum(const Matrix& matrix1, const Matrix& matrix2) {
        if (matrix1.size() != matrix2.size()) {
            throw std::length_error("The length of row are different! Please check them!");
        }
        
        if (matrix1.empty() and matrix2.empty()) {
            return Matrix{};
        } else if (!matrix1.empty() and matrix2.empty()) {
            return Matrix(matrix1);
        } else if (matrix1.empty() and !matrix2.empty()) {
            return Matrix(matrix2);
        }

        Matrix new_mat = zeros(matrix1.size(), matrix1[0].size());
        // Matrix new_mat = zeros(matrix.size(), matrix[0].size());
        for (size_t i = 0; i < matrix1.size(); i++) {
            if (matrix1[i].size() != matrix2[i].size()) {
                throw std::length_error("The length of column are different! Please check them!");
            }
            for (size_t j = 0; j < matrix2.size(); j++) {
                new_mat[i][j] = matrix1[i][j] + matrix2[i][j];
            }
        }
        return new_mat;
    }

    Matrix multiply(const Matrix& matrix, double c) {
        if (matrix.empty()) {
            return Matrix{};
        }

        Matrix new_mat = zeros(matrix.size(), matrix[0].size());
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix[0].size(); j++) {
                new_mat[i][j] = matrix[i][j] * c;
            }
        }
        return new_mat;
    }

    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2) {
        if (matrix1.empty() or matrix2.empty()) {
            return Matrix{};
        }

        if (matrix1[0].size() != matrix2.size()) {
            throw std::length_error("The length of matrix1 column is not equal the length of matrix2 row");
        }
        

        Matrix new_mat = zeros(matrix1.size(), matrix2[0].size());
        for (size_t i = 0; i < matrix1.size(); i++) {
            for (size_t j = 0; j < matrix2[0].size(); j++) {
                for (size_t k = 0; k < matrix2.size(); k++) {
                    new_mat[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
        return new_mat;
    }

    Matrix transpose(const Matrix& matrix) {
        if (matrix.empty()) {
            return Matrix{};
        }

        Matrix new_mat = zeros(matrix[0].size(), matrix.size());
        for (size_t i{}; i < matrix.size(); i++) {
            for (size_t j{}; j < matrix[i].size(); j++) {
                new_mat[j][i] = matrix[i][j];
            }
        }
        return new_mat;
    }

    Matrix minor(const Matrix& matrix, size_t n, size_t m) {
        if (matrix.empty()) {
            return Matrix{};
        }

        if (n >= matrix.size() or m >= matrix[0].size()) {
            throw std::domain_error("parameter error! check n or m");
        }

        Matrix ans = zeros(matrix.size() - 1, matrix[0].size() - 1);
        size_t row_delta = 0;
        for (size_t i{}; i < matrix.size(); i++) {
            if (i == n) {
                row_delta = 1;
                continue;
            }
            size_t col_delta = 0;
            for (size_t j{}; j < matrix[i].size(); j++) {
                if (j == m) {
                    col_delta = 1;
                    continue;
                }
                ans[i - row_delta][j - col_delta] = matrix[i][j];
            }
        }
        return ans;
    }

    double determinant(const Matrix& matrix) {
        if (matrix.size() == 1) {
            return matrix[0][0];
        }

        if (matrix.size() == 0) {
            return 1;
        }
        else if (matrix.size() != matrix[0].size()) {
            throw std::domain_error("matrix error! check col and row");
        }

        double determinant_sum = 0.0;
        for (size_t i{}; i < matrix.size(); i++) {
            double sign = i % 2 == 0 ? 1.0 : -1.0;
            determinant_sum += matrix[i][0] * determinant(minor(matrix, i, 0)) * sign;
        }
        return determinant_sum;
    }

    Matrix inverse(const Matrix& matrix) {
        if (matrix.empty()) {
            return Matrix{};
        } else if (matrix.size() != matrix[0].size()) {
            throw std::domain_error("matrix size error! check col and row");
        }

        if (fabs(determinant(matrix)) <= 1e-6) {
            throw std::domain_error("matrix error! check col and row");
        }

        // generate the adjoint matrix
        Matrix adjoint_matrix = zeros(matrix.size(), matrix[0].size());
        for (size_t i{}; i < matrix.size(); i++) {
            for (size_t j{}; j < matrix[i].size(); j++) {
                adjoint_matrix[i][j] = determinant(minor(matrix, i, j));
            }
        }
        return multiply(transpose(adjoint_matrix), 1.0 / determinant(matrix));
    }

    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis=0) {
        if (matrix1.empty()) {
            return Matrix(matrix2);
        }
        if (matrix2.empty()) {
            return Matrix(matrix1);
        }
        
        if (axis == 0) {
            if (matrix1[0].size() != matrix2[0].size()) {
                throw std::length_error("the matrix should have same col length.");
            }
            Matrix new_mat = zeros(matrix1.size() + matrix2.size(), matrix1[0].size());
            for (size_t i{}; i < matrix1[0].size(); i++) {
                for (size_t j{}; j < matrix1.size(); j++) {
                    new_mat[j][i] = matrix1[j][i];
                }
                for (size_t j{}; j < matrix2.size(); j++) {
                    new_mat[matrix1.size() + j][i] = matrix2[j][i];
                }
            }
            return new_mat;
        } else if (axis == 1) {
            if (matrix1.size() != matrix2.size()) {
                throw std::length_error("the matrix should have same rpw length.");
            }
            Matrix new_mat = zeros(matrix1.size(), matrix1[0].size() + matrix2[0].size());
            for (size_t i{}; i < matrix1.size(); i++) {
                for (size_t j{}; j < matrix1[0].size(); j++) {
                    new_mat[i][j] = matrix1[i][j];
                }
                for (size_t j{}; j < matrix2[0].size(); j++) {
                    new_mat[i][j + matrix1[0].size()] = matrix2[i][j];
                }
            }
            return new_mat;
        } else {
            throw std::domain_error("axis should is either 0 or 1.");
        }
    }

    Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2) {
        if (matrix.empty()) {
            return Matrix{};
        }

        if (r1 >= matrix.size() or r2 >= matrix.size()) {
            throw std::length_error("r1 or r2 length are great than matrix row size!");
        }

        Matrix new_mat(matrix);
        std::vector<double> tmp(matrix[r1]);
        new_mat[r1] = new_mat[r2];
        new_mat[r2] = tmp;
        return new_mat;
    }

    Matrix ero_multiply(const Matrix& matrix, size_t r, double c) {
        if (matrix.empty()) {
            return Matrix{};
        }

        if (r >= matrix.size()) {
            throw std::length_error("r1 or r2 length are great than matrix row size!");
        }

        Matrix new_mat(matrix);
        for (size_t i = 0; i < new_mat[r].size(); i++) {
            new_mat[r][i] *= c;
        }
        return new_mat;
    }

    Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2) {
        if (matrix.empty()) {
            return Matrix{};
        }

        if (r1 >= matrix.size() or r2 >= matrix.size()) {
            throw std::length_error("r1 or r2 length are great than matrix row size!");
        }

        Matrix new_mat(matrix);
        for (size_t i = 0; i < new_mat[r1].size(); i++) {
            new_mat[r2][i] += new_mat[r1][i] * c;
        }
        return new_mat;
    }

    Matrix upper_triangular(const Matrix& matrix) {
        if (matrix.empty()) {
            return Matrix{};
        }

        if (matrix.size() != matrix[0].size()) {
            throw std::length_error("matrix should have same row and col length!");
        }

        Matrix new_mat(matrix);
        for (size_t i{}; i < matrix.size(); i++) {
            if (fabs(new_mat[i][i]) <= 1e-6) {
                // find the first non-zero row and swpa it.
                size_t pos = i + 1;
                while (pos < matrix.size()) {
                    if (fabs(new_mat[pos][pos]) > 1e-6) {
                        break;
                    }
                }
                // If all number is zero, we should jump this col.
                if (pos == matrix.size()) {
                    continue;
                } else {
                    new_mat = ero_swap(new_mat, i, pos);
                }
            }
            for (size_t j = i + 1; j < matrix.size(); j++) {
                new_mat = ero_sum(new_mat, i, -new_mat[j][i] / new_mat[i][i], j);
            }
        }

        return new_mat;
    }
}
