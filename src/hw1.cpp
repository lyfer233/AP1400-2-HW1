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

    Matrix sum(Matrix& matrix, double c) {
        if (matrix.empty()) {
            return Matrix{};
        }

        Matrix new_mat = zeros(matrix.size(), matrix[0].size());
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix.size(); j++) {
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
}
