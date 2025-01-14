#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string.h>

#include "Vector.h"

#define EPSILON 1e-6

namespace Linalg
{
    template <typename T>
    class Matrix
    {
    private:
        std::vector<std::vector<T>> m_matrix;
        int m_num_rows = 0;
        int m_num_columns = 0;

    public:
        friend Matrix<T> operator+(Matrix<T> matrix1, Matrix<T> matrix2)
        {
            if ((matrix1.rows() != matrix2.rows()) || (matrix1.columns() != matrix2.columns()))
            {
                throw std::runtime_error("sum of matrix whose size is (" + std::to_string(matrix1.rows()) + ", " + std::to_string(matrix1.columns()) + ") and matrix whose size is (" + std::to_string(matrix1.rows()) + ", " + std::to_string(matrix1.columns()) + ") is not defined.");
            }

            Matrix<T> matrix;
            matrix.setDimension(matrix1.rows(), matrix2.columns());

            for (int i = 0, j; i < matrix1.rows(); ++i)
            {
                for (j = 0; j < matrix2.columns(); ++j)
                {
                    matrix.substitute(i, j, matrix1.at(i, j) + matrix2.at(i, j));
                }
            }

            return matrix;
        }

        friend Matrix<T> operator-(Matrix<T> matrix1, Matrix<T> matrix2)
        {
            if ((matrix1.rows() != matrix2.rows()) || (matrix1.columns() != matrix2.columns()))
            {
                throw std::runtime_error("difference of matrix whose size is (" + std::to_string(matrix1.rows()) + ", " + std::to_string(matrix1.columns()) + ") and matrix whose size is (" + std::to_string(matrix1.rows()) + ", " + std::to_string(matrix1.columns()) + ") is not defined.");
            }

            Matrix<T> matrix;
            matrix.setDimension(matrix1.rows(), matrix2.columns());

            for (int i = 0, j; i < matrix1.row(); ++i)
            {
                for (j = 0; j < matrix2.columns(); ++j)
                {
                    matrix.substitute(i, j, matrix1.at(i, j) - matrix2.at(i, j));
                }
            }

            return matrix;
        }

        friend Matrix<T> operator*(Matrix<T> matrix1, Matrix<T> matrix2)
        {
            if (matrix1.columns() != matrix2.rows())
            {
                throw std::runtime_error("product of matrix whose size is (" + std::to_string(matrix1.rows()) + ", " + std::to_string(matrix1.columns()) + ") and matrix whose size is (" + std::to_string(matrix1.rows()) + ", " + std::to_string(matrix1.columns()) + ") is not defined.");
            }

            Matrix<T> matrix;
            T S = 0;
            matrix.setDimension(matrix1.rows(), matrix2.columns());

            for (int i = 0, j, k; i < matrix1.rows(); ++i)
            {
                for (j = 0; j < matrix2.columns(); ++j)
                {
                    S = 0;
                    for (k = 0; k < matrix1.columns(); ++k)
                    {
                        S += matrix1.at(i, k) * matrix2.at(k, j);
                    }
                    matrix.substitute(i, j, S);
                }
            }

            return matrix;
        }

        friend Matrix<T> operator*(const T scalar, Matrix<T> matrix)
        {
            Matrix<T> temp_matrix;
            temp_matrix.setDimension(matrix.rows(), matrix.columns());

            for (int i = 0, j; i < matrix.rows(); ++i)
            {
                for (j = 0; j < matrix.columns(); ++j)
                {
                    temp_matrix.substitute(i, j, scalar * matrix.at(i, j));
                }
            }

            return temp_matrix;
        }

        friend Matrix<T> operator*(Matrix<T> matrix, const T scalar)
        {
            Matrix<T> temp_matrix;
            temp_matrix.setDimension(matrix.rows(), matrix.columns());

            for (int i = 0, j; i < matrix.rows(); ++i)
            {
                for (j = 0; j < matrix.columns(); ++j)
                {
                    temp_matrix.substitute(i, j, scalar * matrix.at(i, j));
                }
            }

            return temp_matrix;
        }

        friend Matrix<T> operator/(Matrix<T> matrix, const T num)
        {
            if (num == 0)
            {
                throw std::runtime_error("zero-division is not defined");
            }

            Matrix<T> temp_matrix;
            temp_matrix.setDimension(matrix.rows(), matrix.columns());
            for (int i = 0, j; i < matrix.rows(); ++i)
            {
                for (j = 0; j < matrix.columns(); ++j)
                {
                    temp_matrix.substitute(i, j, matrix.at(i, j) / num);
                }
            }
            return temp_matrix;
        }

        friend Matrix<T> operator+=(Matrix<T> &matrix1, Matrix<T> matrix2)
        {
            matrix1 = matrix1 + matrix2;
            return matrix1;
        }

        friend Matrix<T> operator*=(Matrix<T> &matrix1, Matrix<T> matrix2)
        {
            matrix1 = matrix1 * matrix2;
            return matrix1;
        }

        friend Matrix<T> operator*=(Matrix<T> &matrix1, const T scalar)
        {
            matrix1 = matrix1 * scalar;
            return matrix1;
        }

        friend Matrix<T> operator/=(Matrix<T> &matrix, const T num)
        {
            matrix = matrix / num;
            return matrix;
        }

        /// @brief function that computes n-th power of a matrix
        /// @param n number
        /// @return n-th power of a matrix
        Matrix<T> power(const int n)
        {
            if (!isSquare())
            {
                throw std::runtime_error("power of non-square matrix is not defined");
            }

            Matrix<T> matrix, temp_matrix;
            matrix.setIdentity(m_num_rows);
            temp_matrix.substitute(m_matrix);

            if (n > 0)
            {
                for (int i = 0; i < n; ++i)
                {
                    matrix *= temp_matrix;
                }
            }

            return matrix;
        }

        /// @brief function that prints a matrix
        void print()
        {
            std::cout << "[";
            for (int i = 0, j; i < m_num_rows; ++i)
            {
                if (i == 0)
                {
                    std::cout << "[";
                }
                else
                {
                    std::cout << "\n[";
                }

                for (j = 0; j < m_num_columns; ++j)
                {
                    std::cout << m_matrix.at(i).at(j) << ", ";
                }
                std::cout << "\b\b]" << ", ";
            }
            std::cout << "\b\b]" << std::endl;
        }

        /// @brief function that prints a matrix by assigned format
        /// @param print_mode format to print
        void print(const char *print_mode)
        {
            if (strcmp(print_mode, "python") == 0 || strcmp(print_mode, "py") == 0)
            {
                std::cout << "[";
                for (int i = 0, j; i < m_num_rows; ++i)
                {
                    std::cout << "[";

                    for (j = 0; j < m_num_columns; ++j)
                    {
                        std::cout << m_matrix.at(i).at(j) << ", ";
                    }
                    std::cout << "\b\b]" << ", ";
                }
                std::cout << "\b\b]" << std::endl;
            }
            else if (strcmp(print_mode, "c") == 0 || strcmp(print_mode, "C") == 0 || strcmp(print_mode, "cpp") == 0)
            {
                std::cout << "{";
                for (int i = 0, j; i < m_num_rows; ++i)
                {
                    std::cout << "{";

                    for (j = 0; j < m_num_columns; ++j)
                    {
                        std::cout << m_matrix.at(i).at(j) << ", ";
                    }
                    std::cout << "\b\b}" << ", ";
                }
                std::cout << "\b\b}" << std::endl;
            }
            else if (strcmp(print_mode, "TeX") == 0)
            {
                std::cout << "\\begin{pmatrix}\n";
                for (int i = 0, j; i < m_num_rows; ++i)
                {
                    for (j = 0; j < m_num_columns; ++j)
                    {
                        std::cout << m_matrix.at(i).at(j) << " & ";
                    }

                    if (i == m_num_rows - 1)
                    {
                        std::cout << "\b\b" << "   \n";
                    }
                    else
                    {
                        std::cout << "\b\b\b" << "\\\\\n";
                    }
                }
                std::cout << "\b\b\\end{pmatrix}" << std::endl;
            }
        }

        /// @brief function that returns (``row_index``, ``column_index``)-element of a matrix
        /// @param row_index row index
        /// @param column_index column index
        /// @return (``row_index``, ``column_index``)-element of a matrix
        T at(const int row_index, const int column_index)
        {
            if (row_index >= m_num_rows)
            {
                throw std::runtime_error("index " + std::to_string(row_index) + " is out of range");
            }
            else if (column_index >= m_num_columns)
            {
                throw std::runtime_error("index " + std::to_string(column_index) + " is out of range");
            }

            return m_matrix.at(row_index).at(column_index);
        }

        /// @brief function that sets seed of random by now time
        void setSeed()
        {
            srand(static_cast<unsigned int>(time(NULL)));
        }

        /// @brief function that sets seed of random
        /// @param seed seed
        void setSeed(const unsigned int seed)
        {
            srand(seed);
        }

        /// @brief function that generates random vector
        void random()
        {
            for (int i = 0, j; i < m_num_rows; ++i)
            {
                for (j = 0; j < m_num_columns; ++j)
                {
                    substitute(i, j, rand());
                }
            }
        }

        /// @brief function that generates random vector
        /// @param lower_bound lower bound of elements
        /// @param upper_bound upper bound of elements
        void random(const double lower_bound, const double upper_bound)
        {
            if (lower_bound > upper_bound)
            {
                throw std::runtime_error("lower bound is greater than upper bound");
            }

            for (int i = 0, j; i < m_num_rows; ++i)
            {
                for (j = 0; j < m_num_columns; ++j)
                {
                    substitute(i, j, lower_bound + (upper_bound - lower_bound) * static_cast<double>(rand()) / RAND_MAX);
                }
            }
        }

        /// @brief function that sets row and column
        /// @param num_rows number of rows
        /// @param num_columns number of columns
        void setDimension(const int num_rows, const int num_columns)
        {
            m_num_rows = num_rows;
            m_num_columns = num_columns;
            m_matrix.resize(num_rows);
            for (int i = 0; i < num_rows; ++i)
            {
                m_matrix.at(i).resize(num_columns);
            }
        }

        /// @brief function that sets identity matrix
        /// @param size
        void setIdentity(const int size)
        {
            m_num_rows = size;
            m_num_columns = size;
            m_matrix.resize(size);
            for (int i = 0; i < size; ++i)
            {
                m_matrix.at(i).resize(size);
                m_matrix.at(i).at(i) = 1;
            }
        }

        /// @brief function that sets a diagonal elements of a matrix by a vector
        /// @param diagonal_vector diagonal elements
        void setDiagonal(Vector<T> diagonal_vector)
        {
            setDimension(diagonal_vector.length(), diagonal_vector.length());
            for (int i = 0, j; i < rows(); ++i)
            {
                for (j = 0; j < columns(); ++j)
                {
                    if (i == j)
                    {
                        substitute(i, j, diagonal_vector.at(i));
                    }
                    else
                    {
                        substitute(i, j, 0);
                    }
                }
            }
        }

        /// @brief function that tests that a matrix is regular or not
        /// @return returns true if a matrix is regular, false if is not regular
        bool isRegular()
        {
            if (m_num_rows != m_num_columns)
            {
                return false;
            }
            else if (std::abs(det()) < EPSILON)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        /// @brief function that tests a matrix is square or not
        /// @return true if a matrix is square, false if is not square
        bool isSquare()
        {
            if (m_num_columns == m_num_rows)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        /// @brief function that tests that a matrix is dianonal or not
        /// @return true if a matrix is diagonal, false if is not diagonal
        bool isDiagonal()
        {
            if (!isSquare())
            {
                return false;
            }
            else
            {
                for (int i = 0, j; i < m_num_rows; ++i)
                {
                    for (j = 0; j < m_num_columns; ++j)
                    {
                        if (i != j)
                        {
                            if (std::abs(at(i, j)) >= EPSILON)
                            {
                                return false;
                            }
                        }
                    }
                }
            }
            return true;
        }
        
        /// @brief function that tests that a matrix is uni-modular or not
        /// @return true if a matrix is uni-modular, false if is not uni-modular
        bool isUniModular()
        {
            if (!isRegular)
            {
                return false;
            }

            for (int i = 0, j; i < m_num_rows; ++i)
            {
                for (j = 0; j < m_num_columns; ++j)
                {
                    if (std::abs(round(at(i, j)) - at(i, j)) >= EPSILON)
                    {
                        return false;
                    }
                }
            }

            T temp = determinant();

            if (std::abs(temp - 1) < EPSILON || std::abs(temp + 1) < EPSILON)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        /// @brief function that substitute std::vector<std::vector<T>> to Matirx
        /// @param matrix vector of vector
        void substitute(const std::vector<std::vector<T>> matrix)
        {
            m_num_rows = matrix.size();
            m_num_columns = matrix.at(0).size();
            m_matrix.resize(m_num_rows);
            for (int i = 0, j; i < m_num_rows; ++i)
            {
                m_matrix.at(i).resize(m_num_columns);
                for (j = 0; j < m_num_columns; ++j)
                {
                    m_matrix.at(i).at(j) = matrix.at(i).at(j);
                }
            }
        }

        /// @brief function that substitutes a vector to ``index``-th row of matrix
        /// @param index index number
        /// @param vector a vector to substitute
        void substituteRow(const int index, std::vector<T> vector)
        {
            if (index >= m_num_rows)
            {
                throw std::runtime_error("index " + std::to_string(index) + " is out of range");
            }

            for (int i = 0; i < m_num_columns; ++i)
            {
                m_matrix.at(index).at(i) = vector.at(i);
            }
        }

        /// @brief function that substitutes a vector to ``index``-th row of matrix
        /// @param index index number
        /// @param vector a vector to substitute
        void substituteRow(const int index, Vector<T> vector)
        {
            if (index >= m_num_rows)
            {
                throw std::runtime_error("index " + std::to_string(index) + " is out of range");
            }

            for (int i = 0; i < m_num_columns; ++i)
            {
                m_matrix.at(index).at(i) = vector.at(i);
            }
        }

        /// @brief function that substitutes a vector to ``index``-th column of matrix
        /// @param index index number
        /// @param vector a vector to substitute
        void substituteColumn(const int index, std::vector<T> vector)
        {
            if (index >= m_num_columns)
            {
                throw std::runtime_error("index " + std::to_string(index) + " is out of range");
            }

            for (int i = 0; i < m_num_rows; ++i)
            {
                m_matrix.at(i).at(index) = vector.at(i);
            }
        }

        /// @brief function that substitutes a vector to ``index``-th column of matrix
        /// @param index index number
        /// @param vector a vector to substitute
        void substituteColumn(const int index, Vector<T> vector)
        {
            if (index >= m_num_columns)
            {
                throw std::runtime_error("index " + std::to_string(index) + " is out of range");
            }

            for (int i = 0; i < m_num_rows; ++i)
            {
                m_matrix.at(i).at(index) = vector.at(i);
            }
        }

        /// @brief function that substitutes number to (``row_index``, ``column_index``)-component of matrix
        /// @param row_index index number of row
        /// @param column_index index number of column
        /// @param number number to substitute
        void substitute(const int row_index, const int column_index, const T number)
        {
            if (row_index >= m_num_rows)
            {
                throw std::runtime_error("index " + std::to_string(row_index) + " is out of range");
            }
            else if (column_index >= m_num_columns)
            {
                throw std::runtime_error("index " + std::to_string(column_index) + " is out of range");
            }

            m_matrix.at(row_index).at(column_index) = number;
        }

        /// @brief function that casts matrix
        /// @tparam U type to cast
        /// @return casted matrix
        template <typename U>
        Matrix<U> cast()
        {
            Matrix<U> matrix;
            matrix.setDimension(m_num_rows, m_num_columns);
            for (int i = 0, j; i < m_num_rows; ++i)
            {
                for (j = 0; j < m_num_columns; ++j)
                {
                    matrix.substitute(i, j, static_cast<U>(m_matrix.at(i).at(j)));
                }
            }
            return matrix;
        }

        /// @brief function that returns ``index``-th row of matrix
        /// @param index index
        /// @return ``index``-th row of matrix
        Vector<T> row(const int index)
        {
            if (index >= m_num_rows)
            {
                throw std::runtime_error("index " + std::to_string(index) + " is out of range");
            }

            Vector<T> row_vector;
            row_vector.setLength(m_num_columns);
            for (int i = 0; i < m_num_columns; ++i)
            {
                row_vector.substitute(i, m_matrix.at(index).at(i));
            }
            return row_vector;
        }

        /// @brief function that returns number of rows
        int rows()
        {
            return m_num_rows;
        }

        /// @brief function that returns number of columns
        int columns()
        {
            return m_num_columns;
        }

        /// @brief computes determinant of a matrix
        /// @return determinant of a matrix
        T determinant()
        {
            if (m_num_rows != m_num_columns)
            {
                throw std::runtime_error("determinant of non-square matrix is not defined.");
            }

            T S = 0;
            Matrix<T> temp_matrix;
            temp_matrix.setDimension(m_num_rows - 1, m_num_columns - 1);

            if (m_num_rows == 1)
            {
                return at(0, 0);
            }
            else
            {
                for (int i = 0, j, k; i < m_num_columns; ++i)
                {
                    for (j = 1; j < m_num_rows; ++j)
                    {
                        for (k = 0; k < m_num_columns; ++k)
                        {
                            if (k < i)
                            {
                                temp_matrix.substitute(j - 1, k, at(j, k));
                            }
                            else if (k > i)
                            {
                                temp_matrix.substitute(j - 1, k - 1, at(j, k));
                            }
                        }
                    }

                    if (i % 2 == 1)
                    {
                        S += at(0, i) * temp_matrix.det();
                    }
                    else
                    {
                        S -= at(0, i) * temp_matrix.det();
                    }
                }
            }

            return S;
        }

        /// @brief function that computes determinant of a matirx
        /// @return determinant of a matirx
        T det()
        {
            return determinant();
        }

        /// @brief computes transpose of a matrix
        /// @return transpose of a matrix
        Matrix<T> transpose()
        {
            Matrix<T> matrix;
            matrix.setDimension(m_num_columns, m_num_rows);
            for (int i = 0, j; i < m_num_columns; ++i)
            {
                for (j = 0; j < m_num_rows; ++j)
                {
                    matrix.substitute(i, j, m_matrix.at(j).at(i));
                }
            }
            return matrix;
        }

        /// @brief computes transpose of a matrix
        /// @return transpose of a matrix
        Matrix<T> trans()
        {
            return transpose();
        }

        /// @brief computes trace of a matrix
        /// @return trace of a matrix
        T trace()
        {
            T S = 0;
            for (int i = 0; i < std::min(m_num_rows, m_num_columns); ++i)
            {
                S += at(i, i);
            }
            return S;
        }

        /// @brief function that computes (i, j)-co-factor
        /// @param row_index row index
        /// @param column_index column index
        /// @return (i, j)-co-factor
        T coFactor(const int row_index, const int column_index)
        {
            T temp_determinant;
            Matrix<T> matrix;
            matrix.setDimension(m_num_rows - 1, m_num_columns - 1);

            for (int i = 0, j; i < m_num_rows; ++i)
            {
                for (j = 0; j < m_num_columns; ++j)
                {
                    if (i < row_index)
                    {
                        if (j < column_index)
                        {
                            matrix.substitute(i, j, m_matrix.at(i).at(j));
                        }
                        else if (j > column_index)
                        {
                            matrix.substitute(i, j - 1, m_matrix.at(i).at(j));
                        }
                    }
                    else if (i > row_index)
                    {
                        if (j < column_index)
                        {
                            matrix.substitute(i - 1, j, m_matrix.at(i).at(j));
                        }
                        else if (j > column_index)
                        {
                            matrix.substitute(i - 1, j - 1, m_matrix.at(i).at(j));
                        }
                    }
                }
            }

            temp_determinant = matrix.det();

            if ((row_index + column_index + 1) % 2)
            {
                temp_determinant = -temp_determinant;
            }

            return temp_determinant;
        }

        /// @brief function that computes adjugate of a matrix
        /// @return adjugate matrix
        Matrix<T> adjugate()
        {
            Matrix<T> matrix;
            matrix.setDimension(m_num_rows, m_num_columns);

            for (int i = 0, j; i < m_num_rows; ++i)
            {
                for (j = 0; j < m_num_columns; ++j)
                {
                    matrix.substitute(j, i, coFactor(i, j));
                }
            }
            return matrix;
        }

        /// @brief function that computes inverse of a matrix
        /// @return inverse matrix
        Matrix<double> inverse()
        {
            if (isRegular())
            {
                throw std::runtime_error("inverse matrix of non-regular matrix is not defined");
            }

            Matrix<T> matrix;
            matrix.setDimension(m_num_rows, m_num_columns);
            matrix = adjugate();

            return matrix.cast<double>() * (1.0 / det());
        }

        /// @brief function that computes Gram-Schmidt orthogonalized vector matrix
        /// @return Gram-Schmidt orthogonalized vector matrix
        Matrix<double> GramSchmidt()
        {
            Matrix<T> basis;
            Matrix<double> gso_vector_matrix, gso_coeff_matrix;
            basis.setDimension(m_num_rows, m_num_columns);
            basis.substitute(m_matrix);
            gso_vector_matrix.setDimension(m_num_rows, m_num_columns);
            gso_coeff_matrix.setIdentity(m_num_rows);

            for (int i = 0, j, k; i < m_num_rows; ++i)
            {
                for (j = 0; j < m_num_rows; ++j)
                {
                    gso_vector_matrix.substitute(i, j, basis.cast<double>().at(i, j));
                }

                for (j = 0; j < i; ++j)
                {
                    gso_coeff_matrix.substitute(i, j, basis.cast<double>().row(i).innerProduct(gso_vector_matrix.row(j)) / gso_vector_matrix.row(j).innerProduct(gso_vector_matrix.row(j)));
                    for (k = 0; k < m_num_columns; ++k)
                    {
                        gso_vector_matrix.substitute(i, k, gso_vector_matrix.at(i, k) - gso_coeff_matrix.at(i, j) * gso_vector_matrix.at(j, k));
                    }
                }
            }
            return gso_vector_matrix;
        }
    };
}

#endif // !MATRIX_H
