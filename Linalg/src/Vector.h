#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <time.h>

namespace Linalg
{
    template <typename T>
    class Vector
    {
    private:
        std::vector<T> m_vector;
        int m_length = 0;

    public:
        /// @brief operator of sum of two vectors
        /// @param vector1 vector
        /// @param vector2 vector
        /// @return sum of two vectors
        friend Vector<T> operator+(Vector<T> vector1, Vector<T> vector2)
        {
            return vector1.add(vector2);
        }

        /// @brief operator of defference of two vectors
        /// @param vector1 vector
        /// @param vector2 vector
        /// @return defference of two vectors
        friend Vector<T> operator-(Vector<T> vector1, Vector<T> vector2)
        {
            return vector1.subtract(vector2);
        }

        /// @brief operator that returns scalar multiplied vector
        /// @param scalar number
        /// @param vector vector
        /// @return scalar multiplied vector
        friend Vector<T> operator*(const T scalar, Vector<T> vector)
        {
            return vector.scalarMultiply(scalar);
        }

        /// @brief operator that returns scalar multiplied vector
        /// @param vector vector
        /// @param scalar number
        /// @return scalar multiplied vector
        friend Vector<T> operator*(Vector<T> vector, const T scalar)
        {
            return vector.scalarMultiply(scalar);
        }

        /// @brief operator that computes inner product of two vectors
        /// @param vector1 vector
        /// @param vector2 vector
        /// @return inner product
        friend T operator*(Vector<T> vector1, Vector<T> vector2)
        {
            return vector1.innerProduct(vector2);
        }

        /// @brief operator that divides vector by number
        /// @param vector vector
        /// @param num number
        /// @return vector devided by number
        friend Vector<T> operator/(Vector<T> vector, const T num)
        {
            if (num == 0)
            {
                throw std::runtime_error("zero-division is not defined");
            }
            return vector.scalarMultiply(1.0 / num);
        }

        /// @brief operator that computes element-wise mod
        /// @param vector vector
        /// @param num number
        /// @return vector mod num
        friend Vector<T> operator%(Vector<T> vector, const int num)
        {
            if (num == 0)
            {
                throw std::runtime_error("zero-division is not defined");
            }
            Vector<T> mod_vector;
            mod_vector.setLength(vector.length());
            for (int i = 0; i < vector.length(); ++i)
            {
                mod_vector.substitute(i, vector.at(i) % num);
            }
            return mod_vector;
        }

        friend Vector<T> operator+=(Vector<T> &vector1, Vector<T> vector2)
        {
            vector1 = vector1 + vector2;
            return vector1;
        }

        friend Vector<T> operator-=(Vector<T> &vector1, Vector<T> vector2)
        {
            vector1 = vector1 - vector2;
            return vector1;
        }

        friend Vector<T> operator*=(Vector<T> &vector, const T scalar)
        {
            vector = vector * scalar;
            return vector;
        }

        friend Vector<T> operator/=(Vector<T> &vector, const T num)
        {
            vector = vector / num;
            return vector;
        }

        /// @brief function that prints a vector
        void print()
        {
            for (int i = 0; i < m_length; ++i)
            {
                std::cout << m_vector.at(i) << " ";
            }
            std::cout << std::endl;
        }

        /// @brief function that sets length of a vector
        /// @param length length of vector
        void setLength(const int length)
        {
            m_vector.resize(length);
            m_length = length;
        }

        /// @brief function that returns length of a vector
        /// @return length of vector
        int length()
        {
            return m_length;
        }

        /// @brief function that casts a vector
        /// @tparam U type to cast
        /// @return a vector casted
        template <typename U>
        Vector<U> cast()
        {
            Vector<U> vector;
            vector.setLength(m_length);
            for (int i = 0; i < m_length; ++i)
            {
                vector.substitute(i, static_cast<U>(m_vector.at(i)));
            }
            return vector;
        }

        /// @brief substitute number to vector coordinate
        /// @param index index
        /// @param number_to_sbstitute number to substitute to vector
        void substitute(const int index, const T number_to_sbstitute)
        {
            m_vector.at(index) = number_to_sbstitute;
        }

        /// @brief function that substitute the vector to a vector
        /// @param vector a vector to substitute
        void substitute(std::vector<T> vector)
        {
            m_length = vector.size();
            m_vector.resize(m_length);
            for (int i = 0; i < m_length; ++i)
            {
                m_vector.at(i) = vector.at(i);
            }
        }

        /// @brief function that returns index-th component
        /// @param index index
        /// @return index-th component
        T at(const int index)
        {
            if (index >= m_length)
            {
                throw std::runtime_error("index " + std::to_string(index) + " is out of range");
            }
            return m_vector.at(index);
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
            for (int i = 0; i < m_length; ++i)
            {
                m_vector.at(i) = rand();
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

            for (int i = 0; i < m_length; ++i)
            {
                m_vector.at(i) = lower_bound + (upper_bound - lower_bound) * static_cast<double>(rand()) / RAND_MAX;
            }
        }

        /// @brief function that tests a vector is zero with error ``error``
        /// @param error error
        /// @return returns true if vector is zero with error ``error`` else false
        bool isZero(const double error = 0)
        {
            for (int i = 0; i < m_length; ++i)
            {
                if (std::abs(m_vector.at(i)) > error)
                {
                    return false;
                }
            }
            return true;
        }

        /// @brief function that adds two vectors
        /// @param vector vector to add
        /// @return sum of two vectors
        Vector<T> add(Vector<T> vector)
        {
            if (m_length != vector.length())
            {
                std::string error_massage = "sum of a vector whose length is " + std::to_string(m_length) + " and a vector whose length is " + std::to_string(vector.length()) + " is not defined.";
                throw std::runtime_error(error_massage);
            }

            Vector<T> added_vector;
            added_vector.setLength(m_length);
            for (int i = 0; i < m_length; ++i)
            {
                added_vector.substitute(i, m_vector.at(i) + vector.at(i));
            }
            return added_vector;
        }

        Vector<T> subtract(Vector<T> vector)
        {
            if (m_length != vector.length())
            {
                std::string error_massage = "difference of a vector whose length is " + std::to_string(m_length) + " and a vector whose length is " + std::to_string(vector.length()) + " is not defined.";
                throw std::runtime_error(error_massage);
            }
            Vector<T> subtracted_vector;
            subtracted_vector.setLength(m_length);
            for (int i = 0; i < m_length; ++i)
            {
                subtracted_vector.substitute(i, m_vector.at(i) - vector.at(i));
            }
            return subtracted_vector;
        }

        /// @brief function that multiplies scalar to a vector
        /// @param scalar number
        /// @return a vector multiplied by scalar
        Vector<T> scalarMultiply(const T scalar)
        {
            Vector<T> scalar_multiplied_vector;
            scalar_multiplied_vector.setLength(m_length);
            for (int i = 0; i < m_length; ++i)
            {
                scalar_multiplied_vector.substitute(i, scalar * m_vector.at(i));
            }
            return scalar_multiplied_vector;
        }

        /// @brief function that computes inner product of two vectors
        /// @param vector vector
        /// @return innner product
        T innerProduct(Vector<T> vector)
        {
            if (m_length != vector.length())
            {
                std::string error_massage = "inner product of a vector whose length is " + std::to_string(m_length) + " and a vector whose length is " + std::to_string(vector.length()) + " is not defined.";
                throw std::runtime_error(error_massage);
            }
            T S = 0;
            for (int i = 0; i < m_length; ++i)
            {
                S += m_vector.at(i) * vector.at(i);
            }
            return S;
        }

        /// @brief function that computes lp-norm of a vector
        /// @param p number
        /// @return lp-norm
        double norm(const unsigned int p = 2)
        {
            double S = 0.0;
            for (int i = 0; i < m_length; ++i)
            {
                S += std::pow(m_vector.at(i), p);
            }
            return std::pow(S, 1.0 / p);
        }
    };
}

#endif // !VECTOR_H
