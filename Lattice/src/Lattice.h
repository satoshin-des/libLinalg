#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <vector>

#include "../../Linalg/Core"

namespace LATTICE
{
	template <typename T>
	class Lattice
	{
	private:
		int m_number_of_row = 0;							// 基底行列の行数
		int m_number_of_column = 0;							// 基底行列の列数
		std::vector<std::vector<T>> m_basis;				// 基底行列
		std::vector<std::vector<double>> m_gso_vector;		// GSOベクトル
		std::vector<std::vector<double>> m_gso_coefficient; // GSO係数行列

	public:
		/// @brief 基底行列を出力する関数
		void print()
		{
			for (int i = 0, j; i < m_number_of_row; ++i)
			{
				for (j = 0; j < m_number_of_column; ++j)
				{
					std::cout << m_basis[i][j] << " ";
				}
				std::cout << std::endl;
			}
		}

		/// @brief 格子の階数とかを定義する関数
		/// @param num_rows 基底行列の行数
		/// @param num_cols 基底行列の列数
		void setDimension(const int num_rows, const int num_cols)
		{
			m_number_of_row = num_rows;
			m_number_of_column = num_cols;
			m_basis.resize(num_rows);
			for (int i = 0; i < num_rows; ++i)
			{
				m_basis[i].resize(num_cols);
			}
		}
	};
}
#endif // !LATTICE_H
