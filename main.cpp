#include <iostream>
#include <time.h>
#include "Linalg/Core"

int main()
{
	Linalg::Vector<int> vec1, vec2, vec3;
	Linalg::Matrix<float> mat1, mat2;

	puts("Vector");

	// sets length of vector
	vec1.setLength(3);
	vec2.setLength(3);
	vec3.setLength(4);

	vec1.setSeed();
	vec1.random(-10, 10);
	vec2.random(-20, 20);

	// prints vector
	printf("vec1 = ");
	vec1.print();
	printf("vec2 = ");
	vec2.print();

	// sum
	puts("\nsum");
	vec1 += vec2;
	vec1.print();

	// difference
	puts("\ndiffer");
	vec3 = vec1 - vec2;
	vec3.print();

	// inner product
	puts("\ninner product");
	std::cout << vec1.innerProduct(vec2) << std::endl;

	puts("\nMatrix");

	// sets size of matrix
	mat1.setIdentity(3);
	mat2.setDimension(3, 3);
	mat1.substituteRow(0, vec1.cast<float>());
	mat2.random(-10, 10);
	mat1.print();

	puts("sum");
	(mat1 + mat2).print();

	// determinant
	puts("\ndeterminant:");
	std::cout << mat1.det() << std::endl;

	// transpose
	puts("\ntranspose:");
	mat1.transpose().print();

	// co-Factor
	puts("\nadjugate:");
	mat1.adjugate().print();
	puts("its test");
	(mat1 * mat1.adjugate()).print();

	// inverse
	puts("\ninverse:");
	mat1.inverse().print();
	puts("its test");
	(mat1.cast<double>() * mat1.inverse()).print();

	// power
	puts("\nqubic");
	mat1.power(3).print();
	return 0;
}
