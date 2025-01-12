#include <iostream>
#include "Linalg/Core"

int main()
{
	Linalg::Vector<int> vec1, vec2, vec3;
	Linalg::Matrix<long> mat1;

	puts("Vector");

	// sets length of vector
	vec1.setLength(3);
	vec2.setLength(3);
	vec3.setLength(4);

	// substitution
	for (int i = 0; i < 3; ++i)
	{
		vec1.substitute(i, i + 1);
		vec2.substitute(i, i + i + 2);
	}

	// prints vector
	printf("vec1 = ");
	vec1.print();
	printf("vec2 = ");
	vec2.print();

	// sum
	puts("sum");
	vec1 += vec2;
	vec1.print();

	// d
	puts("differ");
	vec3 = vec1 - vec2;
	vec3.print();

	puts("Matrix");

	// sets size of matrix
	mat1.setIdentity(3);
	mat1.substituteRow(0, vec1.cast<long>());
	mat1.print();

	// determinant
	puts("determinant:");
	std::cout << mat1.det() << std::endl;

	// transpose
	puts("transpose:");
	mat1.transpose().print();

	// co-Factor
	puts("co-factor:");
	mat1.coFactor().print();
	puts("its test");
	(mat1 * mat1.coFactor()).print();

	// inverse
	puts("inverse:");
	mat1.inverse().print();
	puts("its test");
	(mat1.cast<double>() * mat1.inverse()).print();
	return 0;
}
