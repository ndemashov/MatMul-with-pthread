#include "matmul.h"


int main() {
	matrix m1(5, 2), m2(2, 1);
	//std::cout << "M1" << std::endl;
	//m1.print();
	//std::cout << "M2" << std::endl;
	//m2.print();
	matrix* matrices = (matrix*)malloc(2 * sizeof(matrix));
	if (matrices == nullptr) {
		std::cout << "not alloc";
	}
	matrix result(10, 1); 
	result = MatMul::calc(m1, m2);
	std::cout << "Result" << std::endl;
	result.print();
	return 0;
}
