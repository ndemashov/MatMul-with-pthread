#include "matmul.h"
#include <string>

#include <thread>



#include "matmul.h"
#include <string>

int test(const CalcType ct){
    matrix m1(5, 2), m2(2, 3);
    double reference[5][3] = {
		{15, 6, -3},
		{9, 4, -1},
		{3, 2, 1},
		{-3, 0, 3},
		{-9, -2, 5},
	};
	matrix result(5, 3);
	MatMul::calc(ct, m1, m2, result);
	for(unsigned i = 0; i < 5; ++i){
		for(unsigned j = 0; j < 3; ++j){
			if(result.M[i][j] != reference[i][j]){
				return 1;
			}
		}
	}
    return 0;
}

int main() {
	// matrix m1(500, 500), m2(500, 500);

	// std::cout << "M1" << std::endl;
	// m1.print();
	// std::cout << "M2" << std::endl;
	// m2.print();
	// matrix result(500, 500);
	// std::cout << "result" << std::endl;
	// result.print();
	// MatMul::calc(CalcType::ByBlocks, m1, m2, result);
	// std::cout << "Result" << std::endl;
	// result.print();

	if(test(CalcType::Single) == 1){
		std::cout << "ERROR" << std::endl;
	}else{
		std::cout <<"Correct" << std::endl;
	}
	Metric M("output.txt", 20);
	M.eval();
	
	return 0;
}
