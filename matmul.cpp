#include "matmul.h"
#include <cassert>

void* MatMul::mul_by_rows(void* _matrices) {
	matrix* matrices = (matrix*)_matrices;
	double result = 0;
	for(unsigned i = 0; i < matrices[0].m; ++i) {
		result += matrices[0].M[0][i] * matrices[1].M[i][0];
	}
	matrices[2].M[0][0] = result;
	return NULL;
}

matrix MatMul::calc(const CalcType ct, matrix& m1, matrix& m2) {
	assert(m1.m == m2.n);
	matrix* matrices = (matrix*)malloc(3 * sizeof(matrix));
	matrices[0] = m1;
	matrices[1] = m2;
	matrix result(m1.n, m2.m);
	matrices[2] = result;
	if(ct == CalcType::ByRows){
		unsigned thread_count = m1.n;
		pthread_t* thread_handles = (pthread_t*)malloc(thread_count * sizeof(pthread_t));
		for (unsigned thread = 0; thread < thread_count; ++thread) {
			matrix* submatrices = (matrix*)malloc(3 * sizeof(matrix));
			submatrices[0] = matrix(&matrices[0].M[thread], 1, matrices[0].m);
			submatrices[1] = matrices[1];
			submatrices[2] = matrix(&matrices[2].M[thread], 1, 1);
			pthread_create(&thread_handles[thread], NULL, MatMul::mul_by_rows, (void*)submatrices);
		}
		for (unsigned thread = 0; thread < thread_count; ++thread) {
			pthread_join(thread_handles[thread], NULL);
		}
	}
	return matrices[2];
}