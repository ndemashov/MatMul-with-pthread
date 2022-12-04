#include "matmul.h"
#include <cassert>

matrix T(matrix m1){
        matrix t(m1.m, m1.n);
        for (unsigned i = 0; i < m1.n; ++i) {
			for (unsigned j = 0; j < m1.m; ++j) {
                t.M[j][i] = m1.M[i][j];
            }
        }
        return t;
    }
void* MatMul::mul_by_rows(void* _matrices) {
	matrix* matrices = (matrix*)_matrices;
	double result = 0;
	for(unsigned i = 0; i < matrices[1].m; ++i) {
        for (unsigned j = 0; j < matrices[0].m; ++j) {
		    result += matrices[0].M[0][j] * matrices[1].M[j][i];
        }
        matrices[2].M[0][i] = result;
        result = 0;
	}
	return NULL;
}

void* MatMul::mul_by_columns(void* _matrices) {
	matrix* matrices = (matrix*)_matrices;
    std::cout << "new turn with ";
    matrices[0].print();
    matrices[1].print();
    std::cout << "........." << std::endl;
	for(unsigned i = 0; i < matrices[0].m; ++i) {
        for (unsigned j = 0; j < matrices[1].m; ++j) {
		    matrices[3].M[i][j] = matrices[0].M[0][i] * matrices[1].M[0][j];
            
            matrices[2].M[i][j] += matrices[0].M[0][i] * matrices[1].M[0][j];
        }
	}
	return NULL;
}

matrix MatMul::calc(const CalcType ct, matrix& m1, matrix& m2) {
	assert(m1.m == m2.n);
    if(ct == CalcType::ByRows){
        matrix* matrices = (matrix*)malloc(3 * sizeof(matrix));
        matrices[0] = m1;
        matrices[1] = m2;
        matrix result(m1.n, m2.m);
        matrices[2] = result;
		unsigned thread_count = m1.n;   //по количеству строк
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
        return matrices[2];
	} else if(ct == CalcType::ByColumns){
        matrix* matrices = (matrix*)malloc(4 * sizeof(matrix));
        matrices[0] = T(m1);
        matrices[1] = m2;
        matrix result(m1.n, m2.m);
        for(unsigned i = 0; i < result.n; ++i) {
            for (unsigned j = 0; j < result.m; ++j) {
                result.M[i][j] = 0;
            }
	    }
        matrix preresult(m1.n, m2.m);
        matrices[2] = result;
        matrices[3] = preresult;
        unsigned thread_count = m1.m;   //по количеству столбцов
		pthread_t* thread_handles = (pthread_t*)malloc(thread_count * sizeof(pthread_t));
		for (unsigned thread = 0; thread < thread_count; ++thread) {
			matrix* submatrices = (matrix*)malloc(4 * sizeof(matrix));
			submatrices[0] = matrix(&matrices[0].M[thread], 1, matrices[0].m);
			submatrices[1] = matrix(&matrices[1].M[thread], 1, matrices[1].m);
			submatrices[2] = matrices[2];
            submatrices[3] = matrices[3];
			pthread_create(&thread_handles[thread], NULL, MatMul::mul_by_columns, (void*)submatrices);
		}
		for (unsigned thread = 0; thread < thread_count; ++thread) {
			pthread_join(thread_handles[thread], NULL);
		}
        return matrices[2];
    } 
    return m1;
	
}