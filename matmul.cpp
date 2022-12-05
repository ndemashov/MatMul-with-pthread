#include "matmul.h"
#include <cassert>
#include <chrono>
#include <sstream>

void* MatMul::mul_by_rows(void* _matrices) {
	matrix* matrices = (matrix*)_matrices;
	double result = 0;
	for(unsigned i = 0; i < matrices[0].m; ++i) {
		result += matrices[0].M[0][i] * matrices[1].M[i][0];
	}
	matrices[2].M[0][0] = result;
	return NULL;
}

void* MatMul::mul_by_columns(void* _matrices) {
	matrix* matrices = (matrix*)_matrices;
    std::cout << "new turn with ";
    matrices[0].print();
    matrices[1].print();
    std::cout << "........." << std::endl;
	for(unsigned i = 0; i < matrices[0].n; ++i) {
        for (unsigned j = 0; j < matrices[1].m; ++j) {
            matrices[2].M[i][j] = matrices[0].M[i][0] * matrices[1].M[0][j];  //*=
        }
	}
	return NULL;
}

void MatMul::calc(const CalcType ct, const matrix& m1, const matrix& m2, matrix& result) {
	assert(m1.m == m2.n);
	assert(m1.n == result.n && m2.m == result.m);
	if(ct == CalcType::ByRows){
		unsigned thread_count = m1.n * m2.m;
		pthread_t* thread_handles = (pthread_t*)malloc(thread_count * sizeof(pthread_t));
		matrix** submatrices = (matrix**)malloc(thread_count * sizeof(matrix*));
		for (unsigned row = 0; row < m1.n; ++row) {
			for(unsigned column = 0; column < m2.m; ++column){
				unsigned thread = row * m2.m + column;
				submatrices[thread] = (matrix*)malloc(3 * sizeof(matrix));
				submatrices[thread][0] = matrix(m1, row, 0, MatrixParam::row);
				submatrices[thread][1] = matrix(m2, 0, column, MatrixParam::column);
				submatrices[thread][2] = matrix(result, row, column, MatrixParam::cell);
				pthread_create(&thread_handles[thread], NULL, MatMul::mul_by_rows, (void*)submatrices[thread]);
			}
		}
		for (unsigned row = 0; row < m1.n; ++row) {
			for(unsigned column = 0; column < m2.m; ++column) {
				unsigned thread = row * m2.m + column;
				pthread_join(thread_handles[thread], NULL);
				result.M[row][column] = submatrices[thread][2].M[0][0];
				delete submatrices[thread];
			}
		}
		delete submatrices;
		delete thread_handles;
	} else if(ct == CalcType::ByColumns){
        for(unsigned i = 0; i < result.n; ++i) {
            for (unsigned j = 0; j < result.m; ++j) {
                result.M[i][j] = 0;
            }
	    }
        unsigned thread_count = m1.m;   //по количеству столбцов
		pthread_t* thread_handles = (pthread_t*)malloc(thread_count * sizeof(pthread_t));
        matrix** submatrices = (matrix**)malloc(thread_count * sizeof(matrix*));
		for (unsigned thread = 0; thread < thread_count; ++thread) {
            submatrices[thread] = (matrix*)malloc(3 * sizeof(matrix));
			submatrices[thread][0] = matrix(m1, 0, thread,  MatrixParam::column);
			submatrices[thread][1] = matrix(m2, thread, 0, MatrixParam::row);
			submatrices[thread][2] = matrix(result.n, result.m);
			pthread_create(&thread_handles[thread], NULL, MatMul::mul_by_columns, (void*)submatrices[thread]);
		}
		for (unsigned thread = 0; thread < thread_count; ++thread) {
			pthread_join(thread_handles[thread], NULL);
            for (unsigned row = 0; row < m1.n; ++row) {
			    for(unsigned column = 0; column < m2.m; ++column){
                    result.M[row][column] += submatrices[thread][2].M[row][column];
                }
            }
			delete submatrices[thread];
		}

    } 
}

matrix::matrix(const unsigned _n, const unsigned _m) : n(_n), m(_m) {
	M = (double**)malloc(n * sizeof(double*));
	int step = (int)(n * m) / -2;
	for (unsigned i = 0; i < n; ++i) {
		M[i] = (double*)malloc(m * sizeof(double));
		for (unsigned j = 0; j < m; ++j) {
			//M[i][j] = rand();
			M[i][j] = step;
			step++;
		}
	}
}

matrix::matrix(const matrix& mtrx, const unsigned row, const unsigned column, const MatrixParam mp) {
	if(mp == MatrixParam::row){
		n = 1;
		m = mtrx.m;
		M = (double**)malloc(n * sizeof(double*));
		M[0] = (double*)malloc(m * sizeof(double));
		for(unsigned i = 0; i < m; ++i){
			M[0][i] = mtrx.M[row][i];
		}
	}else if(mp == MatrixParam::column){
		n = mtrx.n;
		m = 1;
		M = (double**)malloc(n * sizeof(double*));
		for(unsigned i = 0; i < n; ++i){
			M[i] = (double*)malloc(m * sizeof(double));
			M[i][0] = mtrx.M[i][column];
		}
	}else if(mp == MatrixParam::cell){
		n = 1;
		m = 1;
		M = (double**)malloc(n * sizeof(double*));
		M[0] = (double*)malloc(m * sizeof(double));
		M[0][0] = mtrx.M[row][column];
	}
}

matrix& matrix::operator=(const matrix& mtrx) {
	n = mtrx.n;
	m = mtrx.m;
	M = (double**)malloc(n * sizeof(double*));
	for(unsigned i = 0; i < n; ++i){
		M[i] = (double*)malloc(m * sizeof(double));
		for(unsigned j = 0; j < m; ++j){
			M[i][j] = mtrx.M[i][j];
		}
	}
	return *this;
}

matrix::~matrix(){
	for (unsigned i = 0; i < n; ++i) {
		delete M[i];
	}
	delete M;
}

Metric::Metric(const std::string metric_file, const unsigned _iter_num) : iter_num(_iter_num) {
	file.open(metric_file);
}

Metric::~Metric() {
	file.close();
}

void Metric::eval(){
	for(unsigned iter = 0; iter < iter_num; ++iter) {
		unsigned n = 5 * (iter + 1);
		matrix m1(n, n), m2(n, 1), result(n, 1);
		std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
		MatMul::calc(CalcType::ByRows, m1, m2, result);
		std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
		file << n << ", " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << std::endl;
	}
}
