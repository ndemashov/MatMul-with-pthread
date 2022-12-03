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

void MatMul::calc(const CalcType ct, const matrix& m1, const matrix& m2, matrix& result) {
	assert(m1.m == m2.n);
	assert(m1.n == result.n && m2.m == result.m);
	if(ct == CalcType::ByRows){
		unsigned thread_count = m1.n;
		pthread_t* thread_handles = (pthread_t*)malloc(thread_count * sizeof(pthread_t));
		matrix** submatrices = (matrix**)malloc(thread_count * sizeof(matrix*));
		for (unsigned thread = 0; thread < thread_count; ++thread) {
			submatrices[thread] = (matrix*)malloc(3 * sizeof(matrix));
			submatrices[thread][0] = matrix(m1, thread, MatrixParam::row);
			submatrices[thread][1] = m2;
			submatrices[thread][2] = matrix(result, thread, MatrixParam::row);
			pthread_create(&thread_handles[thread], NULL, MatMul::mul_by_rows, (void*)submatrices[thread]);
		}
		for (unsigned thread = 0; thread < thread_count; ++thread) {
			pthread_join(thread_handles[thread], NULL);
			result.M[thread][0] = submatrices[thread][2].M[0][0];
			delete submatrices[thread];
		}
		delete submatrices;
		delete thread_handles;
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

matrix::matrix(const matrix& mtrx, const unsigned ind, const MatrixParam mp) {
	if(mp == MatrixParam::row){
		n = 1;
		m = mtrx.m;
		M = (double**)malloc(n * sizeof(double*));
		M[0] = (double*)malloc(m * sizeof(double));
		for(unsigned i = 0; i < m; ++i){
			M[0][i] = mtrx.M[ind][i];
		}
	}else{
		n = mtrx.n;
		m = 1;
		M = (double**)malloc(n * sizeof(double*));
		for(unsigned i = 0; i < m; ++i){
			M[i] = (double*)malloc(m * sizeof(double));
			M[i][0] = mtrx.M[ind][i];
		}
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
