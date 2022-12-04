#include "matmul.h"
#include <cassert>
#include <chrono>
#include <sstream>
#include <thread> //для определения количества процессоров
#include <cmath>

void* MatMul::mul_by_rows(void* _matrices) {
	matrix* matrices = (matrix*)_matrices;
	double result = 0;
	for(unsigned i = 0; i < matrices[0].m; ++i) {
		result += matrices[0].M[0][i] * matrices[1].M[i][0];
	}
	matrices[2].M[0][0] = result;
	return NULL;
}

void* MatMul::mul_by_blocks(void* _matrices) {
	matrix* matrices = (matrix*)_matrices;
	for(unsigned i = 0; i < matrices[0].n; ++i) {
		for (unsigned j = 0; j < matrices[1].m; ++j) {
			for (unsigned k = 0; k < matrices[0].m; ++k) {
				matrices[2].M[i][j] += matrices[0].M[i][k] * matrices[1].M[k][j];;
			}
		}
	}
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

	if(ct == CalcType::ByBlocks) {
		unsigned row_block_m1 = 1;
		unsigned col_block_m1 = 1;
		unsigned row_block_m2 = 1;
		unsigned col_block_m2 = 1;
		// Сюда лучше передавать матрицу с наибольшей размерностью, но пока это m1
		unsigned thread_count_m1 = Decomposition::amount_threads(m1, row_block_m1, col_block_m1);
		unsigned thread_count_m2 = Decomposition::amount_threads(m1, row_block_m2, col_block_m2);
		unsigned thread_count = std::max(thread_count_m1, thread_count_m2);
		// row_block * row_ind = thread_count
		pthread_t* thread_handles = (pthread_t*)malloc(thread_count * sizeof(pthread_t));
		matrix** submatrices = (matrix**)malloc(thread_count * sizeof(matrix*));
		unsigned thread = 0;
		unsigned amount_elements_by_row_m1 = ceil(m1.n / row_block_m1);
		unsigned amount_elements_by_col_m1 = ceil(m1.m / col_block_m1);
		unsigned amount_elements_by_row_m2 = ceil(m2.n / row_block_m2);
		unsigned amount_elements_by_col_m2 = ceil(m2.m / col_block_m2);

		for (unsigned row_ind_m1 = 0, row_ind_m2 = 0, i = 0; i < std::max(row_block_m1, row_block_m2) ; ++i) {
			for (unsigned col_ind_m1 = 0, col_ind_m2 = 0, j = 0; j < std::max(col_block_m1, col_block_m2); ++j) {
				if (col_ind_m1 == col_block_m1) {
					row_ind_m1++;
					col_ind_m1 = 0;
				}
				if (col_ind_m2 == col_block_m2) {
					row_ind_m2++;
					col_ind_m2 = 0;
				}
				submatrices[thread] = (matrix*)malloc(3 * sizeof(matrix));
				submatrices[thread][0] = matrix(m1, row_ind_m1, col_ind_m1, amount_elements_by_row_m1, amount_elements_by_col_m1);
				submatrices[thread][1] = matrix(m2, row_ind_m2, col_ind_m2, amount_elements_by_row_m2, amount_elements_by_col_m2);
				submatrices[thread][2] = matrix(result, row_ind_m1, col_ind_m2, amount_elements_by_row_m1, amount_elements_by_col_m2);
				pthread_create(&thread_handles[thread], NULL, MatMul::mul_by_blocks, (void*)submatrices[thread]);
				thread++;
				col_ind_m1++;
				col_ind_m2++;
			}
			row_ind_m1++;
			row_ind_m2++;
			if (row_ind_m1 >= row_block_m1) {
				row_ind_m1 = 0;
			} 
			if (row_ind_m2 >= row_block_m2) {
				row_ind_m2 = 0;
			}
		}



		
		// Не получается сложить блочные матрицы в результирующий массив
		for (unsigned thread = 0; thread < thread_count; ++thread) {
			pthread_join(thread_handles[thread], NULL);
		}
		thread = 0;
		unsigned v = 0;
		unsigned u = 0;

		for (unsigned i = 0; i < std::max(row_block_m1, row_block_m2); ++i) {
			for (unsigned j = 0; j < std::max(col_block_m1, col_block_m2); ++j, ++thread) {
				for (unsigned k = 0; k < submatrices[thread][2].n; ++k) {
					for (unsigned t = 0; t < submatrices[thread][2].m; ++t) {
						result.M[k][t] += submatrices[thread][2].M[k][t];
					}
				}
			}
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

matrix::matrix(const matrix& mtrx, const unsigned row_ind, const unsigned col_ind, 
				const unsigned amount_elements_by_row, const unsigned amount_elements_by_col) {
	n = amount_elements_by_row;
	m = amount_elements_by_col;
	M = (double**)malloc(n * sizeof(double*));

	for (unsigned i = 0, v = row_ind * n + i; i < n; ++i) {
		for (unsigned j = 0, u = col_ind * m + j; j < m; ++j) {
			if (u < mtrx.m && v < mtrx.n) {
				M[i] = (double*)malloc(m * sizeof(double));
				M[i][j] = mtrx.M[v][u];
			}
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

unsigned Decomposition::amount_threads(const matrix& mtrx, unsigned& row_block, unsigned& col_block) {
	// Основано на том, что общее количество процессов представимо в степени 2
	// Определение количества процесcоров
	const auto processor_count = std::thread::hardware_concurrency();
	// Разобьем задачу на поиск блоков по строкам и столбцам. Изначально полагаем, что 
	// количество блоков по строкам равно общему количеству процессов, а по столбцам - 1.
	row_block = processor_count;
	col_block = 1;
	// Для разбиения на блоки по строкам необходимо, чтобы mtrx.n была больше, чем row_block
	while (mtrx.n <= row_block) {
		row_block = row_block >> 1;
	}
	// "Остаток" процессоров записываем в col_block так, чтобы col_block * row_block = process_count
	col_block = processor_count >> (row_block >> 1);
	// Для разбиения на блоки по столбцам необходимо, чтобы mtrx.m была больше, чем col_block
	col_block = col_block ? mtrx.m > col_block : mtrx.m;

	// Произедение col_block * row_block <= processor_count
	return col_block * row_block;
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
