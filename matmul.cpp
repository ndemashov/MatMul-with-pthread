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

	for (unsigned i = 0; i < matrices[2].n; ++i) {
		for (unsigned j = 0; j < matrices[2].m; ++j) {
			matrices[2].M[i][j] = 0.0;
		}
	}

	for (unsigned i = 0; i < matrices[2].n; ++i) {
		for (unsigned k = 0; k < matrices[0].m; ++k) {
			for (unsigned j = 0; j < matrices[2].m; ++j) {
				matrices[2].M[i][j] += matrices[0].M[i][k] * matrices[1].M[k][j];
			}
		}
	}
}


void* MatMul::mul_by_columns(void* _matrices) {
	matrix* matrices = (matrix*)_matrices;
    //std::cout << "new turn with ";
    //matrices[0].print();
    //matrices[1].print();
    //std::cout << "........." << std::endl;
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
	if(ct == CalcType::ByRows) {
		unsigned thread_count = m1.n * m2.n;
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

	}

	if(ct == CalcType::ByBlocks) {

		for (unsigned row_ind = 0; row_ind < result.n; ++row_ind) {
			for (unsigned col_ind = 0; col_ind < result.m; ++col_ind) {
				result.M[row_ind][col_ind] = 0;
			}
		}

		for (unsigned col = 0; col < m2.m; ++col) {

			matrix m2_(m2.n, 1);

			for (unsigned i = 0; i < m2.n; ++i) {
				m2_.M[i][0] = m2.M[i][col];
			}
			//std::cout << "m2_" << std::endl;
			//m2_.print();
			//std::cout << "!!!COL = " << col << std::endl;

			unsigned row_block_m1 = 1;
			unsigned col_block_m1 = 1;
			unsigned row_block_m2 = 1;
			unsigned col_block_m2 = 1;

			// Сюда лучше передавать матрицу с наибольшей размерностью, но пока это m1
			unsigned thread_count = Decomposition::amount_threads(m1, m2_, row_block_m1, col_block_m1, row_block_m2, col_block_m2);
			//unsigned thread_count_m2 = Decomposition::amount_threads(m2, row_block_m2, col_block_m2);

			// std::cout<<"m1 row_block = " << row_block_m1 << std::endl;
			// std::cout<<"m1 col_block = " << col_block_m1 << std::endl;

			// std::cout<<"m2 row_block = " << row_block_m2 << std::endl;
			// std::cout<<"m2 col_block = " << col_block_m2 << std::endl;

			//unsigned thread_count = std::max(thread_count_m1, thread_count_m2);
			//std::cout << "thread_count = " << thread_count << std::endl;
			// row_block * row_ind = thread_count
			pthread_t* thread_handles = (pthread_t*)malloc(thread_count * sizeof(pthread_t));
			matrix** submatrices = (matrix**)malloc(thread_count * sizeof(matrix*));
			unsigned thread = 0;
			unsigned max_amount_elements_by_row_m1 = ceil(m1.n / (double)row_block_m1);
			//std::cout << "max_amount_elements_by_row_m1 = " << max_amount_elements_by_row_m1 << std::endl;
			unsigned max_amount_elements_by_col_m1 = ceil(m1.m / (double)col_block_m1);
			//std::cout << "max_amount_elements_by_col_m1 = " << max_amount_elements_by_col_m1 << std::endl;
			unsigned max_amount_elements_by_row_m2 = ceil(m2_.n / (double)row_block_m2);
			//std::cout << "max_amount_elements_by_row_m2 = " << max_amount_elements_by_row_m2 << std::endl;
			unsigned max_amount_elements_by_col_m2 = ceil(m2_.m / (double)col_block_m2);
			//std::cout << "max_amount_elements_by_col_m2 = " << max_amount_elements_by_col_m2 << std::endl;

			for (unsigned thread = 0; thread < thread_count; ++thread) {
				submatrices[thread] = (matrix*)malloc(3 * sizeof(matrix));
			}

			unsigned amount_elements_by_row_m1_ = max_amount_elements_by_row_m1;
			unsigned amount_elements_by_col_m1_ = max_amount_elements_by_col_m1;
			unsigned row_block_m1_ = row_block_m1;
			unsigned col_block_m1_ = col_block_m1;
			unsigned m_m1 = m1.m;
			unsigned n_m1 = m1.n;
			thread = 0;

			//std::cout << "Create submatrices[0]\n";
			while (thread != thread_count) {
				for (unsigned row_ind = 0; row_ind < row_block_m1; ++row_ind) {
					for (unsigned col_ind = 0; col_ind < col_block_m1; ++col_ind) {
						// std::cout << "row_ind = " << row_ind << std::endl;
						// std::cout << "col_ind = " << col_ind << std::endl;
						submatrices[thread][0] = matrix(m1, row_ind, col_ind, amount_elements_by_row_m1_, amount_elements_by_col_m1_,
														m1.n-n_m1, m1.m - m_m1);

						// std::cout << "Submatrices[0]" << std::endl;
						// submatrices[thread][0].print();
						
						m_m1 = m_m1 - amount_elements_by_col_m1_;
						col_block_m1_ = col_block_m1_ - 1;
						if (col_block_m1_ != 0) {
							amount_elements_by_col_m1_ = ceil(m_m1 / double(col_block_m1_));
						} else {
							//std::cout << "col_block_m1_ = 0" << std::endl;
						}

						thread++;
					}
					// Фикс столбцов
					m_m1 = m1.m;
					col_block_m1_ = col_block_m1;
					amount_elements_by_col_m1_ = max_amount_elements_by_col_m1;
					//Фикс строк
					n_m1 = n_m1 - amount_elements_by_row_m1_;
					row_block_m1_ = row_block_m1_ - 1;
					if (row_block_m1 != 0) {
						amount_elements_by_row_m1_ = ceil(n_m1 / double(row_block_m1_));
					} else {
						//std::cout << "row_block_m1 = 0" << std::endl;
					}
				}
				// Фикс строк
				n_m1 = m1.n;
				row_block_m1_ = row_block_m1;
				amount_elements_by_row_m1_ = max_amount_elements_by_row_m1;
			}

			unsigned amount_elements_by_row_m2_ = max_amount_elements_by_row_m2;
			unsigned amount_elements_by_col_m2_ = max_amount_elements_by_col_m2;
			unsigned row_block_m2_ = row_block_m2;
			unsigned col_block_m2_ = col_block_m2;
			unsigned m_m2 = m2_.m;
			unsigned n_m2 = m2_.n;
			thread = 0;

			//std::cout << "Create submatrices[1]\n";
			while (thread != thread_count) {
				for (unsigned row_ind = 0; row_ind < row_block_m2; ++row_ind) {
					for (unsigned col_ind = 0; col_ind < col_block_m2; ++col_ind) {
						// std::cout << "row_ind = " << row_ind << std::endl;
						// std::cout << "col_ind = " << col_ind << std::endl;
						submatrices[thread][1] = matrix(m2_, row_ind, col_ind, amount_elements_by_row_m2_, amount_elements_by_col_m2_,
														m2_.n - n_m2, m2_.m - m_m2);

						// std::cout << "Submatrices[1]" << std::endl;
						// submatrices[thread][1].print();

						m_m2 = m_m2 - amount_elements_by_col_m2_;
						col_block_m2_ = col_block_m2_ - 1;
						if (col_block_m2_ != 0) {
							amount_elements_by_col_m2_ = ceil(m_m2 / double(col_block_m2_));
						} else {
							//std::cout << "col_block_m2_ = 0" << std::endl;
						}
						// std::cout << "m_m2=" << m_m2 << std::endl;
						// std::cout << "col_block_m2_=" << col_block_m2_ << std::endl;
						// std::cout << "amount_elements_by_col_m2_=" << amount_elements_by_col_m2_ << std::endl;
						thread++;
					}
					// Фикс столбцов
					m_m2 = m2_.m;
					col_block_m2_ = col_block_m2;
					amount_elements_by_col_m2_ = max_amount_elements_by_col_m2;

					//Фикс строк
					n_m2 = n_m2 - amount_elements_by_row_m2_;
					row_block_m2_ = row_block_m2_ - 1;
					if (row_block_m2_ != 0) {
						amount_elements_by_row_m2_ = ceil(n_m2 / double(row_block_m2_));
					} else {
						//std::cout << "row_block_m2_ = 0" << std::endl;
					}
				}
				// Фикс строк
				n_m2 = m2_.n;
				row_block_m2_ = row_block_m2;
				amount_elements_by_row_m2_ = max_amount_elements_by_row_m2;
			}

			unsigned row_block_m3 = std::max(row_block_m1, row_block_m2);
			unsigned col_block_m3 = std::max(col_block_m1, col_block_m2);
			unsigned row_block_m3_ = row_block_m3;
			unsigned col_block_m3_ = col_block_m3;
			thread = 0;

			//std::cout << "Create submatrices[2]\n";
			while (thread != thread_count) {
				for (unsigned row_ind = 0; row_ind < row_block_m3; ++row_ind) {
					for (unsigned col_ind = 0; col_ind < col_block_m3; ++col_ind) {
						// std::cout << "row_ind = " << row_ind << std::endl;
						// std::cout << "col_ind = " << col_ind << std::endl;

						submatrices[thread][2] = matrix(result, row_ind, col_ind, submatrices[thread][0].n, submatrices[thread][1].m,
														0, 0);

						// std::cout << "Submatrices[2]" << std::endl;
						// submatrices[thread][2].print();

						thread++;
					}
				}
			}

			thread = 0;
			while (thread != thread_count) {
				pthread_create(&thread_handles[thread], NULL, MatMul::mul_by_blocks, (void*)submatrices[thread]);
				thread++;
			}

			thread = 0;
			while (thread != thread_count) {
				pthread_join(thread_handles[thread], NULL);
				thread++;
			}

			// thread = 0;
			// while (thread != thread_count) {
			// 	std::cout << "Matrix after thread" << std::endl;
			// 	submatrices[thread][2].print();
			// 	thread++;
			// }

			// std::cout << "Print after 0" << std::endl;
			// result.print();

			unsigned v = 0;
			unsigned u = 0;
			thread = 0;
			unsigned sum = 0;
			for (unsigned row_ind = 0; row_ind < row_block_m3; ++row_ind) {
				for (unsigned col_ind = 0; col_ind < col_block_m3; ++col_ind) {
					for (unsigned i = 0, k = v; i < submatrices[thread][2].n; ++i, ++k) {
						for (unsigned j = 0; j < submatrices[thread][2].m; ++j) {
							//std::cout << "k = " << k << std::endl;
							result.M[k][col] += submatrices[thread][2].M[i][j];
							//std::cout << result.M[k][col] << std::endl;
						}
					}
					thread++;
				}
				v += submatrices[thread - 1][2].n;
				v = v >= result.n ? 0 : v;
			}
	
			for (unsigned thread = 0; thread < thread_count; ++thread) {
				delete submatrices[thread];
			}
			delete submatrices;
			delete thread_handles;

			// std::cout << "++++result++++" << std::endl;
			// result.print();
		}
	} 

	if(ct == CalcType::ByColumns) {
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

	if(ct == CalcType::Single){
		for(int i = 0; i < m1.n; ++i){
			for(int j = 0; j < m2.m; ++j){
				double res = 0;
				for(int z = 0; z < m1.m; ++z){
					res += m1.M[i][z]*m2.M[z][j];
				}
				result.M[i][j] = res;
			}
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

matrix::matrix(const matrix& mtrx, const unsigned row_ind, const unsigned col_ind, 
				const unsigned amount_elements_by_row, const unsigned amount_elements_by_col,
				const unsigned ost_row, const unsigned ost_col) { // ost = m1.n - n_m1
	n = unsigned(amount_elements_by_row);
	m = unsigned(amount_elements_by_col);
	M = (double**)malloc(n * sizeof(double*));

	// std::cout <<"amount_elements_by_row = " << n << std::endl;
	// std::cout <<"amount_elements_by_col = " << m << std::endl;
	// std::cout << "ost_row = " << ost_row << std::endl;
	// std::cout << "ost_col = " << ost_col << std::endl;



	for (unsigned i = 0; i < n; ++i) {
		M[i] = (double*)malloc(m * sizeof(double));
	}

	unsigned v;
	unsigned u;
	for (unsigned i = 0; i < n; ++i) {
		v = ost_row + i;
		//std::cout << "v = " << v << std::endl;
		for (unsigned j = 0; j < m; ++j) {
			u = ost_col + j;
			//std::cout << "u = " << u << std::endl;
			M[i][j] = mtrx.M[v][u];	
		}
	}
	
	//std::cout << "n = " << n << std::endl;
	//std::cout << "m = " << m << std::endl;
	//for (unsigned i = 0; i < n; ++i) {
		//std::cout << "Я тут" << std::endl;
	//	for (unsigned j = 0; j < m; ++j) {
	//		//std::cout<<"Уже ближе!" << std::endl;
	//		std::cout << M[i][j] << std::endl;
	//		//std::cout << "Ура" << std::endl;
	//	}
	//}

	
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

unsigned Decomposition::amount_threads(const matrix& mtrx1, const matrix& mtrx2, unsigned& row_block1, unsigned& col_block1,
									unsigned& row_block2, unsigned& col_block2) {
	// Основано на том, что общее количество процессов представимо в степени 2
	// Определение количества процесcоров
	const auto processor_count = std::thread::hardware_concurrency();
	// Разобьем задачу на поиск блоков по строкам и столбцам. Изначально полагаем, что 
	// количество блоков по строкам равно общему количеству процессов, а по столбцам - 1.
	row_block1 = processor_count;
	col_block1 = 1;
	// Для разбиения на блоки по строкам необходимо, чтобы mtrx.n была больше, чем row_block
	while (mtrx1.n < row_block1) {
		row_block1 = row_block1 >> 1;
	}
	// "Остаток" процессоров записываем в col_block так, чтобы col_block * row_block = process_count
	//std::cout << "rov_block >> 1 = " << (row_block >> 1) << std::endl;
	//std::cout << "processor_count >> (row_block >> 1) = " << (processor_count >> (row_block >> 1)) << std::endl;
	col_block1 = processor_count >> (row_block1 >> 1);
	//std::cout << "col_block in decomposition = " << col_block << std::endl;
	// Для разбиения на блоки по столбцам необходимо, чтобы mtrx.m была больше, чем col_block
	//std::cout << "mtrx.m = " << mtrx.m << std::endl;
	col_block1 = mtrx1.m > col_block1 ? col_block1 : mtrx1.m;
	//std::cout << "col_block ending = " << col_block << std::endl;

	row_block2 = col_block1;
	col_block2 = processor_count >> (row_block2 >> 1);
	col_block2 = mtrx2.m > col_block2 ? col_block2 : mtrx2.m;
	// Произедение col_block * row_block <= processor_count
	return col_block1 * row_block1;
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
		matrix m1(n, n), m2(n, n), result(n, n);
		std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
		MatMul::calc(CalcType::ByColumns, m1, m2, result);
		std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
		file << n << ", " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << std::endl;
	}
}
