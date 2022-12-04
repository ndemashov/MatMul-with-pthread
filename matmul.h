#pragma once
#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>
#include <fstream>

#include<iostream>

enum CalcType{
	ByColumns,
	ByRows,
	ByBlocks
};

enum MatrixParam{
	row,
	column
};

struct matrix {
	unsigned n, m;
	double **M;
	matrix(const unsigned _n, const unsigned _m);
	matrix(double** _M, unsigned _n, unsigned _m);
	matrix(const matrix& mtrx, const unsigned ind, const MatrixParam mp);
	matrix(const matrix& mtrx, const unsigned row_ind, const unsigned col_ind, 
			const unsigned amount_elements_by_row, const unsigned amount_elements_by_col);
	~matrix();
	//matrix(const matrix&) = delete;
	matrix& operator=(const matrix& mtrx);
	
	// temp
	void print() {
		for (unsigned i = 0; i < n; ++i) {
			for (unsigned j = 0; j < m; ++j) {
				std::cout << M[i][j] << ' ';
			}
			std::cout << std::endl;
		}
	}
};

class MatMul {
public:
	static void* mul_by_columns(void* args);
	static void* mul_by_rows(void* args);
	static void* mul_by_blocks(void* args); 
	static void calc(const CalcType ct, const matrix& m1, const matrix& m2, matrix& result);
};

class Decomposition {
public:
	static unsigned amount_threads(const matrix& mtrx, unsigned& row_block, unsigned& col_block);  
};

class Metric{
	std::ofstream file;
	unsigned iter_num;
public:
	Metric(const std::string metric_file, const unsigned _iter_num);
	~Metric();
	void eval();
};
