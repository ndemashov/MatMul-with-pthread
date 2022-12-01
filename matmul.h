#pragma once
#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>
#include <vector>
#include <algorithm>
#include <iterator>

#include<iostream>

struct matrix {
	unsigned n, m;
	double **M;
	matrix(const unsigned _n, const unsigned _m) : n(_n), m(_m) {
		M = (double**)malloc(n * sizeof(double*));
		unsigned step = 1;
		for (unsigned i = 0; i < n; ++i) {
			M[i] = (double*)malloc(m * sizeof(double));
			for (unsigned j = 0; j < m; ++j) {
				//M[i][j] = rand();
				M[i][j] = step;
				step++;
			}
		}
	}
	matrix(double** _M, unsigned _n, unsigned _m) : M(_M), n(_n), m(_m) {}

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
	static matrix calc(matrix& m1, matrix& m2);
};
