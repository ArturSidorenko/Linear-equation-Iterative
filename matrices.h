#pragma once

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<algorithm>
#include<utility>
#include<exception>
#include<iomanip>
#include<cmath>
#include<vector>


const double PI = 4 * atan(1);
typedef std::vector<double> method_tau;
typedef std::vector<double> sol_logs; //logs of system solving
typedef std::vector<size_t> permut;

/*
Custom matrix class
Unlike common approaches, this matrix is stored column-wise, i.e.
the first index is in charge of columns, not rows.
It allows to represent vector-columns in the same class.
*/
class matrix {
	size_t rows_, cols_; //sizes
	double **data_; //contents
public:
	//constructors
	
	matrix(size_t rows, size_t cols); //zero matrix constructor
	matrix(const matrix &other); //copy construnctor
	matrix(matrix &&other); //move constructor
	~matrix();
	

	//operators =
	matrix operator=(const matrix &other);
	matrix operator=(matrix &&other);

	//pluses and multiplies
	//introduce convenient syntax
	matrix operator+(const matrix &M) const;
	matrix operator-(const matrix &M) const;
	matrix operator*(const matrix &M) const;
	matrix operator*(double x) const; //one can write either A*x,...
	friend matrix operator*(double x, const matrix &A);  //...or x*A


	//norms
	double l1_norm() const;
	double linf_norm() const;  

	//access functions
	double get(size_t i, size_t j) const;
	void set(size_t i, size_t j, double x);
	void print(std::ofstream &s) const; //prints matrix to file
	void print(const std::string &s) const; //does the same
	size_t rows() const { return rows_; };
	size_t cols() const { return cols_; };
};

//Identity matrix
matrix id(size_t n);
//Lagrange matrix constructor
matrix lagrange(size_t N); 

//Chebyshev grid
double cheb_grid(unsigned ord, unsigned i, double a, double b);
double inv_cheb_grid(unsigned ord, unsigned k, double a, double b);

//some n-steps method
//it performs sonsecutive multiplications like (I-t*A)*X
//method_tau - array of parameters t
//X is either a vector-column or a set columns in the one matrix
//A must be a square matrix and the number of rows of X must coincide with the size of A
matrix n_steps_method(const matrix &A, const matrix &X, const matrix &B, const method_tau &t);

//copy of the previous function to test convergence
matrix test_n_steps_method(const matrix & A, const matrix & X, const matrix &B, const method_tau & t, const matrix &true_ans);

bool break_loop(const sol_logs &g);

//very clever permutation of length = 2^t
permut clever_met(size_t t);



