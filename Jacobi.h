#ifndef JACOBI_H_INCLUDED
#define JACOBI_H_INCLUDED
#include <fstream>
#include <iostream>
#include "Matrix.h"
#include "math.h"

double Frob_norm(const Matrix<double> &A);

const Matrix<double> G(const double& a_ii, const double& a_ij, const double& a_jj, const int& n, const int& i, const int& j);

Matrix<double> JacobiMethod(const Matrix<double>& M, const double &E);

#endif // JACOBI_H_INCLUDED