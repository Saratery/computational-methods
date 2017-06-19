#ifndef METHODS_H_INCLUDED
#define METHODS_H_INCLUDED

#include <iostream>
#include <stdio.h>
#include "Matrix.h"
#include <vector>
#include <fstream>
#include <string>


Matrix<double> Gauss_full_choice(Matrix<double> A, Matrix<double> b);

double norm_2(const Matrix<double>& x);

Matrix<double> Householder(Matrix<double> A, Matrix<double> b);



#endif // METHODS_H_INCLUDED
