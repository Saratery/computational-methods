#include <fstream>
#include <iostream>
#include "Matrix.h"
#include "math.h"


double Frob_norm(const Matrix<double> &A)
{
    double s = 0;
    for (int i = 1; i <= A.Get_hei(); i++)
        for (int j = 1; j <= A.Get_hei(); j++)
        {
            if (i != j)
            {
                s += A[i][j]*A[i][j];
            }
        }
        return sqrt(s);
}


const Matrix<double> G(const double& a_ii, const double& a_ij, const double& a_jj, const int& n, const int& i, const int& j)
{
    Matrix<double> I(n, n);
    for (int p = 1; p <= n; p++)
        {
            for (int q = 1; q <= n; q++)
            {
                if (p != q)
                {
                    I[p][q] = 0;
                }
                else
                {
                    I[p][q] = 1;
                }
            }
        }
    double tau = (a_ii - a_jj)/(2*a_ij);
    double t = -tau + sqrt(tau*tau + 1);
    double c = sqrt(1/(t*t + 1));
    double s = sqrt(1 - c*c);
    I[i][i] = c;
    I[j][j] = c;
    I[i][j] = -s;
    I[j][i] = s;
    return I;
}

Matrix<double> JacobiMethod(const Matrix<double>& M, const double &E)
{
    Matrix<double> A = M;
    const int n = A.Get_hei();
    double barrier = 1;
    int k = 0;
    int l = 0;
    while (Frob_norm(A) > E)
    {
        l++;
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
            {
                if ((i != j) && (fabs(A[i][j]) > barrier))
                {
                    const Matrix<double> Q = G(A[i][i], A[i][j], A[j][j], n, i, j);
                    Matrix<double> H = Q.Transpose()*A*Q;
                    A = H;
                    k++;
                }
            }
        barrier /= 10;
    }
    cout << "iteration's amount. Jacobi: " << 2*l*n*n + (20 + 2*n*n*n)*k << endl;
    return A;
}
