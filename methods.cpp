#include <iostream>
#include <stdio.h>
#include "Matrix.h"
#include <vector>
#include <fstream>
#include <string>
#include "SLAU solutions.h"
#include <math.h>

const double eps = 1e-20;

Matrix<double> Gauss_full_choice(Matrix<double> A, Matrix<double> b)

{
    if ((A.getHeight() != b.getHeight()) || (b.getLenght() != 1))
    {
        throw Matrix_exception("Impossible system", 5);
    }
    else
    if (A.getHeight() != A.getLenght())
    {
        throw Matrix_exception("A is not a square matrix", 6);
    }
    else
    {
        const short unsigned int n = A.getHeight();
        Matrix<double> x(n, 1);
        int col[n + 1];

        for (int i = 1; i <= n; i++)
        {
            col[i] = i;
        }

        for (int k = 1; k < n; k++)
        {

            int maxi = k;
            int maxj = k;
            double max = fabs(A[k][k]);
            for (int i = k; i <= n; i++)
            {
                for (int j = k; j <= n; j++)
                {
                    if (fabs(A[i][j]) > max)
                    {
                        max = fabs(A[i][j]);
                        maxi = i;
                        maxj = j;
                    }
                }
            }


            A.swapColumns(k, maxj);
            A.swapLines(k, maxi);
            b.swapLines(k, maxi);
            int temp = col[k];
            col[k] = col[maxj];
            col[maxj] = temp;

            if (fabs(A[k][k]) < eps)
            {
                throw Matrix_exception("Division by zero in reverse substitution. Seems like det(A) = 0", 7);
            }

            for (int i = k + 1; i <= n; i++)
            {
                b.Equal_change(i, k, -(A[i][k]/A[k][k]));
                A.Equal_change(i, k, -(A[i][k]/A[k][k]));
                A[i][k] = 0;
            }
           
        }
        
        for (int i = n; i >= 1; i--)
        {
            if (fabs(A[i][i]) < eps)
            {
                throw Matrix_exception("Division by zero in reverse substitution. Seems like det(A) = 0", 7);
            }
            double temp = 0;
            for (int j = n; j > i; j--)
            {
                temp += A[i][j]*x[j][1];
            }
            x[i][1] = (b[i][1] - temp)/A[i][i];
        }


        Matrix<double> y = x;
        for (int i = 1; i <= n; i++)
        {
            y[col[i]] = x[i];
        }
        x = y;
        return x;
    }
}

double norm_2(const Matrix<double>& x)

{
    if (x.getLenght() != 1)
    {
        throw Matrix_exception("Can't count euclidean norm of non-vector", 8);
    }
    else
    {
        return sqrt((x.Transpose()*x)[1][1]);
    }
}

Matrix<double> Householder(Matrix<double> A, Matrix<double> b)

{
    if ((A.getHeight() != b.getHeight()) || (b.getLenght() != 1))
    {
        throw Matrix_exception("Impossible system", 5);
    }
    else
    if (A.getHeight() != A.getLenght())
    {
        throw Matrix_exception("A is not a square matrix", 6);
    }
    else
    {
        const short unsigned int n = A.getHeight();
        Matrix<double> x(n, 1);

        Matrix<double> I(n, n);

        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                if (i != j)
                {
                    I[i][j] = 0;
                }
                else
                {
                    I[i][j] = 1;
                }
            }
        }

        Matrix<double> Q = I;

        for (int i = 1; i < n; i++)
        {
            int k = n - i + 1;

            Matrix<double> e(k, 1);
            for (int j = 2; j <= k; j++)
            {
                e[j][1] = 0;
            }
            e[1][1] = 1;

            Matrix<double> y(k, 1);
            for (int j = 1; j <= k; j++)
            {
                y[j][1] = A[i + j - 1][i];
            }

            double l = norm_2(y);

            Matrix<double> u(k, 1);

            if (norm_2(y + e*norm_2(y)) > norm_2(y - e*norm_2(y)))
            {
                u = y + e*l;
            }
            else
            {
                u = y - e*l;
            }

            u = u*(double(1)/norm_2(u));

            Matrix<double> Ik(k,k);
            for (int i1 = 1; i1 <= k; i1++)
            {
                for (int j1 = 1; j1 <= k; j1++)
                {
                    if (i1!= j1)
                    {
                        Ik[i1][j1] = 0;
                    }
                    else
                    {
                        Ik[i1][j1] = 1;
                    }
                }
            }

            Matrix<double> Hk = Ik - u*u.Transpose()*2;

            Matrix<double> H = I;

            for (int i1 = 1; i1 <= k; i1++)
            {
                for (int j1 = 1; j1 <= k; j1++)
                {
                    H[i + i1 - 1][i + j1 - 1] = Hk[i1][j1];
                }
            }


            A = H*A;
            b = H*b;

            for (int j = i + 1; j <= n; j++)
            {
                A[j][i] = 0;
            }
        }


        for (int i = n; i >= 1; i--)
        {
            if (fabs(A[i][i]) < eps)
            {
                throw Matrix_exception("Division by zero in reverse substitution. Seems like det(A) = 0", 7);
            }
            double temp = 0;
            for (int j = n; j > i; j--)
            {
                temp += A[i][j]*x[j][1];
            }
            x[i][1] = (b[i][1] - temp)/A[i][i];
        }

        return x;
    }
}
