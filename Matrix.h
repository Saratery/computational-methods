#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED
#include <vector>
#include <fstream>
#include <string>

using namespace std;

class Matrix_exception
{
    string message;
    public:
    int code;
    Matrix_exception(const string& s, const int& n) : message(s), code(n) {};
    virtual void Print_message(ostream &output) const {output << message << endl;};
};

template <class T>
class Matrix;

template <class T>
ostream& operator<<(ostream& out, Matrix<T>& m);

template <class T>
istream& operator>>(istream& in, Matrix<T>& m);

template <class T>
class Matrix
{
	private:

    vector<vector<T> > a; 
    unsigned short int height; 
    unsigned short int length; 

    public:
    Matrix(); 
    Matrix(const unsigned short int &n, const unsigned short int &m); 
    Matrix(const Matrix& orig); 

    virtual void Print(ostream &output) const;

    virtual unsigned short int getHeight() const;
    virtual unsigned short int getLength() const;

    void Clear();

    virtual Matrix& operator=(const Matrix& orig);

    vector<T> operator[](const unsigned short int& index) const;
    vector<T>& operator[](const unsigned short int& index);

    virtual Matrix operator+(const Matrix& second); 
    virtual Matrix operator-(const Matrix& second); 

    virtual Matrix operator*(const Matrix& second);
    virtual Matrix operator*(const T& t);

    virtual Matrix Transpose() const;

    virtual Matrix operator-();

    virtual void swapLines(const unsigned short int &i, const unsigned short int &j);
    virtual void swapColumns(const unsigned short int &i, const unsigned short int &j);

    virtual void equalChange(const unsigned short int &i, const unsigned short int &j, const T& t);

    virtual ~Matrix(); 
};

template <class T>
Matrix<T>::Matrix()
:height(0), length(0)
{
    a.clear();
}

template <class T>
Matrix<T>::Matrix(const unsigned short int &n, const unsigned short int &m)
: a(n + 1, vector<T>(m + 1)), height(n), length(m)
{}

template <class T>
Matrix<T>::Matrix(const Matrix& orig)
{
    height = orig.height;
    length = orig.length;
    a = orig.a;
}

template <class T>
void Matrix<T>::Print(ostream &output) const
{
    for (int i = 1; i <= height; i++)
    {
        for (int j = 1; j <= length; j++)
        {
            output << a[i][j] << ' ';
        }
        output << '\n';
    }
}

template <class T>
unsigned short int Matrix<T>::Get_height() const
{
    return height;
}

template <class T>
unsigned short int Matrix<T>::Get_length() const
{
    return length;
}

template <class T>
void Matrix<T>::Clear()
{
    for (int i = 0; i < a.size(); i++)
    {
        a[i].clear();
    }
}

template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix& orig)
{
    if (this == &orig)
    {
        return *this;
    }
    height = orig.height;
    length = orig.length;
    a = orig.a;
    return *this;
}

template <class T>
vector<T> Matrix<T>::operator[](const unsigned short int& index) const
{
    if ((index <= height) && (index > 0))
        {
            return a[index];
        }
        else
        {
            throw Matrix_exception("Incorrect index", 4);
        }
}

template <class T>
vector<T>& Matrix<T>::operator[](const unsigned short int& index)
{
    if ((index <= height) && (index > 0))
    {
        return a[index];
    }
    else
    {
        throw Matrix_exception("Bad index", 4);
    }
}

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix& second)
{
    if ((height != second.height) || (length != second.length))
    {
        throw Matrix_exception("Addition isn't defined for this type of matrices", 1);
    }
    else
    {
        Matrix m(height, length);
        for (int i = 1; i <= height; i++)
        {
            for (int j = 1; j <= length; j++)
            {
                m.a[i][j] = this->a[i][j] + second.a[i][j];
            }
        }
        return m;
    }
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix& second)
{
    if ((height != second.height) || (length != second.length))
    {
        throw Matrix_exception("Subtraction isn't defiend for such matrices", 1);
    }
    else
    {
        Matrix m(height, length);
        for (int i = 1; i <= height; i++)
        {
            for (int j = 1; j <= length; j++)
            {
                m.a[i][j] = this->a[i][j] - second.a[i][j];
            }
        }
        return m;
    }
}

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix& second)
{
    if (length  != second.height)
    {
        throw Matrix_exception("Multiplication isn't defiend for such matrices", 2);
    }
    else
    {
        Matrix m(height, second.length);
        for (int i = 1; i <= height; i++)
        {
            for (int j = 1; j <= second.length; j++)
            {
                T temp = a[i][1]*second.a[1][j];
                for (int k = 2; k <= length; k++)
                {
                    temp += a[i][k]*second.a[k][j];
                }
                m.a[i][j] = temp;
            }
        }
        return m;
    }
}

template <class T>
Matrix<T> Matrix<T>::operator*(const T& t)
{
    Matrix m = *this;
    for (int i = 1; i <= height; i++)
    {
        for (int j = 1; j <= length; j++)
        {
            m.a[i][j] *= t;
        }
    }
    return m;
}

template <class T>
Matrix<T> Matrix<T>::Transpose() const
{
    Matrix m(length, height);
    for (int i = 1; i <= height; i++)
    {
        for (int j = 1; j <= length; j++)
        {
            m.a[j][i] = a[i][j];
        }
    }
    return m;
}

template <class T>
Matrix<T> Matrix<T>::operator-()
{
    Matrix m1(height, length);
    for (int i = 1; i <= height; i++)
    {
        for (int j = 1; j <= length; j++)
        {
            m1.a[i][j] = -a[i][j];
        }
    }
    return m1;
}

template <class T>
void Matrix<T>::swapLines(const unsigned short int &i, const unsigned short int &j)
{
    if ((i > height) || (j > height))
    {
        throw Matrix_exception("Wrong number of line", 3);
    }
    else
    if (i == j)
    {
        return;
    }
    a[i].swap(a[j]);
}

template <class T>
void Matrix<T>::swapColumns(const unsigned short int &i, const unsigned short int &j)
{
    if ((i > length) || (j > length))
    {
        throw Matrix_exception("Wrong number of column", 3);
    }
    else
    if (i == j)
    {
        return;
    }
    T temp;
    for (int k = 1; k <= height; k++)
    {
        temp = a[k][i];
        a[k][i] = a[k][j];
        a[k][j] = temp;
    }
}

template <class T>
void Matrix<T>::equalChange(const unsigned short int &i, const unsigned short int &j, const T& t)
{
    for (int k = 1; k <= length; k++)
    {
        a[i][k] += a[j][k]*t;
    }
}

template <class T>
Matrix<T>::~Matrix()
{
}

template <class T>
ostream& operator<<(ostream& out, Matrix<T>& m)
{
    m.Print(out);
    return out;
}

template <class T>
istream& operator>>(istream& in, Matrix<T>& m)
{
    for (int i = 1; i <= m.Get_height(); i++)
    {
        for (int j = 1; j <= m.Get_length(); j++)
        {
            in >> m[i][j];
        }
    }
    return in;
}


#endif // MATRIX_H_INCLUDED
