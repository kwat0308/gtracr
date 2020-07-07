// Vector class derived from Matrix class
// #include "Matrix.h"
#include <vector>
#include <iostream>
#include "Vector.h"

// Constructor
// default
Vector::Vector()
    : sz{1}
{
    Matrix(0, 1.);
}

// with given size
Vector::Vector(const int size)
    : sz{size}
{
    Matrix(0, sz);
}

// using pointers
Vector::Vector(const int size, const double *dp)
    : sz{size}
{
    Matrix(0, sz, dp);
}

// using std vectors
Vector::Vector(const int size, std::vector<double> vec)
    : sz{size}
{
    Matrix(0, sz);
    for (int i = 0; i < sz; ++i)
    {
        m[i] = vec[i];
    }
}

// dot product with matrix and vector
double Vector::dot(const Vector &v)
{
    inner_prod(v);
}

// sum of vector
double Vector::sum()
{
    double sum = 0.;

    for (int i = 0; i < sz; ++i)
    {
        sum += m[i];
    }
    return sum;
}

// addition between two vectors
Vector Vector::operator+(const Vector &v)
{
    if (dim_equal(v))
    {
        Vector newvec = Vector(sz);
        for (int i = 0; i < sz; ++i)
        {
            newvec.set_value(i, m[i] + v[i]);
        }
        return newvec;
    }
    else
    {
        throw std::runtime_error("Vector sizes must be equal!");
    }
}

// print vector
void Vector::print()
{
    std::cout << "Vector of size " << sz << ": {";
    for (int i = 0; i < sz; ++i)
    {
        std::cout << m[i] << ',';
    }
    std::cout << "}" << std::endl;
}
