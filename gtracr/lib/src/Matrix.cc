/* Matrix class created specifically for this package
*/

#include <vector>
#include <utility>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include "Vector.h"
#include "Matrix.h" // header file

// default constructor, default is 2x2
Matrix::Matrix() 
    :rsz{1}, csz{1}, m{new double[2]}
{
    m[0] = 0.0;
    m[1] = 0.0;
    m[2] = 0.0;
    m[3] = 0.0;
}

// constructor with specified rowsize and columnsize rs, cs
// initialize a rs x cs matrix (filled with zero value)
Matrix::Matrix(const int rs, const int cs)
    : rsz{rs}, csz{cs}, m{new double[rs * cs]}
{
    for (int i = 0; i < rsz; ++i)
    {
        for (int j = 0; j < csz; ++j)
        {
            m[i * csz + j] = 0.0; // assign zero value to each index in memory
        }
    }
}

// constructor with specified rowsize and columnsize rs, cs
// with given data as a 1-d array (2-d array cant be used as parameter unless we know row length)
// initialize a rs x cs matrix
Matrix::Matrix(const int rs, const int cs, const double *dp)
    : rsz{rs}, csz{cs}, m{new double[rs * cs]}
{
    if (sizeof(dp) == sizeof(m)) // if data size == pointer size
    {
        for (int i=0; i<rs*cs; ++i) {
                m[i] = dp[i];
        }
    }
    else
    {
        throw std::runtime_error("Pointer dimensions dont match!");
    }
}

// create from std::vector
Matrix::Matrix(const int rs, const int cs, std::vector<double> vec) 
    : rsz{rs}, csz{cs}, m{new double[rs*cs]}
{
    // assume that vector entries are given as:
    // {0,1,2,3,4,5} == {{0, 1, 2}, {3, 4, 5}} == {{0, 1}, {2, 3}, {4, 5}}
    // for now this is ok, needs to be corrected to account for rowsize and columnsize later on
    for (int k=0; k<rs*cs; ++k) {
        m[k] = vec[k];
        // std::cout<< m[k] << std::endl;
    }
}

// // copy constructor
Matrix::Matrix(const Matrix &mat)
    : rsz{mat.rsz}, csz{mat.csz}, m{new double[mat.rsz * mat.csz]}
{
    if (dim_equal(mat))
    {
        std::copy(mat.m, mat.m + (mat.rsz * mat.csz), m);
    }
    else
    {
        throw std::runtime_error("Dimensions must be equal!");
    }
}

// copy assignment
Matrix &Matrix::operator=(const Matrix &mat)
{
    if (dim_equal(mat))
    {
        double *p = new double[mat.rsz * mat.csz];
        // m = mat.m;
        std::copy(mat.m, mat.m + (mat.rsz * mat.csz), p);
        delete[] m;
        m = p;
        rsz = mat.rsz;
        csz = mat.csz;
        return *this;
    }
    else
    {
        throw std::runtime_error("Dimensions must be equal!");
    }
}

// move constructor
Matrix::Matrix(Matrix &&mat)
    : rsz{mat.rsz}, csz{mat.csz}, m{new double[mat.rsz * mat.csz]}
{
    if (dim_equal(mat))
    {
        mat.rsz = 0;
        mat.csz = 0;
        mat.m = nullptr;
    }

    else
    {
        throw std::runtime_error("Dimensions must be equal!");
    }
}

// move assignment
Matrix &Matrix::operator=(Matrix &&mat)
{
    if (dim_equal(mat))
    {
        delete[] m;
        m = mat.m;
        rsz = mat.rsz;
        csz = mat.csz;

        mat.m = nullptr;
        mat.rsz = 0;
        mat.csz = 0;
        return *this;
    }

    else
    {
        throw std::runtime_error("Dimensions must be equal!");
    }
}

// operator overloading for subscripting
double Matrix::operator()(const int i, const int j)
{
    return m[i * csz + j];
}

// check if two matrices have the same dimensions
bool Matrix::dim_equal(const Matrix &M)
{
    return rsz == M.rows() &&
           csz == M.cols();
}

// inner product between two matrices
double Matrix::inner_prod(const Matrix &M)
{
    double inner_prod{0.};

    if (dim_equal(M))
    {
        for (int i = 0; i < rsz; ++i)
        {
            for (int j = 0; j < csz; ++j)
            {
                inner_prod += m[i * csz + j] * M.get_value(i, j);
            }
        }

        return inner_prod;
    }
    else
    {
        throw std::runtime_error("Dimensions of Matrices not equal!");
    }
}

// norm of matrix
double Matrix::norm()
{
    double norm{0.};

    for (int i = 0; i < rsz; ++i)
    {
        for (int j = 0; j < csz; ++j)
        {
            double re = get_value(i,j);
            norm += re * re;
        }
    }

    return sqrt(norm);
}

// dot product
// matrix dot vector
std::vector<double> Matrix::dot(std::vector<double> vec) 
{
    std::vector<double> newvec = std::vector<double>(csz, 0.);
    if (vec.size()==csz) {
        for (int j=0; j<csz; ++j) {
            double sum = 0.0;
            for (int i=0; i<rsz; ++i) {
                sum += m[i*csz + j]*vec[i];
            }
            newvec[j] = sum;
        }
    }
    return newvec;
}

// addition
Matrix Matrix::operator+(const Matrix& mat) 
{
    if (dim_equal(mat))
    {
        Matrix newmat = Matrix(rsz, csz);

        for (int i = 0; i < rsz; ++i)
        {
            for (int j = 0; j < csz; ++j)
            {
                newmat.set_value(i, j, m[i * csz + j] + mat.get_value(i, j));
            }
        }

        return newmat;
    }
    else
    {
        throw std::runtime_error("Dimensions of Matrices not equal!");
    }
}

// print matrix
void Matrix::print()
{
    int rszlim{50}; // row/column size limit if "too large"

    std::cout << "Matrix (" << rsz << "-by-" << csz << "): " << std::endl;
    std::cout << '{' << std::endl;
    if (rsz > rszlim)
    { // for large rows print using ... notation
        for (int i = 0; i < 3; ++i)
        {
            print_row(i);
        }
        std::cout << "\t..." << std::endl;
        for (int j = 3; j > 0; --j)
        {
            print_row(rsz - j);
        }
    }
    else
    { // otherwise print the whole matrix
        for (int i = 0; i < rsz; ++i)
        {
            print_row(i);
        }
    }
    std::cout << '}' << std::endl;
}

void Matrix::print_row(const int i)
{
    int cszlim{50}; // column size limit

    std::cout << "\t{";
    if (csz > cszlim)
    { // for large columns print using ... notation
        for (int j = 0; j < 3; ++j)
        {
            std::cout << get_value(i, j) << ' ';
        }
        std::cout << "... ";
        for (int j = 3; j > 0; --j)
        {
            std::cout << get_value(i, csz - j) << ' ';
        }
    }
    else
    { // otherwise print the whole matrix
        for (int j = 0; j < csz; ++j)
        {
            std::cout << get_value(i, j) << ' ';
        }
    }

    std::cout << '}' << std::endl;
}
