// header file for Vector class
#ifndef __VECTOR_H_
#define __VECTOR_H_

#include "Matrix.h"

class Vector : public Matrix {
    private:
        int sz;
    public:
        // constructors
        Vector();
        Vector(const int);
        Vector(const int, const double*);
        Vector(const int, std::vector<double>);
        // destructor
        ~Vector() {delete[] m;}
        // getter
        const int size() {return sz;}
        double operator[] (const int i) const {return m[i];} //subscripting
        // setter
        void set_value(const int i, const double& _val) {m[i] = _val;}
        // operations
        double dot(const Vector&);  // dot product
        double sum(); // sum of all elements
        Vector operator+(const Vector&);
        // print vector
        void print();

};

#endif //__VECTOR_H_