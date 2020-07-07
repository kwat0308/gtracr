// contains class declarations of TrajectoryMatrix class
#ifndef __MATRIX_H_
#define __MATRIX_H_

class Vector;

class Matrix
{
    protected:
        //members
        int rsz;   //row size
        int csz;   //column size
        double *m; // pointer to 2-d array

    public:
        //constructor
        Matrix(); // default
        Matrix(const int, const int);  
        Matrix(const int, const int, const double *);  // data from pointer
        Matrix(const int, const int, std::vector<double>);  // from std vector
        // copy
        Matrix(const Matrix &);            // copy constructor
        Matrix &operator=(const Matrix &); //copy assignment

        // move
        Matrix(Matrix &&);            // move constructor
        Matrix &operator=(Matrix &&); // move assignment

        // operator overloading
        double operator()(const int, const int); // subscripting

        // setters
        void set_value(const int i, const int j, const double &newval) { m[i * csz + j] = newval; } // change single value

        // getters
        const int rows() const { return rsz; }                                            // number of rows
        const int cols() const { return csz; }                                            // number of columns
        // const std::pair<int, int> dim() const { return std::pair<int, int> dim (rsz, csz); }   // dimension
        const double get_value(const int i, const int j) const { return m[i * csz + j]; } // value at (i,j) coordinate

        // operations
        bool dim_equal(const Matrix &); // check if dimensions are equal
        double inner_prod(const Matrix &);
        double norm();                                               //Frobenius norm (standard matrix norm)
        std::vector<double> dot(std::vector<double>);  // dot product with vector
        Matrix operator+(const Matrix&); // addition
        // utility functions
        virtual void print();          // print matrix
        void print_row(const int); // auxiliary function for print_mat() to print each row
        // Matrix import_data(const int&, const int&, const std::string&);  // import data (knowing data dimensions)
        // destructor
        ~Matrix() { delete[] m; }
};

#endif //__MATRIX_H_