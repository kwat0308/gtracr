/*
SixVector class for aids in array addition / multiplication
*/

#ifndef __SIXVECTOR_HPP_
#define __SIXVECTOR_HPP_

struct SixVector {

    /*
    Vector class containing the coordinate and momentum components of the trajectory in spherical coordinates.

    Class Members
    -------------
    - r_, theta_, phi_ (double)  : coordinate at an instance in time in trajectory
    - pr_, ptheta_, pphi_ (double) : momentum at any instance in time in trajectory
    */ 
  
    double r, theta, phi, pr, ptheta, pphi;

    /*
    Default constructor. Sets all parameters to default value of 0.
    */    
    SixVector() : r{0.}, theta{0.}, phi{0.}, pr{0.}, ptheta{0.}, pphi{0.} {};

    /*
    Constructor that sets coordinate and momentum values based on given input (double)
    */
    SixVector(double _r, double _theta, double _phi, double _pr, double _ptheta, double _pphi) :
        r{_r}, theta{_theta}, phi{_phi}, pr{_pr}, ptheta{_ptheta}, pphi{_pphi} {};

    /*
    Constructor that initializes with array of elements in std::array
    */
    SixVector(std::array<double, 6> arr) :
        r{arr[0]}, theta{arr[1]}, phi{arr[2]}, pr{arr[3]}, ptheta{arr[4]}, pphi{arr[5]} {};
        

    /* 
    Operator overloading of element-wise addition.
    */
    SixVector &operator+(const SixVector& vec) {
        r += vec.r;
        theta += vec.theta;
        phi += vec.phi;
        pr += vec.pr;
        ptheta += vec.ptheta;
        pphi += vec.pphi;

        return *this;
    }

    /* 
    Operator overloading of element-wise multiplication.
    */
    SixVector &operator*(const SixVector& vec) {
        r *= vec.r;
        theta *= vec.theta;
        phi *= vec.phi;
        pr *= vec.pr;
        ptheta *= vec.ptheta;
        pphi *= vec.pphi;

        return *this;
    }

    /*
    Operator overloading of RIGHT SIDE constant multiplication.
    Note: this only works if scalar is on RIGHT SIDE, those on the LEFT will NOT WORK.
    */
    SixVector &operator*(const double val) {
        r *= val;
        theta *= val;
        phi *= val;
        pr *= val;
        ptheta *= val;
        pphi *= val;

        return *this;
    }

    /*
    Convert struct of SixVector to std::array. Required for pybind11 to interpret values.
    */
    std::array<double, 6> to_array(){return std::array<double, 6>{r, theta, phi, pr, ptheta, pphi};}

    /*
    print results
    */
    void print() {
        std::cout << "r: " << r << "\t theta: " << theta << "\t phi: " << phi << std::endl
                  << "pr: " << pr << "\t ptheta: " << ptheta << "\t pphi: " << pphi << std::endl;
    }

};

#endif /* __SIXVECTOR_HPP_*/