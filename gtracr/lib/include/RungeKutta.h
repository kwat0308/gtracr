//header file for Runge Kutta Integrator
#ifndef __RUNGEKUTTA_H_
#define __RUNGEKUTTA_H_

class MagneticField;

class RungeKutta
{
private:
    MagneticField *bfield;
    double coeff; // mass / charge
    double h;     // step size

public:
    // Constructors
    RungeKutta(const int, const double &);                 // charge and mass
    RungeKutta(const int, const double &, const double &); // charge, mass, and stepsize
    // Destructor
    ~RungeKutta() { delete bfield; };

    //getters
    double stepSize() { return h; }

    // setters
    void set_stepSize(double h_new) { h = h_new; }

    // ODEs
    double dvrdt(double, double, double, double, double, double);

    double dvthetadt(double, double, double, double, double, double);

    double dvphidt(double, double, double, double, double, double);

    // utility functions
    double velocity(double, double, double, double, double);
    double gamma(double, double, double, double, double);
    // evaluator
    // pybind11::array_t<double> evaluate( pybind11::array_t<double>);
    std::vector<double> evaluate(const std::vector<double> &);
};

#endif //__RUNGEKUTTA_H_