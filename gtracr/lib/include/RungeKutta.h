//header file for Runge Kutta Integrator
#ifndef __RUNGEKUTTA_H_
#define __RUNGEKUTTA_H_

class MagneticField;

class RungeKutta
{
private:
    MagneticField *bfield; 
    const double& coeff;  // mass / charge
    const double& h;  // step size

public:
    // Constructors
    RungeKutta(const int, const double &);  // charge and mass
    RungeKutta(const int, const double &, const double&);  // charge, mass, and stepsize
    // Destructor
    ~RungeKutta() {delete bfield;};
    // ODEs
    const double &dvrdt(const double&, const double&, const double&, const double&, const double&, const double&);

    const double &dvthetadt(const double&, const double&, const double&, const double&, const double&, const double&);

    const double &dvphidt(const double&, const double&, const double&, const double&, const double&, const double&);

    // utility functions
    const double &velocity(const double&, const double&, const double&, const double&, const double&);
    const double &gamma(const double&, const double&, const double&, const double&, const double&);
    // evaluator
    // pybind11::array_t<double> evaluate( pybind11::array_t<double>);
    std::vector<double> evaluate( std::vector<double>);
};

#endif //__RUNGEKUTTA_H_