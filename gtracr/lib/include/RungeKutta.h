//header file for Runge Kutta Integrator
#ifndef __RUNGEKUTTA_H_
#define __RUNGEKUTTA_H_
#include "MagneticField.h"

class RungeKutta
{
private:
    MagneticField bfield_;  // the magnetic field
    int charge_; // charge of the particle in e
    double mass_; // mass of the particle in GeV/c^2
    double h;     // step size

public:
    // Constructors
    RungeKutta();                                          // default
    RungeKutta(const int charge, const double &mass);                 // charge and mass
    RungeKutta(const int charge, const double &mass, const double &stepsize); // charge, mass, and stepsize
    // Destructor
    ~RungeKutta() { };
    // copy constructor and assignment operator
    RungeKutta(const RungeKutta &rk);
    RungeKutta &operator=(const RungeKutta &rk);
    //return the step size of the Runge-Kutta integrator
    const double &stepsize() { return h; } 
    //return the charge of the particle associated with the Runge-Kutta integrator
    const int charge() {return charge_; }
    //return the mass of the particle associated with the Runge-Kutta integrator
    const double &mass() { return mass_; }
    // set the step size of the integrator
    void set_stepsize(const double &_h) { h = _h; }
    // set the charge of the particle associated with the integrator
    void set_charge(const int _charge) {charge_ = _charge;}
    // set the mass of the particle associated with the integrator
    void set_mass(const double& _mass) {mass_ = _mass;}
    // velocity in the radial component
    double dr_dt(double pr);
    // velocity in the polar component
    double dtheta_dt(double r, double ptheta);
    // velocity in the azimuthal component
    double dphi_dt(double r, double theta, double pphi);
    // acceleration in r component (time derivative of pr) from Lorentz force
    double dpr_dt(double r, double theta, double phi, double pr, double ptheta, double pphi);
    // acceleration in theta component (time derivative of ptheta) from Lorentz force
    double dptheta_dt(double r, double theta, double phi, double pr, double ptheta, double pphi);
    // acceleration in phi component (time derivative of pphi) from Lorentz force
    double dpphi_dt(double r, double theta, double phi, double pr, double ptheta, double pphi);
    // the Lorentz factor, computed from spherical components
    double gamma(double pr, double ptheta, double pphi);
    // evaluate / perform one step of the Runge-Kutta algorithm
    // takes the 7-vector consisting of the following format:
    // [t, r, theta, phi, pr, ptheta, pphi]
    // and returns the modified 7-vector in that same format
    std::array<double, 7>& evaluate(std::array<double, 7>& vec);
    // pybind11::array_t<double> evaluate( const pybind11::array_t<double>&);
    
};

#endif //__RUNGEKUTTA_H_