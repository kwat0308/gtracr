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
    // set the step size of the integrator
    void set_stepsize(const double &_h) { h = _h; }
    // acceleration in r component (time derivative of vr) from Lorentz force
    double dvr_dt(double r, double theta, double phi, double vr, double vtheta, double vphi);
    // acceleration in theta component (time derivative of vtheta) from Lorentz force
    double dvtheta_dt(double r, double theta, double phi, double vr, double vtheta, double vphi);
    // acceleration in phi component (time derivative of vphi) from Lorentz force
    double dvphi_dt(double r, double theta, double phi, double vr, double vtheta, double vphi);
    // the magnitude of velocity, computed from spherical components
    double velocity(double vr, double vtheta, double vphi, double r, double theta);
    // the Lorentz factor, computed from spherical components
    double gamma(double vr, double vtheta, double vphi, double r, double theta);
    // evaluate / perform one step of the Runge-Kutta algorithm
    // takes the 7-vector consisting of the following format:
    // [t, r, theta, phi, vr, vtheta, vphi]
    // and returns the modified 7-vector in that same format
    std::array<double, 7>& evaluate(std::array<double, 7>& vec);
    // pybind11::array_t<double> evaluate( const pybind11::array_t<double>&);
    
};

#endif //__RUNGEKUTTA_H_