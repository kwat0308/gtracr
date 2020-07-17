//header file for Runge Kutta Integrator
#ifndef __TRAJECTORYTRACER_H_
#define __TRAJECTORYTRACER_H_
#include "MagneticField.h"

class TrajectoryTracer
{
private:
    MagneticField bfield_;  // the magnetic field
    int charge_; // charge of the particle in e
    double mass_; // mass of the particle in GeV/c^2
    double escape_radius_;  // radius in which we set for particle to escape
    double stepsize_;     // step size
    int max_iter_;   // maximum step size

    // arrays to store trajectory coordinates and momentum
    // this should change to pointers for a faster performance in the future
    std::vector<double> time_arr;
    std::vector<double> x_arr;
    std::vector<double> y_arr;
    std::vector<double> z_arr;
    std::vector<double> px_arr;
    std::vector<double> py_arr;
    std::vector<double> pz_arr;

    // // binary value to store if particle has escaped or not
    bool particle_escaped_;

public:
    // default constructor
    TrajectoryTracer();                                         
    // construct using charge and mass of particle
    // optional parameters:
    // - escape_radius: radius from center of Earth in which particle has effectively escaped
    // - step_size: the step size of the integration
    // - max_iter: the maximum number of iterations for integration
    TrajectoryTracer(const int charge, const double &mass, const double& escape_radius, 
                    const double &stepsize, const int max_iter); 
    // TrajectoryTracer(const int charge, const double &mass);                 // charge and mass
    // TrajectoryTracer(const int charge, const double &mass, const double &stepsize); // charge, mass, and stepsize
    // Destructor
    ~TrajectoryTracer() { };
    // copy constructor and assignment operator
    TrajectoryTracer(const TrajectoryTracer &traj_tracer);
    TrajectoryTracer &operator=(const TrajectoryTracer &traj_tracer);
    
    //return the charge of the particle associated with the Runge-Kutta integrator
    const int charge() {return charge_; }
    //return the mass of the particle associated with the Runge-Kutta integrator
    const double &mass() { return mass_; }
    // return the escape radius of the tracer
    const double& escape_radius() { return escape_radius_;}
    //return the step size of the Runge-Kutta integrator
    const double &stepsize() { return stepsize_; } 
    // return the maximum number of steps of the integrator
    const int max_iter() {return max_iter_;}
    // return the boolean if particle has escaped or not
    const bool particle_escaped() {return particle_escaped_;}
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
    std::array<double, 7>& rk_step(std::array<double, 7>& vec);

    // perform the runge kutta integration 
    // inputs required are the initial values [t0, r0, theta0, phi0, pr0, ptheta0, pphi0]
    // and the escape radius in which a particle has effectively escaped
    // returns a map that contains arrays of time, coordinate, 
    // and momentum data in Cartesian coordinates, with keys indicated as follows:
    // ["t", "x", "y", "z", "px", "py", "pz"]
    // also return within this map a 6-vector of the final point with the key "final_values"
    // in spherical coordinates (not Cartesian)
    // this would be passed into Python and converted into a TrajectoryPoint object
    // so that we can get the latitude and longitudes of this final location
    std::map<std::string, std::vector<double> > evaluate(std::array<double, 7>& initial_values);
    
    
};

#endif //__TRAJECTORYTRACER_H_