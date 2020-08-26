// header file for Runge Kutta Integrator
#ifndef __UTRAJECTORYTRACER_HPP_
#define __UTRAJECTORYTRACER_HPP_
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "MagneticField.hpp"
#include "constants.hpp"
#include "igrf.hpp"

/*
Class that traces / tracks the trajectory of a single particle within
Earth's magnetic field.

This is different from the similar class TrajectoryTracer in that the
integration algorithm is performed individually instead of in vector form
(i.e. evaluates each ODE as its own rather as evaluating the ODE as a vector).


Class Members
--------
- bfield_ (MagneticField instance):
    The MagneticField object that contains the information pertaining the
    Earth's magnetic field model (default is dipole model).
- charge_ (double) :
    The charge of the particle in units of Coulombs (default 1.602e-19 C, i.e.
1e).
- mass_ (double) :
    The (rest) mass of the particle in units of kilograms
(default 1.672e-27 kg, i.e. proton mass).
- escape_radius_ (double) :
    The radial distance relative to Earth's center in which the particle is
    set to have "escaped" Earth, i.e. a valid cosmic ray coming from some
    astrophysical source. Units are in kilometers (default 10*Earth's radius)
- stepsize_ (double) :
    The stepsize for the integration process (default 1e-5)
- max_iter_ (int) :
    The maximum number of iterations performed for the integration process
    (default 10000)
- particle_escaped_ (bool) :
    True if particle has effectively "escaped" from Earth (i.e. when the
    radial component of the trajectory > escape_radius_) (default false).
*/
class uTrajectoryTracer {
 private:
  MagneticField bfield_;  // the magnetic field
  double charge_;         // charge of the particle in coulumbs
  double mass_;           // mass of the particle in kg
  double escape_radius_;  // radius in which we set for particle to escape
  double stepsize_;       // step size
  int max_iter_;          // maximum step size
  // binary value to store if particle has escaped or not
  bool particle_escaped_;

  double final_time_;  // the final time of the trajectory
  std::array<double, 6> final_sixvector_; // the final six-vector of the trajectory

  // container for trajectory data throughout integration
  struct {
    double t;
    double r;
    double theta;
    double phi;
    double pr;
    double ptheta;
    double pphi;
  } traj_vector_;

  /* velocity in the radial component */
  inline double dr_dt(double pr);
  /* velocity in the polar component */
  inline double dtheta_dt(double r, double ptheta);
  /* velocity in the azimuthal component */
  inline double dphi_dt(double r, double theta, double pphi);
  /* acceleration in r component (time derivative of pr) from Lorentz force */
  inline double dpr_dt(double r, double theta, double phi, double pr,
                       double ptheta, double pphi);
  /* acceleration in theta component (time derivative of ptheta) from Lorentz
   * force */
  inline double dptheta_dt(double r, double theta, double phi, double pr,
                           double ptheta, double pphi);
  /* acceleration in phi component (time derivative of pphi) from Lorentz
  force */
  inline double dpphi_dt(double r, double theta, double phi, double pr,
                         double ptheta, double pphi);
  /* the Lorentz factor, computed from spherical components */
  inline double gamma(double pr, double ptheta, double pphi);
  /* evaluate / perform one step of the Runge-Kutta algorithm
   takes the 7-vector consisting of the following format:
   [t, r, theta, phi, pr, ptheta, pphi]
   and returns the modified 7-vector in that same format */
  void perform_rkstep();

 public:
  /* Default Constructor for TrajectoryTracer class

  Creates an instance of the TrajectoryTracer, that is, the object that keeps
  track of a single particle trajectory in Earth's magnetic field.

  The default constructor initializes the object with the default values
  provided for the members.

  Parameters
  ------------
  None
*/
  uTrajectoryTracer();
  /* Constructor for TrajectoryTracer class

  Creates an instance of the TrajectoryTracer, that is, the object that keeps
  track of a single particle trajectory in Earth's magnetic field.

  This constructor requires 2 parameters, and the rest may be optional.

  Required Parameters
  -------------------
  - charge (int) :
        The charge of the particle in units of electrons.
  - mass (double) :
        The mass of the particle in units of GeV.

  Optional Parameters
  -------------------
  - escape_radius (double) :
        The radius in which the particle has "escaped" relative to
        Earth's center in units of km (default 10*RE)
  - stepsize (double) :
        The step size of the integration (default 1e-5)
  - max_iter (int) :
        The maximum number of iterations performed in the integration process
        (default 10000)
  - bfield_type (char) :
        The type of Magnetic Field to evaluate the trajectory with. Only types
        'd' (for dipole model) or 'i' (IGRF model) are allowed (default 'd').
  - igrf_params (std::pair<std::string, double>) :
        Parameters required for instantiating the IGRF model. The first entry
        contains the path of the directory in which the .COF data file is
        stored (default to the author's path to the directory). The second entry
        contains the date in which the evaluation of the trajectory is requested
        in decimal date (default 2020.).
  */
  uTrajectoryTracer(const int charge, const double &mass,
                    const double &escape_radius = 10. * constants::RE,
                    const double &stepsize = 1e-5, const int max_iter = 10000,
                    const char bfield_type = 'd',
                    const std::pair<std::string, double> &igrf_params = {
                        "/home/keito/devel/gtracr/data", 2020.});

  // return the charge of the particle associated with the Runge-Kutta
  // integrator
  const double &charge() { return charge_ / constants::ELEMENTARY_CHARGE; }
  // return the mass of the particle associated with the Runge-Kutta integrator
  const double &mass() { return mass_ / constants::KG_PER_GEVC2; }
  // return the escape radius of the tracer
  const double &escape_radius() { return escape_radius_; }
  // return the step size of the Runge-Kutta integrator
  const double &stepsize() { return stepsize_; }
  // return the maximum number of steps of the integrator
  int max_iter() { return max_iter_; }
  // return the boolean if particle has escaped or not
  bool particle_escaped() { return particle_escaped_; }
  // the final time of the trajectory
  const double &final_time() { return final_time_; }
  // the final sixvector of the trajectory
  const std::array<double, 6> &final_sixvector() { return final_sixvector_; }

  /* Evaluates the trajectory of the particle using a 4th-order Runge Kutta
  algorithm.

  Parameters
  -----------
  - t0 (double) :
      the initial time in seconds
  - vec0 (std::array<double, 6>) :
       the initial six-vector (r0, theta0, phi0, pr0, ptheta0, pphi0) at time t0

  Returns
  --------
  None

  */
  void evaluate(double t0, std::array<double, 6> vec0);

  /* Evaluates the trajectory of the particle using a 4th-order Runge Kutta
  algorithm and return a map that contains the information of the particle
  trajectory. This will most often be used for debugging purposes to see the
  actual trajectory.

  Parameters
  -----------
  - t0 (double) :
      the initial time in seconds
  - vec0 (std::array<double, 6>) :
       the initial six-vector (r0, theta0, phi0, pr0, ptheta0, pphi0) at time t0

  Returns
  --------
  - trajectory_data (std::map<std::string, std::vector<double> >) :
      the trajectory information, that is, the time array of the trajectory
      and the six-vector of the trajectory in std::vectors.
      Notes:
      - keys are ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]
      - the final point of the trajectory is also contained in the map
      - std::vectors (dynamic arrays) are used since the length of each
        trajectory is not known at compile time.

  */
  std::map<std::string, std::vector<double>> evaluate_and_get_trajectory(
      double t0, std::array<double, 6> vec0);
};      // end class uTrajectoryTracer
#endif  //__UTRAJECTORYTRACER_HPP_