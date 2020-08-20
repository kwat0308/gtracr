// header file for Runge Kutta Integrator
#ifndef __TRAJECTORYTRACER_HPP_
#define __TRAJECTORYTRACER_HPP_
#include <array>
#include <map>
#include <string>
#include <vector>

#include "MagneticField.hpp"

class TrajectoryTracer {
 private:
  MagneticField bfield_;  // the magnetic field
  double charge_;         // charge of the particle in coulumbs
  double mass_;           // mass of the particle in kg
  double escape_radius_;  // radius in which we set for particle to escape
  double stepsize_;       // step size
  int max_iter_;          // maximum step size
  // binary value to store if particle has escaped or not
  bool particle_escaped_;

 public:
  // default constructor
  TrajectoryTracer();
  // construct using charge and mass of particle
  // optional parameters:
  // - escape_radius: radius from center of Earth in which particle has
  // effectively escaped
  // - step_size: the step size of the integration
  // - max_iter: the maximum number of iterations for integration
  TrajectoryTracer(const int charge, const double &mass,
                   const double &escape_radius, const double &stepsize,
                   const int max_iter, const char bfield_type);
  // TrajectoryTracer(const int charge, const double &mass);                 //
  // charge and mass TrajectoryTracer(const int charge, const double &mass,
  // const double &stepsize); // charge, mass, and stepsize Destructor
  ~TrajectoryTracer(){};
  // copy constructor and assignment operator
  // TrajectoryTracer(const TrajectoryTracer &traj_tracer);
  // TrajectoryTracer &operator=(const TrajectoryTracer &traj_tracer);

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

  /*
  Mathematical operations to speed up the evaluation process
  */
  /*

  */
  /* The ordinary differential equations that describes the motion
  of charge particles in Earth's magnetic field via the Lorentz force
  in spherical coordinates.

  Parameters
  -----------
  - t (double) :
      the time
  - vec (std::array<double, 6>) :
       the six-vector (r, theta, phi, pr, ptheta, pphi) at time t

  Returns
  --------
  - ode_lrz (std::array<double, 6>) :
       the ordinary differential equation for the six vector based on the
  Lorentz force equation
  */
  std::array<double, 6> ode_lrz(const double t,
                                const std::array<double, 6> &vec);

  /*
  The differential equation for the momentum in
  the radial component (dprdt)

  Parameters
  -----------
  - r, theta, phi (double):
        the coordinates at time t
  - br, btheta, bphi (double):
        the magnetic field components at time t
  */
  inline double dprdt(const double &r, const double &theta, const double &phi,
                      const double &br, const double &btheta,
                      const double &bphi);

  /*
  The lorentz factor
  */

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
  void evaluate(double &t0, std::array<double, 6> &vec0);

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
      double &t0, std::array<double, 6> &vec0);
};

#endif  //__TRAJECTORYTRACER_HPP_