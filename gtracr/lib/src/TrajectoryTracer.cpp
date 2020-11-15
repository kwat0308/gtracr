/*
Trajectory Tracer class that traces the trajectory of the particle
by performing a 4th-order Runge Kutta numerical integration algorithm.
*/

#include "TrajectoryTracer.hpp"

// TrajectoryTracer class

/* Default Constructor for TrajectoryTracer class

  Creates an instance of the TrajectoryTracer, that is, the object that keeps
  track of a single particle trajectory in Earth's magnetic field.

  The default constructor initializes the object with the default values
  provided for the members.

  Members
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

  Parameters
  ------------
  None
*/
TrajectoryTracer::TrajectoryTracer()
    : bfield_{MagneticField()},
      charge_{constants::ELEMENTARY_CHARGE},
      mass_{0.938 * constants::KG_PER_GEVC2},
      start_altitude_{100. * (1e3)},
      escape_radius_{10. * constants::RE},
      stepsize_{1e-5},
      max_iter_{10000},
      particle_escaped_{false} {}

/* Constructor for TrajectoryTracer class

Creates an instance of the TrajectoryTracer, that is, the object that keeps
track of a single particle trajectory in Earth's magnetic field.

This constructor requires 2 parameters, and the rest may be optional.

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
TrajectoryTracer::TrajectoryTracer( double charge,  double mass,
                                     double start_altitude,
                                    double escape_radius /*= 10. * constants::RE*/,
                                    double stepsize /*= 1e-5*/,
                                    int max_iter /* = 10000*/,
                                   const char bfield_type /*= 'i'*/,
                                   const std::pair<std::string, double>
                                       &igrf_params /*=
      {"/home/keito/devel/gtracr/data",
      2020.}*/)
    : charge_{charge},
      mass_{mass},
      start_altitude_{start_altitude},
      escape_radius_{escape_radius},
      stepsize_{stepsize},
      max_iter_{max_iter},
      particle_escaped_{false} {
  switch (bfield_type) {
    case 'd':
      bfield_ = MagneticField();
      break;
    case 'i':
      // add file name to DATA_DIR (first component in igrf_params)
      std::string DATA_PATH = igrf_params.first + "/igrf13.json";
      double sdate = igrf_params.second;
      bfield_ = IGRF(DATA_PATH, sdate);
      break;
  }
}

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
SixVector TrajectoryTracer::ode_lrz(const double t,
                                                const SixVector& vec) {

  // get the lorentz factor
  double gmma = lorentz_factor(vec.pr, vec.ptheta, vec.pphi);
  double rel_mass = mass_ * gmma;

  // evaluate B-field
  std::array<double, 3> bf_values = bfield_.values(vec.r, vec.theta, vec.phi);
  double bf_r = bf_values[0];
  double bf_theta = bf_values[1];
  double bf_phi = bf_values[2];

  // get the momentum ODE
  // Note:
  // - charge is inverted to allow backtracking
  // - xx_lrz is the lorentz term in the ODE
  // - xx_sphcmp is the auxiliary terms due to acceleration
  //   in spherical coordiantes in the ODE

  // -q * (ptheta * Bphi - Btheta * pphi)
  double dprdt_lrz = -1. * charge_ * ((vec.ptheta * bf_phi) - (bf_theta * vec.pphi));
  // (ptheta^2 / r) + (pphi^2 / r)
  double dprdt_sphcmp = (((vec.ptheta * vec.ptheta) + (vec.pphi * vec.pphi)) / vec.r);
  // dprdt
  double dprdt = dprdt_lrz + dprdt_sphcmp;

  // q * (pr * Bphi - Br * pphi)
  double dpthetadt_lrz = charge_ * ((vec.pr * bf_phi) - (bf_r * vec.pphi));
  // (pphi^2 * cos(theta) / (r * sin(theta))) - (pr*ptheta / r)
  double dpthetadt_sphcmp =
      ((vec.pphi * vec.pphi * cos(vec.theta)) / (vec.r * sin(vec.theta))) - ((vec.pr * vec.ptheta) / vec.r);
  // dpthetadt
  double dpthetadt = dpthetadt_lrz + dpthetadt_sphcmp;

  // -q * (pr * Btheta - Br * ptheta)
  double dpphidt_lrz = -1. * charge_ * ((vec.pr * bf_theta) - (bf_r * vec.ptheta));
  // (pr * pphi / r) + ((vec.ptheta * pphi * cos(theta)) / (r * sin(theta)))
  double dpphidt_sphcmp =
      ((vec.pr * vec.pphi) / vec.r) + ((vec.ptheta * vec.pphi * cos(vec.theta)) / (vec.r * sin(vec.theta)));
  // dpphidt
  double dpphidt = dpphidt_lrz - dpphidt_sphcmp;

  // create the vector form of the ODE with the position components included
  // as well
  SixVector ode_lrz = SixVector(
      vec.pr, (vec.ptheta / vec.r), (vec.pphi / (vec.r * sin(vec.theta))), dprdt, dpthetadt, dpphidt);

  // SixVector result = ode_lrz * (1. / rel_mass);

  return ode_lrz * (1. / rel_mass);
}
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
void TrajectoryTracer::evaluate(const double &t0, std::array<double, 6> &vec0) {
  double h = stepsize_;  // step size in shorter notation

  // set the initial conditions
  double t = t0;
  SixVector vec = SixVector(vec0);

  vec.print();
  // TODO: make arrays into references for no copying

  // start the loop
  for (int i = 0; i < max_iter_; ++i) {

    k1_vec = ode_lrz(t, vec);
    k2_vec = ode_lrz(t + (0.5 * h), vec + (k1_vec * 0.5 * h));
    k3_vec = ode_lrz(t + (0.5 * h), vec + (k2_vec * 0.5 * h));
    k4_vec = ode_lrz(t + h, vec + k3_vec * h);
    k_vec = (k1_vec + (k2_vec * 2.) + (k3_vec * 2.) + k4_vec) * (h / 6.);
    // increment by weighted sum
    vec = vec + k_vec;
    t += h;  // increase time

    // breaking condition
    // if particle reaches escape radius
    if (vec.r > escape_radius_) {
      particle_escaped_ = true;
      break;
    }  // if (r > escape_radius_)

    // breaking condition
    // if particle reaches back onto Earth's surface again
    if (vec.r < start_altitude_ + constants::RE) {
      break;
    }  // if (r < start_altitude_ + constants::RE)
  }    // for (int i = 0; i < max_iter_; ++i)
  // store the final time and six-vector for checking purposes
  // the last recorded time and six-vector is the final six-vector / time
  final_time_ = t;
  final_sixvector_ = vec.to_array();
}  // evaluate

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
std::map<std::string, std::vector<double>>
TrajectoryTracer::evaluate_and_get_trajectory(double &t0,
                                              std::array<double, 6> &vec0) {
  double h = stepsize_;  // step size in shorter notation

  // set the initial conditions
  double t = t0;
  SixVector vec = SixVector(vec0);

  // initialize array to contain data
  std::vector<double> t_arr, r_arr, theta_arr, phi_arr, pr_arr, ptheta_arr,
      pphi_arr;

  // start the loop
  for (int i = 0; i < max_iter_; ++i) {
    // append the values first
    // how vec looks like:
    // (r, theta, phi, pr, ptheta, pphi) = vec

    t_arr.push_back(t);
    r_arr.push_back(vec.r);
    theta_arr.push_back(vec.theta);
    phi_arr.push_back(vec.phi);
    pr_arr.push_back(vec.pr);
    ptheta_arr.push_back(vec.ptheta);
    pphi_arr.push_back(vec.pphi);

    // evaluate the k-coefficients

    k1_vec = ode_lrz(t, vec);
    k2_vec = ode_lrz(t + (0.5 * h), vec + (k1_vec * 0.5 * h));
    k3_vec = ode_lrz(t + (0.5 * h), vec + (k2_vec * 0.5 * h));
    k4_vec = ode_lrz(t + h, vec + k3_vec * h);
    k_vec = (k1_vec + (k2_vec * 2.) + (k3_vec * 2.) + k4_vec) * (h / 6.);

    // increment by weighted sum
    vec = vec + k_vec;
    t += h;  // increase time

    // breaking condition
    // if particle reaches escape radius
    if (vec.r > escape_radius_) {
      particle_escaped_ = true;
      break;
    }  // if (r > escape_radius_)

    // breaking condition
    // if particle reaches back onto Earth's surface again
    if (vec.r < start_altitude_ + constants::RE) {
      break;
    }  // if (r < start_altitude_ + constants::RE)

  }  // for (int i = 0; i < max_iter_; ++i)

  // store the final time and six-vector for checking purposes
  // the last recorded time and six-vector is the final six-vector / time
  final_time_ = t;
  final_sixvector_ = vec.to_array();

  // create map that contains trajectory data
  std::map<std::string, std::vector<double>> trajectory_data = {
      {"t", t_arr},      {"r", r_arr},   {"theta", theta_arr},
      {"phi", phi_arr},  {"pr", pr_arr}, {"ptheta", ptheta_arr},
      {"pphi", pphi_arr}};

  return trajectory_data;
}  // evaluate_and_get_trajectory

/*
Returns the lorentz factor, evaluated from the momentum

Parameters
----------
- pr (const double &) :
      the momentum in the radial component
- ptheta (const double &) :
      the momentum in the polar direction
- pphi (const double &) :
      the momentum in the azimuthal direction

Returns
-------
- gamma (const double) :
      The lorentz factor for the particular momentum
*/
inline double TrajectoryTracer::lorentz_factor(const double &pr,
                                               const double &ptheta,
                                               const double &pphi) {
  // momentum magnitude
  double pmag = sqrt((pr * pr) + (ptheta * ptheta) + (pphi * pphi));
  // ||p|| / (m*c)
  double pm_ratio = pmag / (mass_ * constants::SPEED_OF_LIGHT);
  // Lorentz factor
  double gamma = sqrt(1. + (pm_ratio * pm_ratio));
  return gamma;
};  // lorentz_factor
