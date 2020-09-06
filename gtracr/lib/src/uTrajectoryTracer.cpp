// Runge Kutta integrator class

#include "uTrajectoryTracer.hpp"

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

uTrajectoryTracer::uTrajectoryTracer()
    : bfield_{MagneticField()},
      charge_{1. * constants::ELEMENTARY_CHARGE},
      mass_{0.938 * constants::KG_PER_GEVC2},
      escape_radius_{10. * constants::RE},
      stepsize_{1e-5},
      max_iter_{10000},
      particle_escaped_{false} {}

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
uTrajectoryTracer::uTrajectoryTracer(const int charge, const double &mass,
                                     const double &
                                         escape_radius /*= 10. * constants::RE*/
                                     ,
                                     const double &stepsize /*= 1e-5*/,
                                     const int max_iter /*= 10000*/,
                                     const char bfield_type /*= 'd'*/,
                                     const std::pair<std::string, double>
                                         &igrf_params /*=
        {"/home/keito/devel/gtracr/data",
        2020.}*/)
    : charge_{charge * constants::ELEMENTARY_CHARGE},
      mass_{mass * constants::KG_PER_GEVC2},
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
void uTrajectoryTracer::evaluate(double t0, std::array<double, 6> vec0) {
  // assign initial values from array to trajectory vector structure
  traj_vector_.t = t0;
  traj_vector_.r = vec0[0];
  traj_vector_.theta = vec0[1];
  traj_vector_.phi = vec0[2];
  traj_vector_.pr = vec0[3];
  traj_vector_.ptheta = vec0[4];
  traj_vector_.pphi = vec0[5];

  // start the integration process
  for (int i = 0; i < max_iter_; ++i) {
    // evaluate a runge kutta step
    // return the next iteration of values
    // traj_vector = rk_step(traj_vector);
    perform_rkstep();
    // break condition depending on value of r
    // this is set based on if particle has "escaped"
    // or if the particle has reached back to earth
    // i.e. an allowed or forbidden trajectory

    // an allowed trajectory
    if (traj_vector_.r > escape_radius_) {
      particle_escaped_ = true;
      // std::cout << "Allowed Trajectory!" << std::endl;
      break;
    }

    // a forbidden trajectory
    if (traj_vector_.r < constants::RE) {
      // std::cout << "Forbidden Trajectory!" << std::endl;
      break;
    }
  }
  // get final time and six-vector
  final_time_ = traj_vector_.t;
  final_sixvector_ = {traj_vector_.r, traj_vector_.theta,
          traj_vector_.phi,    traj_vector_.pr,
          traj_vector_.ptheta, traj_vector_.pphi};
}

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
uTrajectoryTracer::evaluate_and_get_trajectory(double t0,
                                               std::array<double, 6> vec0) {
  // a container that holds the variables throughout each step
  // variables are given in the order [t, r, theta, phi, pr, ptheta,
  //   pphi]
  // this can also be replaced with 7 doubles by moving the function
  //   rk_step
  // directly into the for loop
  traj_vector_.t = t0;
  traj_vector_.r = vec0[0];
  traj_vector_.theta = vec0[1];
  traj_vector_.phi = vec0[2];
  traj_vector_.pr = vec0[3];
  traj_vector_.ptheta = vec0[4];
  traj_vector_.pphi = vec0[5];

  // set up the vectors for the values on trajectory
  // in spherical coordinates
  std::vector<double> time_arr;
  std::vector<double> r_arr;
  std::vector<double> theta_arr;
  std::vector<double> phi_arr;
  std::vector<double> pr_arr;
  std::vector<double> ptheta_arr;
  std::vector<double> pphi_arr;

  // start the integration process
  for (int i = 0; i < max_iter_; ++i) {
    // append to arrays first

    time_arr.push_back(traj_vector_.t);
    r_arr.push_back(traj_vector_.r);
    theta_arr.push_back(traj_vector_.theta);
    phi_arr.push_back(traj_vector_.phi);
    pr_arr.push_back(traj_vector_.pr);
    ptheta_arr.push_back(traj_vector_.ptheta);
    pphi_arr.push_back(traj_vector_.pphi);

    // evaluate a runge kutta step
    // return the next iteration of values
    // traj_vector = rk_step(traj_vector);
    perform_rkstep();

    // break condition depending on value of r
    // this is set based on if particle has "escaped"
    // or if the particle has reached back to earth
    // i.e. an allowed or forbidden trajectory

    // an allowed trajectory
    if (traj_vector_.r > escape_radius_) {
      particle_escaped_ = true;
      // std::cout << "Allowed Trajectory!" << std::endl;
      break;
    }

    // a forbidden trajectory
    if (traj_vector_.r < constants::RE) {
      // std::cout << "Forbidden Trajectory!" << std::endl;
      break;
    }
  }

  // get final time and six-vector
  final_time_ = traj_vector_.t;
  final_sixvector_ = {traj_vector_.r, traj_vector_.theta,
          traj_vector_.phi,    traj_vector_.pr,
          traj_vector_.ptheta, traj_vector_.pphi};

  // convert the final values of the trajectory into a std::vector
  // to put this into our map
  // dont want the time component, so start from 2nd component of
  // trajectory vector
  // std::vector<double> final_values{traj_vector_.r,      traj_vector_.theta,
  //                                  traj_vector_.phi,    traj_vector_.pr,
  //                                  traj_vector_.ptheta, traj_vector_.pphi};

  std::map<std::string, std::vector<double>> trajectory_data = {
      {"t", time_arr},      {"r", r_arr},
      {"theta", theta_arr}, {"phi", phi_arr},
      {"pr", pr_arr},       {"ptheta", ptheta_arr},
      {"pphi", pphi_arr}};

  return trajectory_data;
}

// ODEs based on relativistic Lorentz force equation with auxiliary terms
// acceleration in spherical coordinates drdt
inline double uTrajectoryTracer::dr_dt(double pr) { return pr; }

// dtheta dt
inline double uTrajectoryTracer::dtheta_dt(double r, double ptheta) {
  return ptheta / r;
}
// dphidt
inline double uTrajectoryTracer::dphi_dt(double r, double theta, double pphi) {
  return pphi / (r * sin(theta));
}
// dvrdt
inline double uTrajectoryTracer::dpr_dt(double r, double theta, double phi,
                                        double pr, double ptheta, double pphi) {
  double lorentz_term =
      (-1. * charge_) * (ptheta * bfield_.values(r, theta, phi)[2] -
                         bfield_.values(r, theta, phi)[1] * pphi);
  double auxiliary_terms = ((ptheta * ptheta) / r) + ((pphi * pphi) / r);
  double dpr_dt = lorentz_term + auxiliary_terms;
  return dpr_dt;
}

// dpthetadt
inline double uTrajectoryTracer::dptheta_dt(double r, double theta, double phi,
                                            double pr, double ptheta,
                                            double pphi) {
  double lorentz_term = (charge_) * (bfield_.values(r, theta, phi)[2] * pr -
                                     pphi * bfield_.values(r, theta, phi)[0]);
  double auxiliary_terms =
      ((pphi * pphi * cos(theta)) / (r * sin(theta))) - ((pr * ptheta) / r);
  double dptheta_dt = lorentz_term + auxiliary_terms;
  return dptheta_dt;
}

// dpphidt
inline double uTrajectoryTracer::dpphi_dt(double r, double theta, double phi,
                                          double pr, double ptheta,
                                          double pphi) {
  double lorentz_term =
      (-1. * charge_) * (pr * bfield_.values(r, theta, phi)[1] -
                         bfield_.values(r, theta, phi)[0] * ptheta);
  double auxiliary_terms =
      ((pr * pphi) / r) + ((ptheta * pphi * cos(theta)) / (r * sin(theta)));
  double dpphi_dt = lorentz_term - auxiliary_terms;
  return dpphi_dt;
}

// lorentz factor
inline double uTrajectoryTracer::gamma(double pr, double ptheta, double pphi) {
  double momentum =
      sqrt((pr * pr) + (ptheta * ptheta) + (pphi * pphi));  // momentum
  double momentum_ratio =
      momentum / (mass_ * constants::SPEED_OF_LIGHT);  // p / mc
  double gamma = sqrt(1. + (momentum_ratio * momentum_ratio));
  return gamma;
}

void uTrajectoryTracer::perform_rkstep() {
  double t = traj_vector_.t;
  double r = traj_vector_.r;
  double theta = traj_vector_.theta;
  double phi = traj_vector_.phi;
  double pr = traj_vector_.pr;
  double ptheta = traj_vector_.ptheta;
  double pphi = traj_vector_.pphi;

  // evaluate relativistic mass here
  double rel_mass = mass_ * gamma(pr, ptheta, pphi);

  double r_k1 = stepsize_ * dr_dt(pr);
  double theta_k1 = stepsize_ * dtheta_dt(r, ptheta);
  double phi_k1 = stepsize_ * dphi_dt(r, theta, pphi);
  double pr_k1 = stepsize_ * dpr_dt(r, theta, phi, pr, ptheta, pphi);
  double ptheta_k1 = stepsize_ * dptheta_dt(r, theta, phi, pr, ptheta, pphi);
  double pphi_k1 = stepsize_ * dpphi_dt(r, theta, phi, pr, ptheta, pphi);

  double r_k2 = stepsize_ * dr_dt(pr + 0.5 * pr_k1);
  double theta_k2 =
      stepsize_ * dtheta_dt(r + 0.5 * r_k1, ptheta + 0.5 * ptheta_k1);
  double phi_k2 = stepsize_ * dphi_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1,
                                      pphi + 0.5 * pphi_k1);
  double pr_k2 =
      stepsize_ * dpr_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1,
                         phi + 0.5 * phi_k1, pr + 0.5 * pr_k1,
                         ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);
  double ptheta_k2 =
      stepsize_ * dptheta_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1,
                             phi + 0.5 * phi_k1, pr + 0.5 * pr_k1,
                             ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);
  double pphi_k2 =
      stepsize_ * dpphi_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1,
                           phi + 0.5 * phi_k1, pr + 0.5 * pr_k1,
                           ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);

  double r_k3 = stepsize_ * dr_dt(pr + 0.5 * pr_k2);
  double theta_k3 =
      stepsize_ * dtheta_dt(r + 0.5 * r_k2, ptheta + 0.5 * ptheta_k2);
  double phi_k3 = stepsize_ * dphi_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2,
                                      pphi + 0.5 * pphi_k2);
  double pr_k3 =
      stepsize_ * dpr_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2,
                         phi + 0.5 * phi_k2, pr + 0.5 * pr_k2,
                         ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);
  double ptheta_k3 =
      stepsize_ * dptheta_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2,
                             phi + 0.5 * phi_k2, pr + 0.5 * pr_k2,
                             ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);
  double pphi_k3 =
      stepsize_ * dpphi_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2,
                           phi + 0.5 * phi_k2, pr + 0.5 * pr_k2,
                           ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);

  double r_k4 = stepsize_ * dr_dt(pr + pr_k3);
  double theta_k4 = stepsize_ * dtheta_dt(r + r_k3, ptheta + ptheta_k3);
  double phi_k4 =
      stepsize_ * dphi_dt(r + r_k3, theta + theta_k3, pphi + pphi_k3);
  double pr_k4 =
      stepsize_ * dpr_dt(r + r_k3, theta + theta_k3, phi + phi_k3, pr + pr_k3,
                         ptheta + ptheta_k3, pphi + pphi_k3);
  double ptheta_k4 =
      stepsize_ * dptheta_dt(r + r_k3, theta + theta_k3, phi + phi_k3,
                             pr + pr_k3, ptheta + ptheta_k3, pphi + pphi_k3);
  double pphi_k4 =
      stepsize_ * dpphi_dt(r + r_k3, theta + theta_k3, phi + phi_k3, pr + pr_k3,
                           ptheta + ptheta_k3, pphi + pphi_k3);

  traj_vector_.r +=
      (1. / (6. * rel_mass)) * (r_k1 + 2. * r_k2 + 2. * r_k3 + r_k4);
  traj_vector_.theta += (1. / (6. * rel_mass)) *
                        (theta_k1 + 2. * theta_k2 + 2. * theta_k3 + theta_k4);
  traj_vector_.phi +=
      (1. / (6. * rel_mass)) * (phi_k1 + 2. * phi_k2 + 2. * phi_k3 + phi_k4);
  traj_vector_.pr +=
      (1. / (6. * rel_mass)) * (pr_k1 + 2. * pr_k2 + 2. * pr_k3 + pr_k4);
  traj_vector_.ptheta += (1. / (6. * rel_mass)) * (ptheta_k1 + 2. * ptheta_k2 +
                                                   2. * ptheta_k3 + ptheta_k4);
  traj_vector_.pphi += (1. / (6. * rel_mass)) *
                       (pphi_k1 + 2. * pphi_k2 + 2. * pphi_k3 + pphi_k4);
  traj_vector_.t += stepsize_;

  // for (double val:vec) {
  //     std::cout << val << ' ';
  // }
  // std::cout << std::endl;

  // return vec;
}
