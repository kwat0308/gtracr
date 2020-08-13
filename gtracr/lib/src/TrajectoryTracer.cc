// Runge Kutta integrator class

#include "TrajectoryTracer.h"

#include <cmath>

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "MagneticField.h"
#include "constants.h"

/*
Operator overloading between std::array
These are free functions and are not members of
the TrajectoryTracer class.
*/

/* Element-wise addition of 2 std::array<double, 6>.
    Used for compact notations when evaluating the ODE.

    Example: vec = vec1 + vec2 is the same thing as:
    for (int i=0; i<vec1.size(); ++i) {
      vec[i] = vec1[i] + vec2[i];
    }

  Parameters
  ------------
  - vec1 (std::array<double, 6>) : the first array
  - vec2 (std::array<double, 6>) : the second array
  Returns
  --------
  - vec_sum (vec2) (std::array<double, 6>) :
      the element-wise sum of vec1 and vec2
*/
inline std::array<double, 6> operator+(std::array<double, 6> lh_vec,
                                       std::array<double, 6> rh_vec) {
  std::transform(lh_vec.begin(), lh_vec.end(), rh_vec.begin(), lh_vec.begin(),
                 std::plus<double>());
  return lh_vec;
}

/* Element-wise multiplication of 2 std::array<double, 6>.
    Used for compact notations when evaluating the ODE

    Example: vec = vec1 * vec2 is the same thing as:
    for (int i=0; i<vec1.size(); ++i) {
      vec[i] = vec1[i] * vec2[i];
    }

  Parameters
  ------------
  - vec1 (std::array<double, 6>) : the first array
  - vec2 (std::array<double, 6>) : the second array
  Returns
  --------
  - vec_mult (vec2) (std::array<double, 6>) :
      the element-wise multiplication of vec1 and vec2
*/
inline std::array<double, 6> operator*(std::array<double, 6> lh_vec,
                                       std::array<double, 6> rh_vec) {
  std::transform(lh_vec.begin(), lh_vec.end(), rh_vec.begin(), lh_vec.begin(),
                 std::multiplies<double>());
  return lh_vec;
}

/* Scalar multiplication between a value and a std::array
    Used for compact notations when evaluating the ODE

    Example: vec_smult = val * vec is the same thing as:
    for (int i = 0; i < 6; ++i) {
    vec_smult[i] = val * vec[i];
  }

  Parameters
  ------------
  - val (double) : the scalar
  - vec (std::array<double, 6>) : the array
  Returns
  --------
  - vec_smult (std::array<double, 6>) :
      the scalar multiplication of val and vec
*/
inline std::array<double, 6> operator*(const double lh_val,
                                       std::array<double, 6> rh_vec) {
  std::transform(
      rh_vec.cbegin(), rh_vec.cend(), rh_vec.begin(),
      std::bind(std::multiplies<double>(), std::placeholders::_1, lh_val));
  return rh_vec;
}

// TrajectoryTracer class

// Constructor for TrajectoryTracer
// default: proton
TrajectoryTracer::TrajectoryTracer()
    : bfield_{MagneticField()}, charge_{1. * constants::ELEMENTARY_CHARGE},
      mass_{0.938 * constants::KG_PER_GEVC2}, escape_radius_{10. *
                                                             constants::RE},
      stepsize_{1e-5}, max_iter_{10000}, particle_escaped_{false} {}

// Requires the charge and mass of the particle
TrajectoryTracer::TrajectoryTracer(const int charge, const double &mass,
                                   const double &escape_radius = 10. *
                                                                 constants::RE,
                                   const double &stepsize = 1e-5,
                                   const int max_iter = 10000,
                                   const char bfield_type = 'd')
    : charge_{charge * constants::ELEMENTARY_CHARGE},
      mass_{mass * constants::KG_PER_GEVC2}, escape_radius_{escape_radius},
      stepsize_{stepsize}, max_iter_{max_iter}, particle_escaped_{false} {
  switch (bfield_type) {
  case 'd':
    bfield_ = MagneticField();
  case 'i':
    // bfield = IGRF13();
    bfield_ = MagneticField();
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
std::array<double, 6>
TrajectoryTracer::ode_lrz(const double t, const std::array<double, 6> &vec) {

  // unpack array for readability
  double r = vec[0];
  double theta = vec[1];
  double phi = vec[2];
  double pr = vec[3];
  double ptheta = vec[4];
  double pphi = vec[5];

  // get the lorentz factor
  // momentum magnitude
  double pmag = sqrt((pr * pr) + (ptheta * ptheta) + (pphi * pphi));
  // ||p|| / (m*c)
  double pm_ratio = pmag / (mass_ * constants::SPEED_OF_LIGHT);
  // Lorentz factor
  double gamma = sqrt(1. + (pm_ratio * pm_ratio));
  // std::cout << gamma << std::endl;
  double rel_mass = mass_ * gamma;

  // evaluate B-field
  double bf_r = bfield_.Br(r, theta, phi);
  double bf_theta = bfield_.Btheta(r, theta, phi);
  double bf_phi = bfield_.Bphi(r, theta, phi);

  // std::cout << bf_r << ' ' << bf_theta << ' ' << bf_phi << std::endl;

  // get the momentum ODE
  // Note:
  // - charge is inverted to allow backtracking
  // - xx_lrz is the lorentz term in the ODE
  // - xx_sphcmp is the auxiliary terms due to acceleration
  //   in spherical coordiantes in the ODE

  // -q * (ptheta * Bphi - Btheta * pphi)
  double dprdt_lrz = -1. * charge_ * ((ptheta * bf_phi) - (bf_theta * pphi));
  // (ptheta^2 / r) + (pphi^2 / r)
  double dprdt_sphcmp = (((ptheta * ptheta) + (pphi * pphi)) / r);
  // dprdt
  double dprdt = dprdt_lrz + dprdt_sphcmp;

  // q * (pr * Bphi - Br * pphi)
  double dpthetadt_lrz = charge_ * ((pr * bf_phi) - (bf_r * pphi));
  // (pphi^2 * cos(theta) / (r * sin(theta))) - (pr*ptheta / r)
  double dpthetadt_sphcmp =
      ((pphi * pphi * cos(theta)) / (r * sin(theta))) - ((pr * ptheta) / r);
  // dpthetadt
  double dpthetadt = dpthetadt_lrz + dpthetadt_sphcmp;

  // -q * (pr * Btheta - Br * ptheta)
  double dpphidt_lrz = -1. * charge_ * ((pr * bf_theta) - (bf_r * ptheta));
  // (pr * pphi / r) + ((ptheta * pphi * cos(theta)) / (r * sin(theta)))
  double dpphidt_sphcmp =
      ((pr * pphi) / r) + ((ptheta * pphi * cos(theta)) / (r * sin(theta)));
  // dpphidt
  double dpphidt = dpphidt_lrz - dpphidt_sphcmp;

  // create the vector form of the ODE with the position components included
  // as well
  std::array<double, 6> ode_lrz = {
      pr, (ptheta / r), (pphi / (r * sin(theta))), dprdt, dpthetadt, dpphidt};

  return (1. / rel_mass) * ode_lrz;
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
void TrajectoryTracer::evaluate(double &t0, std::array<double, 6> &vec0) {
  double h = stepsize_; // step size in shorter notation

  // set the initial conditions
  double t = t0;
  std::array<double, 6> vec = vec0;

  // start the loop
  for (int i = 0; i < max_iter_; ++i) {
    // evaluate the k-coefficients

    std::array<double, 6> k1_vec = h * ode_lrz(t, vec);
    std::array<double, 6> k2_vec =
        h * ode_lrz(t + (0.5 * h), vec + (0.5 * k1_vec));
    std::array<double, 6> k3_vec =
        h * ode_lrz(t + (0.5 * h), vec + (0.5 * k2_vec));
    std::array<double, 6> k4_vec = h * ode_lrz(t + h, vec + k3_vec);

    std::array<double, 6> k_vec =
        (1. / 6.) * (k1_vec + (2. * k2_vec) + (2. * k3_vec) + k4_vec);

    // increment by weighted sum
    vec = vec + k_vec;
    t += h; // increase time

    const double &r = vec[0]; // radius, redefine for readability

    // breaking condition
    // if particle reaches escape radius
    if (r > constants::RE + escape_radius_) {
      particle_escaped_ = true;
      break;
    }

    // breaking condition
    // if particle reaches back onto Earth's surface again
    if (r < constants::RE) {
      break;
    }
  }
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
TrajectoryTracer::evaluate_and_get_trajectory(double &t0,
                                              std::array<double, 6> &vec0) {
  double h = stepsize_; // step size in shorter notation

  // set the initial conditions
  double t = t0;
  std::array<double, 6> vec = vec0;

  // for (auto val : vec) {
  //   std::cout << val << '\t';
  // }
  // std::cout << std::endl;

  // initialize array to contain data
  std::vector<double> t_arr, r_arr, theta_arr, phi_arr, pr_arr, ptheta_arr,
      pphi_arr;

  // start the loop
  for (int i = 0; i < max_iter_; ++i) {

    // append the values first
    // how vec looks like:
    // (r, theta, phi, pr, ptheta, pphi) = vec

    // for (auto val : vec) {
    //   std::cout << val << '\t';
    // }
    // std::cout << std::endl;

    t_arr.push_back(t);
    r_arr.push_back(vec[0]);
    theta_arr.push_back(vec[1]);
    phi_arr.push_back(vec[2]);
    pr_arr.push_back(vec[3]);
    ptheta_arr.push_back(vec[4]);
    pphi_arr.push_back(vec[5]);

    // evaluate the k-coefficients

    std::array<double, 6> k1_vec = h * ode_lrz(t, vec);
    std::array<double, 6> k2_vec =
        h * ode_lrz(t + (0.5 * h), vec + (0.5 * k1_vec));
    std::array<double, 6> k3_vec =
        h * ode_lrz(t + (0.5 * h), vec + (0.5 * k2_vec));
    std::array<double, 6> k4_vec = h * ode_lrz(t + h, vec + k3_vec);

    std::array<double, 6> k_vec =
        (1. / 6.) * (k1_vec + (2. * k2_vec) + (2. * k3_vec) + k4_vec);

    // for (auto kval : k_vec) {
    //   std::cout << kval << ' ';
    // }

    // increment by weighted sum
    vec = vec + k_vec;
    t += h; // increase time

    const double &r = vec[0]; // radius, redefine for readability

    // breaking condition
    // if particle reaches escape radius
    if (r > constants::RE + escape_radius_) {
      particle_escaped_ = true;
      break;
    }

    // breaking condition
    // if particle reaches back onto Earth's surface again
    if (r < constants::RE) {
      break;
    }
  }

  // convert final six-vector from std::array into std::vector
  std::vector<double> final_vec(vec.cbegin(), vec.cend());

  // create map that contains trajectory data
  std::map<std::string, std::vector<double>> trajectory_data = {
      {"t", t_arr},         {"r", r_arr},
      {"theta", theta_arr}, {"phi", phi_arr},
      {"pr", pr_arr},       {"ptheta", ptheta_arr},
      {"pphi", pphi_arr},   {"final_vector", final_vec}};

  return trajectory_data;
}

/* Element-wise addition between itself and another
    std::array.
    Used for compact notations when evaluating the ODE

    Example: vec += vec1 is the same thing as:
    for (int i=0; i<vec.size(); ++i) {
      vec[i] = vec[i] + vec1[i];
    }

  Parameters
  ------------
  - lh_vec (std::array<double, 6>) : the array being summed to
  - rh_vec (std::array<double, 6>) : the array that will sum
  Returns
  --------
  - lh_vec (std::array<double, 6>) : the array being summed to
*/
// inline std::array<double, 6> operator+=(std::array<double, 6> lh_vec,
//                                         std::array<double, 6> rh_vec) {
//   std::transform(lh_vec.begin(), lh_vec.end(), rh_vec.begin(),
//   lh_vec.begin(),
//                  std::plus<double>());
//   // for (int i = 0; i < 6; ++i) {
//   //   vec[i] = vec[i] + vec1[i];
//   // }
//   return lh_vec;
// }

// // ODEs based on relativistic Lorentz force equation with auxiliary terms
// from
// // acceleration in spherical coordinates drdt
// inline double TrajectoryTracer::dr_dt(double pr) { return pr; }

// // dtheta dt
// inline double TrajectoryTracer::dtheta_dt(double r, double ptheta) {
//   return ptheta / r;
// }
// // dphidt
// inline double TrajectoryTracer::dphi_dt(double r, double theta, double pphi)
// {
//   return pphi / (r * sin(theta));
// }
// // dvrdt
// inline double TrajectoryTracer::dpr_dt(double r, double theta, double phi,
//                                        double pr, double ptheta, double pphi)
//                                        {
//   double lorentz_term = (-1. * charge_ * constants::ELEMENTARY_CHARGE) *
//                         (ptheta * bfield_.Bphi(r, theta, phi) -
//                          bfield_.Btheta(r, theta, phi) * pphi);
//   double auxiliary_terms = ((ptheta * ptheta) / r) + ((pphi * pphi) / r);
//   double dpr_dt = lorentz_term + auxiliary_terms;
//   return dpr_dt;
// }

// // dpthetadt
// inline double TrajectoryTracer::dptheta_dt(double r, double theta, double
// phi,
//                                            double pr, double ptheta,
//                                            double pphi) {
//   double lorentz_term =
//       (charge_ * constants::ELEMENTARY_CHARGE) *
//       (bfield_.Bphi(r, theta, phi) * pr - pphi * bfield_.Br(r, theta, phi));
//   double auxiliary_terms =
//       ((pphi * pphi * cos(theta)) / (r * sin(theta))) - ((pr * ptheta) / r);
//   double dptheta_dt = lorentz_term + auxiliary_terms;
//   return dptheta_dt;
// }

// // dpphidt
// inline double TrajectoryTracer::dpphi_dt(double r, double theta, double phi,
//                                          double pr, double ptheta,
//                                          double pphi) {
//   double lorentz_term =
//       (-1. * charge_ * constants::ELEMENTARY_CHARGE) *
//       (pr * bfield_.Btheta(r, theta, phi) - bfield_.Br(r, theta, phi) *
//       ptheta);
//   double auxiliary_terms =
//       ((pr * pphi) / r) + ((ptheta * pphi * cos(theta)) / (r * sin(theta)));
//   double dpphi_dt = lorentz_term - auxiliary_terms;
//   return dpphi_dt;
// }

// // lorentz factor
// inline double TrajectoryTracer::gamma(double pr, double ptheta, double pphi)
// {
//   double momentum =
//       sqrt((pr * pr) + (ptheta * ptheta) + (pphi * pphi)); // momentum
//       magnitude
//   double momentum_ratio =
//       momentum / (mass_ * constants::SPEED_OF_LIGHT); // p / mc
//   double gamma = sqrt(1. + (momentum_ratio * momentum_ratio));
//   return gamma;
// }

// void TrajectoryTracer::perform_rkstep() {
//   double t = traj_vector_.t;
//   double r = traj_vector_.r;
//   double theta = traj_vector_.theta;
//   double phi = traj_vector_.phi;
//   double pr = traj_vector_.pr;
//   double ptheta = traj_vector_.ptheta;
//   double pphi = traj_vector_.pphi;

//   // evaluate relativistic mass here
//   double rel_mass = mass_ * gamma(pr, ptheta, pphi) *
//   constants::KG_PER_GEVC2;

//   // std::cout << gamma(pr, ptheta, pphi) << std::endl;

//   // std::cout << t << ' ' << r << ' ' << theta << ' ' << phi << ' ' << pr <<
//   '
//   // ' << ptheta << ' ' << pphi << ' ' << std::endl;

//   double r_k1 = stepsize_ * dr_dt(pr);
//   double theta_k1 = stepsize_ * dtheta_dt(r, ptheta);
//   double phi_k1 = stepsize_ * dphi_dt(r, theta, pphi);
//   double pr_k1 = stepsize_ * dpr_dt(r, theta, phi, pr, ptheta, pphi);
//   double ptheta_k1 = stepsize_ * dptheta_dt(r, theta, phi, pr, ptheta, pphi);
//   double pphi_k1 = stepsize_ * dpphi_dt(r, theta, phi, pr, ptheta, pphi);

//   // std::cout << r_k1 << ' ' << theta_k1 << ' ' << phi_k1 << ' ' << pr_k1 <<
//   '
//   // ' << ptheta_k1 << ' ' << pphi_k1 << std::endl;

//   double r_k2 = stepsize_ * dr_dt(pr + 0.5 * pr_k1);
//   double theta_k2 =
//       stepsize_ * dtheta_dt(r + 0.5 * r_k1, ptheta + 0.5 * ptheta_k1);
//   double phi_k2 = stepsize_ * dphi_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1,
//                                       pphi + 0.5 * pphi_k1);
//   double pr_k2 =
//       stepsize_ * dpr_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1,
//                          phi + 0.5 * phi_k1, pr + 0.5 * pr_k1,
//                          ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);
//   double ptheta_k2 =
//       stepsize_ * dptheta_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1,
//                              phi + 0.5 * phi_k1, pr + 0.5 * pr_k1,
//                              ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);
//   double pphi_k2 =
//       stepsize_ * dpphi_dt(r + 0.5 * r_k1, theta + 0.5 * theta_k1,
//                            phi + 0.5 * phi_k1, pr + 0.5 * pr_k1,
//                            ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1);

//   // std::cout << r_k2 << ' ' << theta_k2 << ' ' << phi_k2 << ' ' << pr_k2 <<
//   '
//   // ' << ptheta_k2 << ' ' << pphi_k2 << std::endl;

//   double r_k3 = stepsize_ * dr_dt(pr + 0.5 * pr_k2);
//   double theta_k3 =
//       stepsize_ * dtheta_dt(r + 0.5 * r_k2, ptheta + 0.5 * ptheta_k2);
//   double phi_k3 = stepsize_ * dphi_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2,
//                                       pphi + 0.5 * pphi_k2);
//   double pr_k3 =
//       stepsize_ * dpr_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2,
//                          phi + 0.5 * phi_k2, pr + 0.5 * pr_k2,
//                          ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);
//   double ptheta_k3 =
//       stepsize_ * dptheta_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2,
//                              phi + 0.5 * phi_k2, pr + 0.5 * pr_k2,
//                              ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);
//   double pphi_k3 =
//       stepsize_ * dpphi_dt(r + 0.5 * r_k2, theta + 0.5 * theta_k2,
//                            phi + 0.5 * phi_k2, pr + 0.5 * pr_k2,
//                            ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2);

//   // std::cout << r_k3 << ' ' << theta_k3 << ' ' << phi_k3 << ' ' << pr_k3 <<
//   '
//   // ' << ptheta_k3 << ' ' << pphi_k3 << std::endl;

//   double r_k4 = stepsize_ * dr_dt(pr + pr_k3);
//   double theta_k4 = stepsize_ * dtheta_dt(r + r_k3, ptheta + ptheta_k3);
//   double phi_k4 =
//       stepsize_ * dphi_dt(r + r_k3, theta + theta_k3, pphi + pphi_k3);
//   double pr_k4 =
//       stepsize_ * dpr_dt(r + r_k3, theta + theta_k3, phi + phi_k3, pr +
//       pr_k3,
//                          ptheta + ptheta_k3, pphi + pphi_k3);
//   double ptheta_k4 =
//       stepsize_ * dptheta_dt(r + r_k3, theta + theta_k3, phi + phi_k3,
//                              pr + pr_k3, ptheta + ptheta_k3, pphi + pphi_k3);
//   double pphi_k4 =
//       stepsize_ * dpphi_dt(r + r_k3, theta + theta_k3, phi + phi_k3, pr +
//       pr_k3,
//                            ptheta + ptheta_k3, pphi + pphi_k3);

//   traj_vector_.r +=
//       (1. / (6. * rel_mass)) * (r_k1 + 2. * r_k2 + 2. * r_k3 + r_k4);
//   traj_vector_.theta += (1. / (6. * rel_mass)) *
//                         (theta_k1 + 2. * theta_k2 + 2. * theta_k3 +
//                         theta_k4);
//   traj_vector_.phi +=
//       (1. / (6. * rel_mass)) * (phi_k1 + 2. * phi_k2 + 2. * phi_k3 + phi_k4);
//   traj_vector_.pr +=
//       (1. / (6. * rel_mass)) * (pr_k1 + 2. * pr_k2 + 2. * pr_k3 + pr_k4);
//   traj_vector_.ptheta += (1. / (6. * rel_mass)) * (ptheta_k1 + 2. * ptheta_k2
//   +
//                                                    2. * ptheta_k3 +
//                                                    ptheta_k4);
//   traj_vector_.pphi += (1. / (6. * rel_mass)) *
//                        (pphi_k1 + 2. * pphi_k2 + 2. * pphi_k3 + pphi_k4);
//   traj_vector_.t += stepsize_;

//   // for (double val:vec) {
//   //     std::cout << val << ' ';
//   // }
//   // std::cout << std::endl;

//   // return vec;
// }

// // copy constructor
// TrajectoryTracer::TrajectoryTracer(const TrajectoryTracer &traj_tracer)
//     : bfield_{traj_tracer.bfield_},
//       charge_{traj_tracer.charge_},
//       mass_{traj_tracer.mass_},
//       escape_radius_{traj_tracer.escape_radius_},
//       stepsize_{traj_tracer.stepsize_},
//       max_iter_{traj_tracer.max_iter_},
//       particle_escaped_{false} {}

// // copy assignment operator
// TrajectoryTracer &TrajectoryTracer::operator=(
//     const TrajectoryTracer &traj_tracer)
// {
//   bfield_ = traj_tracer.bfield_;
//   charge_ = traj_tracer.charge_;
//   mass_ = traj_tracer.mass_;
//   escape_radius_ = traj_tracer.escape_radius_;
//   stepsize_ = traj_tracer.stepsize_;
//   max_iter_ = traj_tracer.max_iter_;
//   particle_escaped_ = false;
//   return *this;
// }

// // a container that holds the variables throughout each step
// // variables are given in the order [t, r, theta, phi, pr, ptheta, pphi]
// // this can also be replaced with 7 doubles by moving the function
// rk_step
// // directly into the for loop
// void TrajectoryTracer::evaluate(std::array<double, 7> &initial_values) {
//   // std::array<double, 7> traj_vector = initial_values;
//   // assign initial values from array to trajectory vector structure
//   traj_vector_.t = initial_values[0];
//   traj_vector_.r = initial_values[1];
//   traj_vector_.theta = initial_values[2];
//   traj_vector_.phi = initial_values[3];
//   traj_vector_.pr = initial_values[4];
//   traj_vector_.ptheta = initial_values[5];
//   traj_vector_.pphi = initial_values[6];

//   // start the integration process
//   for (int i = 0; i < max_iter_; ++i) {
//     // append to arrays first
//     // to do this we need to convert spherical to cartesian

//     // first rename the variables for readability
//     // these can probably be const but lets leave that for now
//     // double t = traj_vector[0];
//     // double r = traj_vector[1];
//     // double theta = traj_vector[2];
//     // double phi = traj_vector[3];
//     // double pr = traj_vector[4];
//     // double ptheta = traj_vector[5];
//     // double pphi = traj_vector[6];

//     // evaluate a runge kutta step
//     // return the next iteration of values
//     // traj_vector = rk_step(traj_vector);
//     perform_rkstep();
//     // break condition depending on value of r
//     // this is set based on if particle has "escaped"
//     // or if the particle has reached back to earth
//     // i.e. an allowed or forbidden trajectory

//     // const double &radius = traj_vector[1];

//     // an allowed trajectory
//     if (traj_vector_.r > constants::RE + escape_radius_) {
//       particle_escaped_ = true;
//       // std::cout << "Allowed Trajectory!" << std::endl;
//       break;
//     }

//     // a forbidden trajectory
//     if (traj_vector_.r < constants::RE) {
//       // std::cout << "Forbidden Trajectory!" << std::endl;
//       break;
//     }
//   }
// }

// // evaluate the full runge kutta algorithm (i.e. this contains the loop)
// // input is a 7-vector that contains the initial values for the
// integration
// // processs in the form [t0, r0, theta0, phi0, pr0, ptheta0, pphi0]
// std::map<std::string, std::vector<double>>
// TrajectoryTracer::evaluate_and_get_trajectories(
//     std::array<double, 7> &initial_values) {
//   // arrays to store trajectory coordinates and momentum
//   // this should change to pointers for a faster performance in the
//   future
//   // std::vector<double> time_arr;
//   // std::vector<double> x_arr;
//   // std::vector<double> y_arr;
//   // std::vector<double> z_arr;
//   // std::vector<double> px_arr;
//   // std::vector<double> py_arr;
//   // std::vector<double> pz_arr;

//   // a container that holds the variables throughout each step
//   // variables are given in the order [t, r, theta, phi, pr, ptheta,
//   pphi]
//   // this can also be replaced with 7 doubles by moving the function
//   rk_step
//   // directly into the for loop
//   // std::array<double, 7> traj_vector = initial_values;
//   traj_vector_.t = initial_values[0];
//   traj_vector_.r = initial_values[1];
//   traj_vector_.theta = initial_values[2];
//   traj_vector_.phi = initial_values[3];
//   traj_vector_.pr = initial_values[4];
//   traj_vector_.ptheta = initial_values[5];
//   traj_vector_.pphi = initial_values[6];

//   // set up the vectors for the values on trajectory
//   // in spherical coordinates
//   std::vector<double> time_arr;
//   std::vector<double> r_arr;
//   std::vector<double> theta_arr;
//   std::vector<double> phi_arr;
//   std::vector<double> pr_arr;
//   std::vector<double> ptheta_arr;
//   std::vector<double> pphi_arr;

//   // start the integration process
//   for (int i = 0; i < max_iter_; ++i) {
//     // append to arrays first
//     // to do this we need to convert spherical to cartesian

//     // first rename the variables for readability
//     // these can probably be const but lets leave that for now
//     double t = traj_vector_.t;
//     double r = traj_vector_.r;
//     double theta = traj_vector_.theta;
//     double phi = traj_vector_.phi;
//     double pr = traj_vector_.pr;
//     double ptheta = traj_vector_.ptheta;
//     double pphi = traj_vector_.pphi;

//     // convert the coordinates
//     // double x = r * sin(theta) * cos(phi);
//     // double y = r * sin(theta) * sin(phi);
//     // double z = r * cos(theta);

//     // // convert the momentum
//     // double px = pr * sin(theta) * cos(phi) +
//     //             r * ptheta * cos(theta) * cos(phi) -
//     //             r * pphi * sin(theta) * sin(phi);
//     // double py = pr * sin(theta) * sin(phi) +
//     //             r * ptheta * cos(theta) * sin(phi) +
//     //             r * pphi * sin(theta) * cos(phi);
//     // double pz = pr * cos(theta) - r * ptheta * sin(theta);

//     // // finally append the values

//     // time_arr.push_back(t);
//     // x_arr.push_back(x);
//     // y_arr.push_back(y);
//     // z_arr.push_back(z);
//     // px_arr.push_back(px);
//     // py_arr.push_back(py);
//     // pz_arr.push_back(pz);

//     time_arr.push_back(t);
//     r_arr.push_back(r);
//     theta_arr.push_back(theta);
//     phi_arr.push_back(phi);
//     pr_arr.push_back(pr);
//     ptheta_arr.push_back(ptheta);
//     pphi_arr.push_back(pphi);

//     // evaluate a runge kutta step
//     // return the next iteration of values
//     // traj_vector = rk_step(traj_vector);
//     perform_rkstep();

//     // break condition depending on value of r
//     // this is set based on if particle has "escaped"
//     // or if the particle has reached back to earth
//     // i.e. an allowed or forbidden trajectory

//     // const double &radius = traj_vector[1];

//     // an allowed trajectory
//     if (traj_vector_.r > constants::RE + escape_radius_) {
//       particle_escaped_ = true;
//       // std::cout << "Allowed Trajectory!" << std::endl;
//       break;
//     }

//     // a forbidden trajectory
//     if (traj_vector_.r < constants::RE) {
//       // std::cout << "Forbidden Trajectory!" << std::endl;
//       break;
//     }
//   }

//   // std::cout << "Total Number of iterations: " << i << std::endl;

//   // convert the final values of the trajectory into a std::vector
//   // to put this into our map
//   // dont want the time component, so start from 2nd component of
//   // trajectory vector
//   // std::vector<double> final_values(traj_vector.begin() + 1,
//   // traj_vector.end());
//   std::vector<double> final_values{traj_vector_.r, traj_vector_.theta,
//                                    traj_vector_.phi,    traj_vector_.pr,
//                                    traj_vector_.ptheta,
//                                    traj_vector_.pphi};

//   // create a map and return the map
//   //   std::map<std::string, std::vector<double>> value_map = {
//   //       {"t", time_arr}, {"x", x_arr},
//   //       {"y", y_arr},    {"z", z_arr},
//   //       {"px", px_arr},  {"py", py_arr},
//   //       {"pz", pz_arr},  {"final_values", final_values}};

//   //   return value_map;
//   // }

//   std::map<std::string, std::vector<double>> value_map = {
//       {"t", time_arr},      {"r", r_arr},
//       {"theta", theta_arr}, {"phi", phi_arr},
//       {"pr", pr_arr},       {"ptheta", ptheta_arr},
//       {"pphi", pphi_arr},   {"final_values", final_values}};

//   return value_map;
// }
