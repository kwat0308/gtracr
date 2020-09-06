/*
Trajectory Tracer class that traces the trajectory of the particle
by performing a 4th-order Runge Kutta numerical integration algorithm.
*/

#include "TrajectoryTracer.hpp"

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
TrajectoryTracer::TrajectoryTracer(const int charge, const double &mass,
                                   const double
                                       &escape_radius /*= 10. * constants::RE*/,
                                   const double &stepsize /*= 1e-5*/,
                                   const int max_iter /* = 10000*/,
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
std::array<double, 6> TrajectoryTracer::ode_lrz(
    const double t, const std::array<double, 6> &vec) {
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
  // double bf_r = bfield_.Br(r, theta, phi);
  // double bf_theta = bfield_.Btheta(r, theta, phi);
  // double bf_phi = bfield_.Bphi(r, theta, phi);
  std::array<double, 3> bf_values = bfield_.values(r, theta, phi);
  double bf_r = bf_values[0];
  double bf_theta = bf_values[1];
  double bf_phi = bf_values[2];
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
  double h = stepsize_;  // step size in shorter notation

  // set the initial conditions
  double t = t0;
  std::array<double, 6> vec = vec0;
  // TODO: make arrays into references for no copying

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
    t += h;  // increase time

    const double &r = vec[0];  // radius, redefine for readability

    // breaking condition
    // if particle reaches escape radius
    if (r > escape_radius_) {
      particle_escaped_ = true;
      break;
    }  // if (r > escape_radius_)

    // breaking condition
    // if particle reaches back onto Earth's surface again
    if (r < constants::RE) {
      break;
    }  // if (r < constants::RE)
  }  // for (int i = 0; i < max_iter_; ++i)
  // store the final time and six-vector for checking purposes
  // the last recorded time and six-vector is the final six-vector / time
  final_time_ = t;
  final_sixvector_ = vec;
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
    t += h;  // increase time

    const double &r = vec[0];  // radius, redefine for readability

    // breaking condition
    // if particle reaches escape radius
    if (r > escape_radius_) {
      particle_escaped_ = true;
      break;
    } // if (r > escape_radius_)

    // breaking condition
    // if particle reaches back onto Earth's surface again
    if (r < constants::RE) {
      break;
    } // if (r < constants::RE)
    
  } // for (int i = 0; i < max_iter_; ++i)

  // store the final time and six-vector for checking purposes
  // the last recorded time and six-vector is the final six-vector / time
  final_time_ = t;
  final_sixvector_ = vec;

  // create map that contains trajectory data
  std::map<std::string, std::vector<double>> trajectory_data = {
      {"t", t_arr},         {"r", r_arr},
      {"theta", theta_arr}, {"phi", phi_arr},
      {"pr", pr_arr},       {"ptheta", ptheta_arr},
      {"pphi", pphi_arr}};

  return trajectory_data;
} // evaluate_and_get_trajectory
