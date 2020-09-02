// runs ~300 iterations of the TrajectoryTracer object
// to evaluate the performance of the C++ portion
// of the package
#include <array>
// #include <map>
// #include <string>
#include <vector>

#include "TrajectoryTracer.hpp"
#include "constants.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

std::map<std::string, std::vector<double>> convert_trajectory_to_cartesian(
    std::map<std::string, std::vector<double>> result_map) {
  std::vector<double> time_arr = result_map["t"];
  std::vector<double> x_arr;
  std::vector<double> y_arr;
  std::vector<double> z_arr;
  std::vector<double> px_arr;
  std::vector<double> py_arr;
  std::vector<double> pz_arr;

  for (int i = 0; i < time_arr.size(); ++i) {
    double r = result_map["r"][i] / constants::RE;
    double theta = result_map["theta"][i];
    double phi = result_map["phi"][i];
    double pr = result_map["pr"][i];
    double ptheta = result_map["ptheta"][i];
    double pphi = result_map["pphi"][i];

    // convert the coordinates
    double x = r * sin(theta) * cos(phi);
    double y = r * sin(theta) * sin(phi);
    double z = r * cos(theta);

    // convert the momentum
    double px = pr * sin(theta) * cos(phi) +
                r * ptheta * cos(theta) * cos(phi) -
                r * pphi * sin(theta) * sin(phi);
    double py = pr * sin(theta) * sin(phi) +
                r * ptheta * cos(theta) * sin(phi) +
                r * pphi * sin(theta) * cos(phi);
    double pz = pr * cos(theta) - r * ptheta * sin(theta);

    // finally append the values
    x_arr.push_back(x);
    y_arr.push_back(y);
    z_arr.push_back(z);
    px_arr.push_back(px);
    py_arr.push_back(py);
    pz_arr.push_back(pz);
  }

  std::map<std::string, std::vector<double>> value_map = {
      {"t", time_arr}, {"x", x_arr},   {"y", y_arr},  {"z", z_arr},
      {"px", px_arr},  {"py", py_arr}, {"pz", pz_arr}};

  return value_map;
}  // namespace
   //
// std::map<std::string, std::vector<double>> convert_trajectory_to_cartesian(
//     std::map<std::string, std::vector<double>> result_map)

void plot_trajectory(
    std::map<std::string, std::vector<double>> trajectory_map) {
  std::vector<double> t_arr = trajectory_map["t"];
  std::vector<double> x_arr = trajectory_map["x"];
  std::vector<double> y_arr = trajectory_map["y"];
  std::vector<double> z_arr = trajectory_map["z"];

  plt::plot(x_arr, y_arr);
  // plt::save("test.png");
  plt::show();
}

int main() {
  //   int max_iter = 1;  // max number of iterations
  //   // contain everything within a for loop
  //   // since we imitate whats done in Python but just doing it all in C++
  //   for (int i = 0; i < max_iter; ++i) {
  // initialize trajectory tracer
  // use default (proton, stepsize=1e-5, maxsteps=10000)
  TrajectoryTracer traj_tracer =
      TrajectoryTracer(1, 0.938, 10 * constants::RE, 1e-5, 10000, 'i');

  // define initial values
  double t0 = 0.;
  std::array<double, 6> vec0{constants::RE + 100000. + 100000.,  // r0
                             constants::pi / 2.,                 // theta0
                             constants::pi / 2.,                 // phi0
                             0.,                                 // pr0
                             21. * constants::KG_M_S_PER_GEVC,   // ptheta0
                             0.};                                // pphi0
  // obtain the result as a map (as it is currently)
  std::map<std::string, std::vector<double>> result_map =
      traj_tracer.evaluate_and_get_trajectory(t0, vec0);
  // // here now we just evaluate and dont get the trajectory data
  // // traj_tracer.evaluate(initial_values);

  std::map<std::string, std::vector<double>> cartesian_map =
      convert_trajectory_to_cartesian(result_map);

  plot_trajectory(cartesian_map);
}