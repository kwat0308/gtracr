#include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
#include "pybind11/stl.h" // for STL container type conversions
#include "constants.h"
#include "Matrix.h"
#include "Vector.h"
#include "Particle.h"
#include "Location.h"
#include "RungeKutta.h"
#include "TrajectoryPoint.h"
#include "Trajectory.h"

namespace py = pybind11;

PYBIND11_MODULE(Trajectory, M)
{
  // matrix class
  M.doc() = "Trajectory class that manages a single particle trajectory";
  py::class_<Trajectory>(M, "Trajectory")
      .def(py::init<>())
      .def(py::init<const std::string &, const std::string &, const double &, const double &, const double &, const double &, const double &, const int>(), py::arg("plabel"), py::arg("locname"), py::arg("zenithAngle")=0.0, py::arg("azimuthAngle")=0.0, py::arg("energy") = 0.0, py::arg("rigidity") = 0.0, py::arg("escapeAltitude") = 1000., py::arg("maxBuffer") = 100000)
      .def(py::init<const std::string &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const int>(), py::arg("plabel"), py::arg("latitude")=0.0, py::arg("longitude")=0.0, py::arg("altitude")=0.0, py::arg("zenithAngle")=0.0, py::arg("azimuthAngle")=0.0, py::arg("energy") = 0.0, py::arg("rigidity") = 0.0, py::arg("escapeAltitude") = 1000., py::arg("maxBuffer") = 100000)
      .def("getTrajectory", &Trajectory::getTrajectory, py::arg("maxStep")=10000, py::arg("stepSize")=0.01)
      .def("getPlottingVariables", &Trajectory::getPlottingVariables)
      .def("print", &Trajectory::print);
  //   .def_property("stepSize", &Trajectory::stepSize, &Trajectory::set_stepSize)
  //   .def("evaluate", &Trajectory::evaluate);
}