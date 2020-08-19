#include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
#include "pybind11/stl.h" // for STL container type conversions
#include "TrajectoryTracer.h"   // TrajectoryTracer header file
#include "MagneticField.h"

namespace py = pybind11;

PYBIND11_MODULE(_trajectorytracer, M)
{
  // matrix class
  M.doc() = "Evaluates the trajectory of a particle with a given mass and charge using a 4th-order Runge Kutta algorithm.";
  py::class_<TrajectoryTracer>(M, "TrajectoryTracer")
      .def(py::init<>())
      .def(py::init<const int, const double &, const double&, 
                    const double &, const int>())
      // .def_property("charge", &TrajectoryTracer::charge, &TrajectoryTracer::set_charge)
      // .def_property("mass", &TrajectoryTracer::mass, &TrajectoryTracer::set_mass)
      // .def_property("escape_radius", &TrajectoryTracer::escape_radius, &TrajectoryTracer::set_escape_radius)
      // .def_property("step_size", &TrajectoryTracer::stepsize, &TrajectoryTracer::set_stepsize)
      // .def_property("max_iter", &TrajectoryTracer::max_iter, &TrajectoryTracer::set_max_iter)
      .def_property_readonly("charge", &TrajectoryTracer::charge)
      .def_property_readonly("mass", &TrajectoryTracer::mass)
      .def_property_readonly("escape_radius", &TrajectoryTracer::escape_radius)
      .def_property_readonly("step_size", &TrajectoryTracer::stepsize)
      .def_property_readonly("max_iter", &TrajectoryTracer::max_iter)
      .def_property_readonly("particle_escaped", &TrajectoryTracer::particle_escaped)
      .def("evaluate", &TrajectoryTracer::evaluate)
      .def("evaluate_and_get_trajectories", &TrajectoryTracer::evaluate_and_get_trajectories);
}