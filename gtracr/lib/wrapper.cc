#include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
#include "Location.h" // Location header file
#include "MagneticField.h"
#include "Particle.h"         // Particle header file
#include "TrajectoryTracer.h" // TrajectoryTracer header file
#include "uTrajectoryTracer.h" // uTrajectoryTracer header file
#include "pybind11/stl.h"     // for STL container type conversions

namespace py = pybind11;

PYBIND11_MODULE(_gtracr, M) {
  // matrix class
  //   M.doc() = "Evaluates the trajectory of a particle with a given mass and
  //   charge using a 4th-order Runge Kutta algorithm.";
  py::class_<TrajectoryTracer>(M, "TrajectoryTracer")
      .def(py::init<>())
      .def(py::init<const int, const double &, const double &, const double &,
                    const int, const char>())
      // .def_property("charge", &TrajectoryTracer::charge,
      // &TrajectoryTracer::set_charge) .def_property("mass",
      // &TrajectoryTracer::mass, &TrajectoryTracer::set_mass)
      // .def_property("escape_radius", &TrajectoryTracer::escape_radius,
      // &TrajectoryTracer::set_escape_radius) .def_property("step_size",
      // &TrajectoryTracer::stepsize, &TrajectoryTracer::set_stepsize)
      // .def_property("max_iter", &TrajectoryTracer::max_iter,
      // &TrajectoryTracer::set_max_iter)
      .def_property_readonly("charge", &TrajectoryTracer::charge)
      .def_property_readonly("mass", &TrajectoryTracer::mass)
      .def_property_readonly("escape_radius", &TrajectoryTracer::escape_radius)
      .def_property_readonly("step_size", &TrajectoryTracer::stepsize)
      .def_property_readonly("max_iter", &TrajectoryTracer::max_iter)
      .def_property_readonly("particle_escaped",
                             &TrajectoryTracer::particle_escaped)
      .def("evaluate", &TrajectoryTracer::evaluate)
      .def("evaluate_and_get_trajectory",
           &TrajectoryTracer::evaluate_and_get_trajectory),

  py::class_<uTrajectoryTracer>(M, "uTrajectoryTracer")
      .def(py::init<>())
      .def(py::init<const int, const double &, const double &, const double &,
                    const int, const char>())
      // .def_property("charge", &uTrajectoryTracer::charge,
      // &uTrajectoryTracer::set_charge) .def_property("mass",
      // &uTrajectoryTracer::mass, &uTrajectoryTracer::set_mass)
      // .def_property("escape_radius", &uTrajectoryTracer::escape_radius,
      // &uTrajectoryTracer::set_escape_radius) .def_property("step_size",
      // &uTrajectoryTracer::stepsize, &uTrajectoryTracer::set_stepsize)
      // .def_property("max_iter", &uTrajectoryTracer::max_iter,
      // &uTrajectoryTracer::set_max_iter)
      .def_property_readonly("charge", &uTrajectoryTracer::charge)
      .def_property_readonly("mass", &uTrajectoryTracer::mass)
      .def_property_readonly("escape_radius", &uTrajectoryTracer::escape_radius)
      .def_property_readonly("step_size", &uTrajectoryTracer::stepsize)
      .def_property_readonly("max_iter", &uTrajectoryTracer::max_iter)
      .def_property_readonly("particle_escaped",
                             &uTrajectoryTracer::particle_escaped)
      .def("evaluate", &uTrajectoryTracer::evaluate)
      .def("evaluate_and_get_trajectory",
           &uTrajectoryTracer::evaluate_and_get_trajectory);

  // py::class_<Particle>(M, "Particle")
  //     .def(py::init<>())
  //     .def(py::init<const std::string &, const int, const double &, const
  //     int, const std::string &>()) .def(py::init<const std::string &, const
  //     int, const double &, const int, const std::string &, const double &,
  //     const double &>()) .def_property_readonly("name", &Particle::name)
  //     .def_property_readonly("mass", &Particle::mass)
  //     .def_property_readonly("charge", &Particle::charge)
  //     .def_property_readonly("pdgid", &Particle::pdgid)
  //     .def_property_readonly("label", &Particle::label)
  //     .def_property("momentum", &Particle::momentum, &Particle::set_momentum)
  //     .def_property("velocity", &Particle::velocity, &Particle::set_velocity)
  //     .def_property("rigidity", &Particle::rigidity, &Particle::set_rigidity)
  //     .def("set_from_energy", &Particle::set_from_energy)
  //     .def("set_from_momentum", &Particle::set_from_momentum)
  //     .def("set_from_rigidity", &Particle::set_from_rigidity)
  //     .def("set_from_velocity", &Particle::set_from_velocity)
  //     .def("get_energy_rigidity", &Particle::get_energy_rigidity)
  //     .def("print", &Particle::print);
  // py::class_<Location>(M, "Location")
  //     .def(py::init<>())
  //     .def(py::init<const double &, const double &, const double &>())
  //     .def(py::init<const std::string &, const double &, const double &,
  //     const double &>()) .def_property_readonly("name", &Location::name)
  //     .def_property_readonly("latitude", &Location::latitude)
  //     .def_property_readonly("longitude", &Location::longitude)
  //     .def_property_readonly("altitude", &Location::altitude);
}