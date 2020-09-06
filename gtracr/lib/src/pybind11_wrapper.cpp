
#include "MagneticField.hpp"
#include "TrajectoryTracer.hpp"  // TrajectoryTracer header file
#include "igrf.hpp"              // IGRF model
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"         // for STL container type conversions
#include "uTrajectoryTracer.hpp"  // uTrajectoryTracer header file

namespace py = pybind11;

PYBIND11_MODULE(_libgtracr, M) {
  /*
  Extension module for gtracr to C++.
  */
  // M.def("Extension module for gtracr to C++.");
  py::class_<TrajectoryTracer>(M, "TrajectoryTracer",
                               py::module_local())  // TrajectoryTracer class
      .def(py::init<>())
      .def(py::init<const int, const double &, const double &, const double &,
                    const int, const char,
                    const std::pair<std::string, double> &>())
      .def_property_readonly("charge", &TrajectoryTracer::charge)
      .def_property_readonly("mass", &TrajectoryTracer::mass)
      .def_property_readonly("escape_radius", &TrajectoryTracer::escape_radius)
      .def_property_readonly("step_size", &TrajectoryTracer::stepsize)
      .def_property_readonly("max_iter", &TrajectoryTracer::max_iter)
      .def_property_readonly("particle_escaped",
                             &TrajectoryTracer::particle_escaped)
     .def_property_readonly("final_time",
                             &TrajectoryTracer::final_time)
     .def_property_readonly("final_sixvector",
                             &TrajectoryTracer::final_sixvector)
      .def("evaluate", &TrajectoryTracer::evaluate)
      .def("evaluate_and_get_trajectory",
           &TrajectoryTracer::evaluate_and_get_trajectory),

      py::class_<uTrajectoryTracer>(
          M, "uTrajectoryTracer",
          py::module_local())  // TrajectoryTracer class unvectorized
          .def(py::init<>())
          .def(py::init<const int, const double &, const double &,
                        const double &, const int, const char,
                        const std::pair<std::string, double> &>())
          .def_property_readonly("charge", &uTrajectoryTracer::charge)
          .def_property_readonly("mass", &uTrajectoryTracer::mass)
          .def_property_readonly("escape_radius",
                                 &uTrajectoryTracer::escape_radius)
          .def_property_readonly("step_size", &uTrajectoryTracer::stepsize)
          .def_property_readonly("max_iter", &uTrajectoryTracer::max_iter)
          .def_property_readonly("particle_escaped",
                                 &uTrajectoryTracer::particle_escaped)
          .def_property_readonly("final_time",
                             &uTrajectoryTracer::final_time)
          .def_property_readonly("final_sixvector",
                             &uTrajectoryTracer::final_sixvector)
          .def("evaluate", &uTrajectoryTracer::evaluate)
          .def("evaluate_and_get_trajectory",
               &uTrajectoryTracer::evaluate_and_get_trajectory),

     py::class_<MagneticField>(M, "MagneticField", py::module_local())  // Dipole Field class
          .def(py::init<>())
          .def("values", &MagneticField::values),

      py::class_<IGRF>(M, "IGRF", py::module_local())  // IGRF class
          .def(py::init<const std::string &, const double>())
          .def_property_readonly("sdate", &IGRF::sdate)
          .def_property_readonly("nmax", &IGRF::nmax)
          .def_property_readonly("cartesian_values", &IGRF::cartesian_values)
          .def("values", &IGRF::values);

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