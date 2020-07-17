#include <pybind11/pybind11.h>
#include "Particle.h" // Particle header file

namespace py = pybind11;

PYBIND11_MODULE(_particle, M)
{
    M.doc() = "Particle Library with kinematic and charge properties.";
    py::class_<Particle>(M, "Particle")
        .def(py::init<>())
        .def(py::init<const std::string &, const int, const double &, const int, const std::string &>())
        .def(py::init<const std::string &, const int, const double &, const int, const std::string &, const double &, const double &>())
        .def_property_readonly("name", &Particle::name)
        .def_property_readonly("mass", &Particle::mass)
        .def_property_readonly("charge", &Particle::charge)
        .def_property_readonly("pdgid", &Particle::pdgid)
        .def_property_readonly("label", &Particle::label)
        .def_property("momentum", &Particle::momentum, &Particle::set_momentum)
        .def_property("velocity", &Particle::velocity, &Particle::set_velocity)
        .def_property("rigidity", &Particle::rigidity, &Particle::set_rigidity)
        .def("set_from_energy", &Particle::set_from_energy)
        .def("set_from_momentum", &Particle::set_from_momentum)
        .def("set_from_rigidity", &Particle::set_from_rigidity)
        .def("set_from_velocity", &Particle::set_from_velocity)
        .def("get_energy_rigidity", &Particle::get_energy_rigidity)
        .def("print", &Particle::print);  
}