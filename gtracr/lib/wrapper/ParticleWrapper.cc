#include <pybind11/pybind11.h>
#include "Particle.h" // Particle header file

namespace py = pybind11;

PYBIND11_MODULE(Particle, M)
{
    M.doc() = "Particle Library with kinematic and charge properties.";
    py::class_<Particle>(M, "Particle")
        .def(py::init<>())
        .def(py::init<const std::string&, const int, const double &, const int, const std::string&>())
        .def(py::init<const std::string&, const int, const double &, const int, const std::string&, const double &>())
        .def("name", &Particle::name)
        .def("mass", &Particle::mass)
        .def("charge", &Particle::charge)
        .def("pdgid", &Particle::pdgid)
        .def("label", &Particle::label)
        .def("momentum", &Particle::momentum)
        .def("velocity", &Particle::velocity)
        .def("rigidity", &Particle::rigidity)
        .def("set_from_energy", &Particle::set_from_energy)
        .def("set_from_momentum", &Particle::set_from_momentum)
        .def("set_from_rigidity", &Particle::set_from_rigidity)
        .def("set_from_velocity", &Particle::set_from_velocity)
        .def("get_energy_rigidity", &Particle::get_energy_rigidity);

}