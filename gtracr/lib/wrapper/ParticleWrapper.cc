#include <pybind11/pybind11.h>
#include "Particle.h" // Particle header file

namespace py = pybind11;

PYBIND11_MODULE(libgtracr, M)
{
    M.doc() = "Particle Library with kinematic and charge properties.";
    py::class_<Particle>(M, "Particle")
        .def(py::init<>())
        .def(py::init<const char *, const int, const double &, const int, const char *>())
        .def(py::init<const char *, const int, const double &, const int, const char *, const double &>())
        .def_readonly("name", &Particle::name)
        .def_readonly("mass", &Particle::mass)
        .def_readonly("charge", &Particle::charge)
        .def_readonly("pdgid", &Particle::pdgid)
        .def_readonly("label", &Particle::label)
        .def_readonly("momentum", &Particle::momentum)
        .def_readonly("velocity", &Particle::velocity)
        .def_readonly("rigidity", &Particle::rigidity);
}