#include <pybind11/pybind11.h>
#include "Particle.h" // Particle header file

namespace py = pybind11;

PYBIND11_MODULE(Particle, M)
{
    M.doc() = "Particle Library with kinematic and charge properties.";
    py::class_<Particle>(M, "Particle")
        .def(py::init<>())
        .def(py::init<const std::string &, const int, const double &, const int, const std::string &>())
        .def(py::init<const std::string &, const int, const double &, const int, const std::string &, const double &, const double &>())
        .def_property("name", &Particle::name, &Particle::set_name)
        .def_property("mass", &Particle::mass, &Particle::set_mass)
        .def_property("charge", &Particle::charge, &Particle::set_charge)
        .def_property("pdgid", &Particle::pdgid, &Particle::set_pdgid)
        .def_property("label", &Particle::label, &Particle::set_label)
        .def_property("momentum", &Particle::momentum, &Particle::set_momentum)
        .def_property("velocity", &Particle::velocity, &Particle::set_velocity)
        .def_property("rigidity", &Particle::rigidity, &Particle::set_rigidity)
        .def("set_from_energy", &Particle::set_from_energy)
        .def("set_from_momentum", &Particle::set_from_momentum)
        .def("set_from_rigidity", &Particle::set_from_rigidity)
        .def("set_from_velocity", &Particle::set_from_velocity)
        .def("get_energy_rigidity", &Particle::get_energy_rigidity)
        .def("print", &Particle::print);
        // .def("__repr__",
        // [](const Particle &part) {
        //     return "<Particle.Particle object of " + part.label() + " with momentum " + part.momentum() ">";
        // })
        // .def("__str__",
        // [](const Particle &part) {
        //     return part.name() + "(" + part.label() + "): PDGID: " + part.pdgid() + " Charge [e]: " + part.charge() + " Mass [GeV/c^2]: " + part.mass()
        //             + " Momentum [GeV/c]: " + part.momentum() + " Velocity [c]: " + part.velocity() + " Rigidity [GV]: " + part.rigidity() + "\n";
        // });   
}