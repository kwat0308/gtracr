#include <pybind11/pybind11.h>
// #include "pybind11/stl.h" // for STL container type conversions
#include "Location.h" // Location header file

namespace py = pybind11;

PYBIND11_MODULE(Location, M)
{
    // matrix class
    M.doc() = "Location library set from geodesic coordinates.";
    py::class_<Location>(M, "Location")
        .def(py::init<>())
        .def(py::init<const char *, const double &, const double &>())
        .def(py::init<const char *, const double &, const double &, const double &>())
        .def_readonly("name", &Particle::name)
        .def_readonly("latitude", &Particle::latitude)
        .def_readonly("longitude", &Particle::longitude)
        .def_readonly("altitude", &Particle::altitude);
}