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
        .def(py::init<const std::string &, const double &, const double &>())
        .def(py::init<const std::string &, const double &, const double &, const double &>())
        .def("name", &Location::name)
        .def("latitude", &Location::latitude)
        .def("longitude", &Location::longitude)
        .def("altitude", &Location::altitude);
}