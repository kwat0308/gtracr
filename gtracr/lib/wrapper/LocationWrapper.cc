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
        .def(py::init<const double &, const double &, const double &>())
        .def(py::init<const std::string &, const double &, const double &, const double &>())
        .def_property("name", &Location::name, &Location::set_name)
        .def_property("latitude", &Location::latitude, &Location::set_latitude)
        .def_property("longitude", &Location::longitude, &Location::set_longitude)
        .def_property("altitude", &Location::altitude, &Location::set_altitude);
}