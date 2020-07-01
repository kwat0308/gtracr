#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "pybind11/stl.h" // for STL container type conversions
#include "RungeKutta.h" // RungeKutta header file
#include "MagneticField.h"

namespace py = pybind11;

PYBIND11_MODULE(RungeKutta, M)
{
  // matrix class
  M.doc() = "4th-order Runge Kutta Integrator";
  py::class_<RungeKutta>(M, "RungeKutta", py::buffer_protocol())
      .def(py::init<const int, const double &>())
      .def(py::init<const int, const double &, const double&>())
      .def("evaluate", &RungeKutta::evaluate);
}