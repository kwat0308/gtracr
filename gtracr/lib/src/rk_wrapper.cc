#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // for STL container type conversions

namespace py = pybind11;

// declarations
// const std::vector<double>& evaluate(const double&, const int, const double&, const double&, const double&, const double& , const double&, const double&, const double&, const double&);

// Python binding
PYBIND11_MODULE(RungeKutta, m)
{
    m.doc() = "Runge Kutta Evaluator";
    m.def("evaluate", &evaluate, "evaluates runge_kutta");
    // return m.ptr();
}