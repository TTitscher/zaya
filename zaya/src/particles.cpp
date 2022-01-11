#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <iostream>

#include <Eigen/Core>
#include "helper.h"

namespace py = pybind11;


PYBIND11_MODULE(_particles, m)
{
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("dx_sphere", &DxSphere, "Compute sphere movement of FBA algorithm");
    m.def("rsa", &RSA, "Random sequential addition algorithm");
}
