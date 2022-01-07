#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <Eigen/Core>

namespace py = pybind11;

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

void tripple(Eigen::Ref<RowMatrixXd> x)
{
    x *= 3;
}

PYBIND11_MODULE(_particles, m)
{
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("tripple", &tripple, "Function that tripples a vector");
}
