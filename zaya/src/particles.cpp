#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <iostream>

#include <Eigen/Core>
#include "helper.h"

namespace py = pybind11;

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

Eigen::Vector3d PeriodicDistance(const Eigen::Vector3d& xi, const Eigen::Vector3d& xj)
{
    const Eigen::Vector3d rji_unperiodic = xi - xj;
    // we need to account for the periodicity of the structure
    // and also consider the "ghost" neighbors, i.e. this situation
    //
    //    +-----------------+
    //    |                 |
    //    |->O          O>--|
    //    |                 |
    //    +-----------------+
    //
    // where the vector through the boundary can be shorter than the
    // vector within the box.
    const Eigen::Vector3d correction = Eigen::round(rji_unperiodic.array());
    return rji_unperiodic - correction;
}

double DxSphere(const RowMatrixXd& x, const Eigen::VectorXd& r, double dr, double rho, Eigen::Ref<RowMatrixXd> dx)
{
    dx.setZero();
    double dr_in = 1000.;

    helper::SubBoxes boxes(1., r.maxCoeff() + dr);
    for (int i = 0; i < r.rows(); ++i)
        boxes.Add(i, x.row(i), r[i] + dr);

#pragma omp parallel shared(dx, r, x, boxes)
    {
        double min_dr = dr_in;

#pragma omp for
        for (int i = 0; i < r.rows(); ++i)
        {
            const double ri = r[i];
            auto dxi = dx.row(i);

            for (int j : boxes.Neighbors(x.row(i), ri + dr, i))
            {
                const double r_out_i = ri + dr;
                const double r_out_j = r[j] + dr;
                const double sigma = r_out_i + r_out_j;

                Eigen::Vector3d rji = PeriodicDistance(x.row(i), x.row(j));
                const double abs_rji2 = rji.squaredNorm();

                if (abs_rji2 > sigma * sigma)
                    continue;

                const double abs_rji = sqrt(abs_rji2);

                const double allowed_distance = ri + r[j];
                min_dr = std::min(min_dr, 0.5 * (abs_rji - allowed_distance));

                const double inv_sigma2 = 1 / sigma / sigma;
                const double p_ij = 4 * r_out_i * r_out_j * (1 - abs_rji * abs_rji * inv_sigma2) * inv_sigma2;
                dxi += rho * p_ij / ri * rji / abs_rji;
            }
        }
#pragma omp critical
        dr_in = std::min(dr_in, min_dr);
    }
    return dr_in;
}

RowMatrixXd RSA(Eigen::VectorXd r, int max_tries = 1e5)
{
    RowMatrixXd x(r.rows(), 3);

    // The emptier the box, the more likely a large particle will fit.
    // So we add them in descending volume.
    if (not std::is_sorted(r.data(), r.data() + r.size(), std::greater<double>()))
        throw 0;
    helper::SubBoxes boxes(1., r[0]);

    for (int i = 0; i < r.rows(); ++i)
    {
        const double ri = r[i];
        // update the subboxes, if they changed
        if (boxes.Init(ri))
            for (int j = 0; j < i; ++j)
                boxes.Add(j, x.row(j), r[j]);

        int i_try = 0;
        while (true)
        {
            i_try++;

            Eigen::Vector3d xi = 0.5 * Eigen::Vector3d::Random(); // -0.5, 0.5
            xi.array() += 0.5; // 0, 1

            bool overlap = false;
            for (int j : boxes.Neighbors(xi, ri, i))
            {
                const double abs_rji = PeriodicDistance(xi, x.row(j)).norm();

                if (abs_rji < (ri + r[j]))
                {
                    overlap = true;
                    break;
                }
            }
            if (not overlap)
            {
                x.row(i) = xi;
                boxes.Add(i, xi, ri);
                break;
            }

            if (i_try >= max_tries)
                throw 5.;
        }
    }
    return x;
}

PYBIND11_MODULE(_particles, m)
{
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("dx_sphere", &DxSphere, "Compute sphere movement of FBA algorithm");
    m.def("rsa", &RSA, "Random sequential addition algorithm");
}
