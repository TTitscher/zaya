#include <cmath>
#include <iostream>
#include <vector>

#include "src/helper.h"

#define BOXES 1


double FactorOut(const Eigen::VectorXd& r, double eta = 1.)
{
    const double V_box = 1.;
    double d3 = 0.;
    for (int i = 0; i < r.rows(); ++i)
    {
        d3 += std::pow(2. * r[i], 3.);
    }
    return std::cbrt(V_box * eta * 6. / M_PI / d3);
}

double Volume(const Eigen::VectorXd& r, double factor)
{
    double r3 = 0.;
    for (int i = 0; i < r.rows(); ++i)
    {
        r3 += std::pow(factor * r[i], 3.);
    }
    return 4. / 3. * M_PI * r3;
}

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

void DxSphere(const RowMatrixXd& x, const Eigen::VectorXd& r, double f_out, double rho, double& f_in,
              Eigen::Ref<RowMatrixXd> dx)
{
    dx.setZero();

    helper::SubBoxes boxes(1., r.maxCoeff() * f_out);
    for (int i = 0; i < r.rows(); ++i)
        boxes.Add(i, x.row(i), r[i] * f_out);

    for (int i = 0; i < r.rows(); ++i)
    {
        const Eigen::Vector3d& xi = x.row(i);
        const double ri = r[i];
        auto dxi = dx.row(i);

        for (int j : boxes.Neighbors(xi, ri * f_out))
        {
            if (j == i)
                continue;

            const double r_out_i = f_out * ri;
            const double r_out_j = f_out * r[j];
            const double sigma = r_out_i + r_out_j;

            const Eigen::Vector3d& xj = x.row(j);
            Eigen::Vector3d rji = xi - xj;
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

            for (int k = 0; k < 3; ++k)
            {
                auto abs_diff = std::abs(rji[k]);
                int sign = static_cast<int>((0. < rji[k]) - (rji[k] < 0.));

                rji[k] -= (abs_diff > 1 - abs_diff) * sign;
            }
            const double abs_rji = rji.norm();

            const double allowed_distance = ri + r[j];
            f_in = std::min(f_in, abs_rji / allowed_distance);

            if (abs_rji > sigma)
                continue;

            const double inv_sigma2 = 1 / sigma / sigma;
            const double p_ij = 4 * r_out_i * r_out_j * (1 - abs_rji * abs_rji * inv_sigma2) * inv_sigma2;
            dxi += rho * p_ij / ri * rji / abs_rji;
        }
    }
}

void FBA(RowMatrixXd x, const Eigen::VectorXd& r, double rho, double tau)
{
    const double f_out0 = FactorOut(r, 1.);
    double f_out = f_out0;

    RowMatrixXd dx = RowMatrixXd::Zero(r.rows(), 3);

    for (int iteration = 0; iteration < 10000; ++iteration)
    {
        // std::cout << f_out << std::endl;


        // std::cout << x.minCoeff() << std::endl;
        // std::cout << x.maxCoeff() << std::endl;
        double f_in = 1000.;
        DxSphere(x, r, f_out, rho, f_in, dx);


        const double V_real = Volume(r, f_in);
        const double V_virt = Volume(r, f_out);

        const double nu = ceil(-log10(V_virt - V_real));
        std::cout << iteration << ": " << V_real << " | " << V_virt << " | " << nu << " | " << std::endl;


        if (f_in > f_out)
            return;


        x += dx;
        for (int irow = 0; irow < r.rows(); ++irow)
        {
            for (int icol = 0; icol < 3; ++icol)
            {
                while (x(irow, icol) < 0.)
                    x(irow, icol) += 1.;
                while (x(irow, icol) >= 1.)
                    x(irow, icol) -= 1.;

                // if x(irow, icol) >= 0.);
                // assert(x(irow, icol) < 1.);
            }
        }
        // std::cout << x << std::endl;
        // const double num = 1.;
        // auto mod_matrix = (x.array() - (num * (x.array() / num))).matrix();
        // std::cout << x << std::endl;
        // x = mod_matrix;
        // auto y = x.unaryExpr(std::fmod);
        // x = y;
        f_out -= pow(0.5, nu) * f_out / (2 * tau);
    }
}


int main(int argc, char* argv[])
{
    const double L = 1.;
    int N = 10;
    try
    {
        N = std::stoi(argv[1]);
    }
    catch (std::logic_error)
    {
    }
    // int N = 100;
    // Eigen::MatrixX3d x(2, 3);
    // x << 0.9, 0.1, 0.5, 0.2, 0.8, 0.5;
    // std::cout << x << std::endl;

    RowMatrixXd x = Eigen::MatrixX3d::Random(N, 3) * 0.5;
    x.array() += 0.5;
    Eigen::VectorXd r = Eigen::VectorXd::Random(x.rows());
    r.array() += 1.1;

    FBA(x, r, 0.01, 1000);


    // std::cout << FactorIn(x, r) << std::endl;
    //
    // std::cout << DrSphere(x, r, f_out, 0.01) << std::endl;
    // std::cout << average_checks << std::endl;
    // std::cout << average_size << std::endl;

    return 0;
}
