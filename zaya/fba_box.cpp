#include <iostream>
#include <vector>
#include <chrono>

#include "src/helper.h"

#define NDEBUG

double DrOut(Eigen::VectorXd r, Eigen::VectorXd box, double eta = 1.)
{
    const double V_box = box.prod();
    return std::cbrt(V_box * eta * 3. / 4. / M_PI / r.rows()) - r[0];
}

double Volume(const Eigen::VectorXd& r, double dr)
{
    double r3 = 0.;
    for (int i = 0; i < r.rows(); ++i)
    {
        r3 += std::pow(r[i] + dr, 3.);
    }
    return 4. / 3. * M_PI * r3;
}

void FBA(helper::RowMatrixXd x, const Eigen::VectorXd& r, Eigen::Vector3d box, double rho, double tau)
{
    Eigen::setNbThreads(0);
    Eigen::initParallel();


    const double dr0 = DrOut(r, box, 1.);
    double dr = dr0;

    std::cout << dr << std::endl;
    std::cout << Volume(r, dr) / box.prod() << std::endl;

    helper::RowMatrixXd dx = helper::RowMatrixXd::Zero(r.rows(), 3);
    helper::RowMatrixXd x_periodic(r.rows(), 3);

    helper::SubBoxes boxes(box, r.maxCoeff() + dr);
    for (int i = 0; i < r.rows(); ++i)
        boxes.Add(i, x.row(i), r[i] + dr);

    for (int iteration = 0; iteration < 10000; ++iteration)
    {

        double dr_in = helper::DxSphere(x, r, box, dr, rho, dx, boxes);

        const double V_real = Volume(r, dr_in) / box.prod();
        const double V_virt = Volume(r, dr) / box.prod();

        const double nu = ceil(-log10(V_virt - V_real));
        std::cout << iteration << ": " << V_real << " | " << V_virt << " | " << nu << " | " << std::endl;

        if (dr_in > dr)
            break;

        x += dx;
        for (int i : {0, 1, 2})
            x_periodic.col(i) = box[i] * Eigen::floor(x.col(i).array() / box[i]);
        x -= x_periodic;

        dr -= pow(0.5, nu) * dr0 / (2. * tau);
    }
}


int main(int argc, char* argv[])
{
    const double L = 1.;
    int N = 100;
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

    Eigen::Vector3d box(1, 2, 3);

    Eigen::VectorXd r = 0.01 * Eigen::VectorXd::Ones(N);
    std::cout << Volume(r, 0.) / box.prod() << std::endl;

    auto x = helper::RSA(r, box);


    FBA(x, r, box, 0.0001, 1000.);


    // std::cout << FactorIn(x, r) << std::endl;
    //
    // std::cout << DrSphere(x, r, f_out, 0.01) << std::endl;
    // std::cout << average_checks << std::endl;
    // std::cout << average_size << std::endl;

    return 0;
}
