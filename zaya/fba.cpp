#include <eigen3/Eigen/Dense>
#include <iostream>

template <typename T>
void Shape(const T& what) {
    std::cout << "(" << what.rows() << "," << what.cols() << ")\n";
}

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

Eigen::Vector3d PeriodicDistance(Eigen::Vector3d a, Eigen::Vector3d b) {
    auto diff = a - b;
    Eigen::Vector3d fix;
    for (unsigned k = 0; k < 3; ++k) {
        auto abs_diff = std::abs(diff[k]);
        int sign = sgn(diff[k]);
        fix[k] = (abs_diff > 1 - abs_diff) * sign;
    }
    return diff - fix;
}

double FactorIn(const Eigen::MatrixX3d& r, const Eigen::VectorXd& d) {
    double d_spheres = 1000.;
    const unsigned int N = d.size();
    for (unsigned int i = 0; i < N - 1; ++i) {
        for (unsigned int j = i + 1; j < N; ++j) {
            const double current_distance =
                PeriodicDistance(r.row(i), r.row(j)).norm();
            const double radius = (d[j] + d[i]) / 2;

            const double factor = current_distance / radius;
            d_spheres = std::min(d_spheres, factor);
        }
    }
    return d_spheres;
}

int main(int argc, char* argv[]) {
    const double L = 1.;
    int N = std::stoi(argv[1]);
    // int N = 100;

    Eigen::Vector3d a(0.9, 0.1, 0.5);
    Eigen::Vector3d b(0.2, 0.8, 0.5);

    std::cout << PeriodicDistance(a, b) << std::endl;
    // Eigen::MatrixX3d r(2, 3);
    // r << 0.9, 0.1, 0.5, 0.2, 0.8, 0.5;
    // std::cout << r << std::endl;

    Eigen::MatrixX3d r = Eigen::MatrixX3d::Random(N, 3);
    Eigen::VectorXd d = Eigen::VectorXd::Ones(r.rows()) * 0.1;

    std::cout << FactorIn(r, d) << std::endl;

    return 0;
}
