#pragma once
#include <algorithm>
#include <Eigen/Core>
#include <iostream>

namespace helper
{

class SubBoxes
{
public:
    SubBoxes(Eigen::Vector3d l, double r)
        : _l(l)
    {
        Init(r);
    }

    bool Init(double r)
    {
        Eigen::Vector3i n = Eigen::floor(_l.array() / (2 * r)).cast<int>();
        if (_n == n)
            return false;

        if (n.prod() > 1e6)
            throw std::runtime_error("Would require > 1e6 boxes. Stopping.");
        _boxes.resize(n.prod(), {});
        _n = n;
        return true;
    }


    inline int Id(int i, int j, int k) const
    {
        return i * _n.y() * _n.z() + j * _n.z() + k;
    }

    inline std::tuple<int, int, int> Idk(Eigen::Vector3d x) const
    {
        return {floor(x.x() * _n.x() / _l.x()), floor(x.y() * _n.y() / _l.y()), floor(x.z() * _n.z() / _l.z())};
    }

    void Add(int id, Eigen::Vector3d x)
    {
        auto [idx, idy, idz] = Idk(x);
        const int box_id = Id(idx, idy, idz);
        _boxes.at(box_id).push_back(id);
        // int s = _boxes.at(box_id).size();
        // if (s > 1)
        // std::coumax< s << std::endl;
    }

    void Neighbors(Eigen::Vector3d x, double r, int to_del, std::vector<int>& n)
    {
        n.clear();

        auto [xs, ys, zs] = Idk(x.array() - r);
        auto [xe, ye, ze] = Idk(x.array() + r);

        for (int x = xs; x <= xe; ++x)
        {
            for (int y = ys; y <= ye; ++y)
            {
                for (int z = zs; z <= ze; ++z)
                {
                    // wrap coords around _n
                    const int xmod = (x + _n.x()) % _n.x();
                    const int ymod = (y + _n.y()) % _n.y();
                    const int zmod = (z + _n.z()) % _n.z();

                    const int box_id = Id(xmod, ymod, zmod);
                    for (const auto& box : _boxes.at(box_id))
                        if (box != to_del)
                            n.push_back(box);
                }
            }
        }
    }

    Eigen::Vector3i _n;
    const Eigen::Vector3d _l;
    std::vector<std::vector<int>> _boxes;
};

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

inline Eigen::Vector3d PeriodicDistance(const Eigen::Vector3d& xi, const Eigen::Vector3d& xj,
                                        const Eigen::Vector3d& box)
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
    const Eigen::Vector3d correction = box.array() * Eigen::round(rji_unperiodic.array() / box.array());
    return rji_unperiodic - correction;
}
class FBA
{
public:
    FBA(Eigen::Ref<RowMatrixXd> x, const Eigen::VectorXd& r, Eigen::Vector3d box, double dr0)
        : _x(x)
        , _r(r)
        , _box(box)
        , _boxes(box, r.maxCoeff() + dr0)
    {
        for (int i = 0; i < _r.rows(); ++i)
            _boxes.Add(i, x.row(i));
    }

    double Step(Eigen::Ref<RowMatrixXd> dx, double dr, double rho)
    {
        double dr_in = 1000.;
        const double rmax = _r[0];

        _boxes.Init(rmax + dr);
        for (auto& box : _boxes._boxes)
            box.clear();

        for (int i = 0; i < _r.rows(); ++i)
            _boxes.Add(i, _x.row(i));

#pragma omp parallel
        {
            double min_dr = dr_in;
            int my_n_neighbors = 0;
            std::vector<int> neighbors;
            neighbors.reserve(100);

#pragma omp for
            for (int i = 0; i < _r.rows(); ++i)
            {
                const double ri = _r[i];
                Eigen::Vector3d dxi = Eigen::Vector3d::Zero();
                const double r_out_i = ri + dr;

                _boxes.Neighbors(_x.row(i), ri + rmax + 2 * dr, i, neighbors);
                my_n_neighbors += neighbors.size();

                for (int j : neighbors)
                {
                    const double r_out_j = _r[j] + dr;
                    const double sigma = r_out_i + r_out_j;
                    const double sigma2 = sigma * sigma;
                    Eigen::Vector3d rji = PeriodicDistance(_x.row(i), _x.row(j), _box);
                    const double abs_rji2 = rji.squaredNorm();


                    if (abs_rji2 > sigma2)
                        continue;

                    const double abs_rji = std::sqrt(abs_rji2);

                    const double allowed_distance = _r[j] + ri;

                    min_dr = std::min(min_dr, (abs_rji - allowed_distance) / 2.);

                    const double inv_sigma2 = 1 / sigma2;
                    const double p_ij = 4 * r_out_i * r_out_j * (1 - abs_rji2 * inv_sigma2) * inv_sigma2;
                    dxi += p_ij / ri * rji / abs_rji;
                }
                dx.row(i) = rho * dxi;
            }
#pragma omp critical
            dr_in = std::min(dr_in, min_dr);
            n_neighbors += my_n_neighbors;
        }
        return dr_in;
    }


    Eigen::Ref<RowMatrixXd> _x;
    const Eigen::VectorXd& _r;
    Eigen::Vector3d _box;
    SubBoxes _boxes;

    unsigned long n_neighbors = 0;
};

RowMatrixXd RSA(Eigen::VectorXd r, Eigen::Vector3d box, int max_tries = 1e5)
{
    RowMatrixXd x(r.rows(), 3);

    // The emptier the box, the more likely a large particle will fit.
    // So we add them in descending volume.
    if (not std::is_sorted(r.data(), r.data() + r.size(), std::greater<double>()))
        throw 0;
    SubBoxes boxes(box, r[0]);

    std::vector<int> neighbors;

    for (int i = 0; i < r.rows(); ++i)
    {
        const double ri = r[i];
        // update the subboxes, if they changed
        if (boxes.Init(ri))
            for (int j = 0; j < i; ++j)
                boxes.Add(j, x.row(j));

        int i_try = 0;
        while (true)
        {
            i_try++;

            Eigen::Vector3d xi = 0.5 * Eigen::Vector3d::Random(); // -0.5, 0.5
            xi.array() += 0.5;
            xi.array() *= box.array();

            bool overlap = false;
            boxes.Neighbors(xi, ri + r.maxCoeff(), i, neighbors);

            for (int j : neighbors)
            {
                const double abs_rji = PeriodicDistance(xi, x.row(j), box).norm();

                if (abs_rji < (ri + r[j]))
                {
                    overlap = true;
                    break;
                }
            }
            if (not overlap)
            {
                x.row(i) = xi;
                boxes.Add(i, xi);
                break;
            }

            if (i_try >= max_tries)
                throw 5.;
        }
    }
    return x;
}

} // namespace helper
