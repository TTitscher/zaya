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

    inline std::vector<int> Ids(Eigen::Vector3d x, double r) const
    {
        const int xs = floor((x.x() - r) * _n.x() / _l.x());
        const int ys = floor((x.y() - r) * _n.y() / _l.y());
        const int zs = floor((x.z() - r) * _n.z() / _l.z());

        const int xe = floor((x.x() + r) * _n.x() / _l.x());
        const int ye = floor((x.y() + r) * _n.y() / _l.y());
        const int ze = floor((x.z() + r) * _n.z() / _l.z());

        std::vector<int> ids;
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
                    // std::cout << xmod << " " << ymod << " " << zmod << std::endl;
                    ids.push_back(Id(xmod, ymod, zmod));
                }
            }
        }


        return ids;
    }
    void Add(int id, Eigen::Vector3d x, double r)
    {
        for (int box_id : Ids(x, r))
            _boxes.at(box_id).push_back(id);
    }

    std::vector<int> Neighbors(Eigen::Vector3d x, double r, int to_del) const
    {
        std::vector<int> n;
        // n.reserve(20);

        for (int id : Ids(x, r))
        {
            const auto& box = _boxes.at(id);
            int size = n.size();
            n.insert(n.end(), box.begin(), box.end());
            std::inplace_merge(n.begin(), n.begin() + size, n.end());
        }
        n.erase(std::unique(n.begin(), n.end()), n.end());
        n.erase(std::remove(n.begin(), n.end(), to_del), n.end());
        return n;
    }

    Eigen::Vector3i _n;
    const Eigen::Vector3d _l;
    std::vector<std::vector<int>> _boxes;
};

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

Eigen::Vector3d PeriodicDistance(const Eigen::Vector3d& xi, const Eigen::Vector3d& xj, const Eigen::Vector3d& box)
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

double DxSphere(const RowMatrixXd& x, const Eigen::VectorXd& r, const Eigen::Vector3d box, double dr, double rho,
                Eigen::Ref<RowMatrixXd> dx)
{
    dx.setZero();
    double dr_in = 1000.;

    SubBoxes boxes(box, r.maxCoeff() + dr);
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

            auto neighbors = boxes.Neighbors(x.row(i), ri + dr, i);

            const auto xj = x(neighbors, Eigen::indexing::all);
            const Eigen::VectorXd rj = r(neighbors);

            const double r_out_i = ri + dr;
            const Eigen::VectorXd r_out_j = rj.array() + dr;
            const Eigen::VectorXd allowed_distances = rj.array() + ri;
            const Eigen::VectorXd sigma = r_out_i + r_out_j.array();

            RowMatrixXd rjis(neighbors.size(), 3);
            for (int k = 0; k < 3; ++k)
            {
                Eigen::VectorXd rji_unperiodic_k = x(i, k) - xj.col(k).array();
                Eigen::VectorXd correction_k = box[k] * Eigen::round(rji_unperiodic_k.array() / box[k]);
                rjis.col(k) = rji_unperiodic_k - correction_k;
            }
            const Eigen::VectorXd abs_rji2s = rjis.rowwise().squaredNorm();
            const Eigen::VectorXd abs_rjis = abs_rji2s.cwiseSqrt();


            min_dr = std::min(min_dr, 0.5 * (abs_rjis - allowed_distances).minCoeff());

            const Eigen::VectorXd inv_sigma2s = 1 / (sigma.array() * sigma.array());

            const Eigen::VectorXd p_ijs =
                    (4 * r_out_i * r_out_j.array() * (1 - abs_rjis.array() * abs_rjis.array() * inv_sigma2s.array()) *
                     inv_sigma2s.array())
                            .cwiseMax(0.);

            dxi[0] = (rho * p_ijs.array() / ri * rjis.col(0).array() / abs_rjis.array()).sum();
            dxi[1] = (rho * p_ijs.array() / ri * rjis.col(1).array() / abs_rjis.array()).sum();
            dxi[2] = (rho * p_ijs.array() / ri * rjis.col(2).array() / abs_rjis.array()).sum();

            // for (int k = 0; k < neighbors.size(); ++k)
            //// for (int j : neighbors)
            //{
            //    int j = neighbors[k];
            //
            //    Eigen::Vector3d rji = rjis.row(k);
            //    const double abs_rji2 = abs_rji2s[k];
            //
            //    // if (abs_rji2 > sigma[k] * sigma[k])
            //    // continue;
            //
            //    const double abs_rji = abs_rjis[k];
            //
            //
            //    const double inv_sigma2 = inv_sigma2s[k];
            //    const double p_ij = p_ijs[k];
            //    dxi += rho * p_ij / ri * rji / abs_rji;
            //}
        }
#pragma omp critical
        dr_in = std::min(dr_in, min_dr);
    }
    return dr_in;
}

RowMatrixXd RSA(Eigen::VectorXd r, Eigen::Vector3d box, int max_tries = 1e5)
{
    RowMatrixXd x(r.rows(), 3);

    // The emptier the box, the more likely a large particle will fit.
    // So we add them in descending volume.
    if (not std::is_sorted(r.data(), r.data() + r.size(), std::greater<double>()))
        throw 0;
    SubBoxes boxes(box, r[0]);

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
            xi.array() += 0.5;
            xi.array() *= box.array();

            bool overlap = false;
            for (int j : boxes.Neighbors(xi, ri, i))
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
                boxes.Add(i, xi, ri);
                break;
            }

            if (i_try >= max_tries)
                throw 5.;
        }
    }
    return x;
}

} // namespace helper
