#pragma once
#include <algorithm>
#include <eigen3/Eigen/Core>
#include <iostream>

namespace helper
{

class SubBoxes
{
public:
    SubBoxes(double l, double r)
        : _l(l)
    {
        Init(r);
    }

    bool Init(double r)
    {
        int n = std::floor(_l / (2 * r));
        if (_n == n)
            return false;

        if (n > 200)
            throw std::runtime_error("n (" + std::to_string(n) + ") > 200. Too big.");
        _boxes.resize(n * n * n, {});
        _n = n;
        return true;
    }


    inline int Id(int i, int j, int k) const
    {
        return i * _n * _n + j * _n + k;
    }
    inline int Idk(double x) const
    {
        return floor(x * _n / _l);
    }

    inline std::vector<int> Ids(Eigen::Vector3d x, double r) const
    {
        const int xs = Idk(x.x() - r);
        const int ys = Idk(x.y() - r);
        const int zs = Idk(x.z() - r);

        const int xe = Idk(x.x() + r);
        const int ye = Idk(x.y() + r);
        const int ze = Idk(x.z() + r);

        std::vector<int> ids;
        for (int x = xs; x <= xe; ++x)
        {
            for (int y = ys; y <= ye; ++y)
            {
                for (int z = zs; z <= ze; ++z)
                {
                    // wrap coords around _n
                    const int xmod = (x + _n) % _n;
                    const int ymod = (y + _n) % _n;
                    const int zmod = (z + _n) % _n;
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

    int _n = 0;
    const double _l;
    std::vector<std::vector<int>> _boxes;
};

} // namespace helper
