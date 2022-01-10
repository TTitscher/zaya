#pragma once
#include <algorithm>
#include <eigen3/Eigen/Core>

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
        int n = _l / (2 * r);
        if (_n == n)
            return false;

        if (n > 200)
            throw std::runtime_error("n (" + std::to_string(n) + ") > 200. Too big.");
        _boxes.resize(n * n * n, {});
        _n = n;
        return true;
    }


    int Id(int i, int j, int k) const
    {
        return i * _n * _n + j * _n + k;
    }
    int Idk(double x) const
    {
        return floor(x * _n / _l);
    }

    std::vector<int> Ids(Eigen::Vector3d x, double r) const
    {
        // std::cout << "x = " << x.transpose() << " r = " << r << std::endl;
        int xs = Idk(x.x() - r);
        int ys = Idk(x.y() - r);
        int zs = Idk(x.z() - r);

        int xe = Idk(x.x() + r);
        int ye = Idk(x.y() + r);
        int ze = Idk(x.z() + r);

        std::vector<int> ids;
        for (int x = xs; x <= xe; ++x)
        {
            for (int y = ys; y <= ye; ++y)
            {
                for (int z = zs; z <= ze; ++z)
                {
                    // wrap coords around _n
                    int xmod = (x + _n) % _n;
                    int ymod = (y + _n) % _n;
                    int zmod = (z + _n) % _n;
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

    const std::vector<int>& Get(int id) const
    {
        return _boxes.at(id);
    }

    std::vector<int> Neighbors(Eigen::Vector3d x, double r, int to_del = 0) const
    {
        std::vector<int> n;

        for (int id : Ids(x, r))
        {
            const auto& box = Get(id);
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
