import matplotlib.pyplot as plt
import numpy as np
import zaya
import pytest
from scipy.spatial import distance, KDTree
import bisect


class AdaptivePlot:
    def __init__(self, N):
        plt.ion()  # interactive on
        self.ax = plt.gca()
        self.ax.set_xlim((0, 1))
        self.ax.set_ylim((0, 1))
        self.ax.set_aspect("equal")
        self.cs = [plt.Circle((0, 0), 0, lw=1, fill=False) for _ in range(N)]
        self.Fs = [plt.Arrow((0, 0), 0, lw=1, fill=False) for _ in range(N)]

        for c in self.cs:
            self.ax.add_patch(c)

    def __call__(self, x, r, F, text=""):
        for c, x in zip(self.cs, x):
            c.set_center(x)
            c.set_radius(r)
        self.ax.set_title(text)
        # self.ax.relim()
        self.ax.autoscale_view(True, True, True)

        self.ax.figure.canvas.flush_events()

    def close(self):
        plt.ioff()  # interactive off
        plt.close()

    def keep(self):
        """
        Deactivates the interactive mode to show the plot.
        """
        plt.ioff()  # interactive off
        plt.show()


def show_circles_and_force(xs, d_in, d_out=None, Fs_wall=None, Fs_sphere=None):
    if Fs_wall is None:
        Fs_wall = np.zeros_like(xs)
    if Fs_sphere is None:
        Fs_sphere = np.zeros_like(xs)
    assert len(xs) == len(Fs_wall) == len(Fs_sphere)
    ax = plt.gca()
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 1))
    ax.set_aspect("equal")

    for x, Fwall, Fsphere in zip(xs, Fs_wall, Fs_sphere):
        ax.add_patch(plt.Circle(x, d_in / 2, fill=True, color="gray"))
        if d_out is not None:
            ax.add_patch(plt.Circle(x, d_out / 2, fill=False, lw=2, color="gray"))
        ax.arrow(
            *x, *(Fwall + Fsphere), width=0.01, color="red", length_includes_head=True
        )
        ax.arrow(*x, *Fwall, width=0.005, color="green", length_includes_head=True)
        ax.arrow(*x, *Fsphere, width=0.005, color="blue", length_includes_head=True)

    plt.show()


def distance_to_unit_cube(xs):
    x0 = xs[:, 0]
    x1 = 1 - xs[:, 0]
    y0 = xs[:, 1]
    y1 = 1 - xs[:, 1]
    return x0, x1, y0, y1


def d_in(xs, tree):
    with zaya.TTimer("wall distance"):
        d = distance_to_unit_cube(xs)
        d_wall = min(np.min(x) for x in d) * 2

    with zaya.TTimer("sphere distance"):
        d, i = tree.query(x, k=2)
        d_spheres = np.min(d[:, 1])

    return min(d_wall, d_spheres)


def test():
    x = np.array([[0.1, 0.2], [0.3, 0.4]])
    x0, x1, y0, y1 = distance_to_unit_cube(x)

    assert x0[0] == pytest.approx(0.1)
    assert x0[1] == pytest.approx(0.3)
    assert x1[0] == pytest.approx(0.9)
    assert x1[1] == pytest.approx(0.7)

    assert y0[0] == pytest.approx(0.2)
    assert y0[1] == pytest.approx(0.4)
    assert y1[0] == pytest.approx(0.8)
    assert y1[1] == pytest.approx(0.6)

import numba

def RSA(N, R, n_tries=1000):
    """
    Places N spheres of radius R in a unit square. Stops, if n_tries are
    exceeded.
    """
    spheres = np.zeros((N, 2))
    for i in range(N):
        n = 0
        while True:
            x = np.random.uniform(R, 1 - R, size=2)
            d = np.sum((spheres - x) ** 2, axis=1)
            if np.all(d > (2 * R) ** 2):
                # print("place sphere {i} at {x}")
                spheres[i] = x
                break

            n += 1
            if n > n_tries:
                return None
    return spheres


def estimate_d(eta):
    return 2 * L * (eta / (np.pi * N)) ** 0.5


test()


# x = np.array([[0.3, 0.5], [0.8, 0.2]])
# r = 0.1



    # print(all_neighbors)
    # for i, p in enumerate(x):
        # q = tree.query_ball_point(p, r=2*d_out)
        # print(i in q)

def sphere_force2(x, sigma, kappa, all_neighbors):
    F_sphere = np.zeros_like(x)

    return F_sphere

@numba.njit(fastmath=True)
def sphere_force(x, sigma, kappa, values, idx):
    F_sphere = np.zeros_like(x)

    for i_sphere in range(N):
        indices = values[idx[i_sphere]:idx[i_sphere+1]]

        diffs = x[i_sphere] - x[indices]

        deltas = (diffs[:,0]**2 + diffs[:,1]**2)**0.5
        # deltas = np.linalg.norm(diffs, axis=1)

        Ps = np.maximum((sigma**2- deltas**2) / sigma**2, 0)
        # assert np.all(Ps < 1)
        factor = kappa / sigma* Ps / deltas
        factor = factor.reshape(-1, 1)

        F = np.sum(diffs * factor, axis=0)
        F_sphere[i_sphere] = F
    return F_sphere


@profile
def fba(x):
    d_out0 = estimate_d(1.0)
    d_out = d_out0

    # show_circles_and_force(x, r0, d_out)

    tau = 2000

    kappa = 0.05 / N


    for iteration in range(100000):
        with zaya.TTimer("KDTree| build"):
            if iteration < 10 or iteration % 100 == 0:
                tree = KDTree(x)

    

        d = d_in(x, tree)
        if d < 0:
            show_circles_and_force(x, d, d_out, F_wall, F_sphere)
            raise RuntimeError(f"{d = } !?")
        # V_real = N * np.pi / (6 * L ** 3) * d ** 3
        # V_virt = N * np.pi / (6 * L ** 3) * d_out_old ** 3
        V_real = N * np.pi / 4 * d ** 2 / L ** 2
        V_virt = N * np.pi / 4 * d_out ** 2 / L ** 2

        dV = V_virt - V_real
        if dV < 0:
            break

        print(iteration, V_real, V_virt, dV, flush=True)

        nu = np.ceil(-np.log10(dV))
        d_out -= 0.5 ** nu * d_out0 / (2 * tau)


        with zaya.TTimer("Force wall"):
            F_wall = np.zeros_like(x)  # something like an "overlap force"
            distance_wall = distance_to_unit_cube(x)

            sigma_wall = d_out / 2.0

            Pwall = [
                np.maximum((sigma_wall**2 - distance**2) / sigma_wall**2, 0)
                for distance in distance_wall
            ]
            # print(Pwall[0])
            # print(Pwall[1])

            for pp in Pwall:
                assert np.all(pp < 1)

            F_wall[:, 0] += kappa / sigma_wall * Pwall[0]
            F_wall[:, 0] -= kappa / sigma_wall * Pwall[1]

            F_wall[:, 1] += kappa / sigma_wall * Pwall[2]
            F_wall[:, 1] -= kappa / sigma_wall * Pwall[3]

        if iteration < 10 or iteration % 100 == 0:

            with zaya.TTimer("Force spheres"):
                sigma_sphere = d_out
                with zaya.TTimer("KDTree| querry all"):
                    all_neighbors = tree.query_ball_point(x, r=3*sigma_sphere)

                with zaya.TTimer("KDTree| remove"):
                    idx = [0]
                    for i, neighbors in enumerate(all_neighbors):
                        del neighbors[bisect.bisect_left(neighbors, i)]
                        idx.append(idx[i] + len(neighbors))
                    values = np.concatenate(all_neighbors, dtype=int, casting="unsafe")
                    idx = np.array(idx)


        F_sphere = sphere_force(x, sigma_sphere, kappa, values, idx)
                    
        with zaya.TTimer("Update positions"):
            x += F_wall
            x += F_sphere

        # kappa = min(kappa*1.01, 0.1)
    
    return x, d, d_out

np.random.seed(6175)
N = 500
L = 1  # box length
# np.random.seed(6174)
with zaya.TTimer("RSA"):
    r0 = estimate_d(0.1) / 2
    x = RSA(N, r0)
    assert x is not None

x, d, d_out = fba(x)
show_circles_and_force(x, d)

zaya.list_timings()

# p = AdaptivePlot(N)
# p(x, r)
# p.keep()
