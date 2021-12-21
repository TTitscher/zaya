import matplotlib.pyplot as plt
import numpy as np
import zaya
from scipy.spatial import distance


class AdaptivePlot:
    def __init__(self, N):
        plt.ion()  # interactive on
        self.ax = plt.gca()
        self.ax.set_xlim((0, 1))
        self.ax.set_ylim((0, 1))
        self.ax.set_aspect("equal")
        self.cs = [plt.Circle((0, 0), 0, lw=1, fill=False) for _ in range(N)]

        for c in self.cs:
            self.ax.add_patch(c)

    def __call__(self, x, r, text=""):
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


def d_in(xs):
    distance_to_wall_0 = np.min(x)
    distance_to_wall_1 = np.min(1 - x)

    return min(np.min(distance.pdist(x)), distance_to_wall_0, distance_to_wall_1)


N = 20
L = 1  # box length
np.random.seed(6174)
x = np.random.random((N, 2))
# r = 0.1

eta = 0.9  # nominal packing density
d_out_old = 2 * L * (3 * eta / (4 * np.pi * N)) ** (1 / 3)

tau = 1.0

kappa = 0.2

for iteration in range(10):
    d = d_in(x)
    V_real = N * np.pi / (6 * L ** 3) * d ** 3
    V_virt = N * np.pi / (6 * L ** 3) * d_out_old ** 3

    dV = V_virt - V_real

    print(iteration, V_real, V_virt, dV)

    nu = np.ceil(-np.log10(dV))
    d_out_new = d_out_old * (1 - 0.5 ** nu / (N * tau))

    x_new = np.copy(x)
    for i in range(N):
        xi = x[i]

        def add(diff, sigma):
            delta = np.linalg.norm(diff)
            P = (sigma - delta) / sigma if delta < sigma else 0
            return kappa / d * P * np.asarray(diff) / delta

        # wall terms
        # x_new[i] += add((xi[0], 0), sigma=d / 2)
        # x_new[i] += add((0, xi[1]), sigma=d / 2)
        # x_new[i] += add((1 - xi[0], 0), sigma=d / 2)
        # x_new[i] += add((0, 1 - xi[1]), sigma=d / 2)

        for j in range(N):
            if i == j:
                continue

            x_new[i] += add(xi - x[j], sigma=d)

    # print(x_new - x)
    x = np.copy(x_new)

    d_out_old = d_out_new


zaya.list_timings()

# p = AdaptivePlot(N)
# p(x, r)
# p.keep()
