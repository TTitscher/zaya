import matplotlib.pyplot as plt
import numpy as np
import zaya
import pytest
from scipy.spatial import distance, KDTree
import bisect

from sphere_visu import SphereVisualizer


def distance_to_unit_cube(xs):
    x0 = xs[:, 0]
    x1 = 1 - xs[:, 0]
    y0 = xs[:, 1]
    y1 = 1 - xs[:, 1]
    z0 = xs[:, 2]
    z1 = 1 - xs[:, 2]
    return x0, x1, y0, y1, z0, z1


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
            r = np.random.uniform(R, 1 - R, size=2)
            d = np.sum((spheres - r) ** 2, axis=1)
            if np.all(d > (2 * R) ** 2):
                # print("place sphere {i} at {r}")
                spheres[i] = r
                break

            n += 1
            if n > n_tries:
                return None
    return spheres


def estimate_d_factor_out(eta):
    return 2 * L * (eta / (np.pi * N)) ** 0.5



np.random.seed(6174)
# N = 100
# L = 1  # box length
# np.random.seed(6174)
# with zaya.TTimer("RSA"):
#     r0 = estimate_d(0.1) / 2
#     r = RSA(N, r0)
#     assert r is not None

N = 2
r = np.array([[0.3, 0.5, 0.5], [0.8, 0.5, 0.5]])
d = np.array([0.3, 0.3])


r = np.random.uniform(0.2, 0.8, (3**3,3))
d = np.ones(len(r))
# d = np.random.uniform(0.2, 0.8, 8)

N = len(r)

visu = SphereVisualizer(len(r))
visu.add_box(1, 1, 1)
# visu.update_data(r, d)
# visu.show()


def factor_in(r, d):
    """
    Calculate the factor_in in 
        d_i_in = factor_in * d_i 
    such that the spheres (centers `r` and diameters `d`) do not overlap with 
    themselves or the walls.
    """
    with zaya.TTimer("wall distance"):
        distances = distance_to_unit_cube(r)
        d_wall = min(np.min(dd / (d / 2)) for dd in distances)

    with zaya.TTimer("sphere distance"):
        d_spheres = d_wall
        for i in range(len(d)):
            for j in range(len(d)):
                if i == j:
                    continue
                current_distance = np.linalg.norm(r[i] - r[j])
                radii= (d[i] + d[j])/2
                factor = current_distance / radii
                d_spheres = min(d_spheres, factor)

    return min(d_wall, d_spheres)

def factor_out(r, d, eta=1.0):
    """
    Calculate the factor_out in 
        d_i_out = factor_out * d_i 
    such that the spheres (centers `r` and diameters `d`) have a volume 
    fraction of `eta`, where the volume fraction is defined as

        eta = V_spheres / V_box = sum(pi/6 d³) / V_box
            = pi/6 factor_out³ * sum(d_i³)

    """
    V_box = 1.
    return (V_box * eta * 6 / np.pi / np.sum(d**3))**(1/3)

f_in = factor_in(r, d)
assert factor_in(r, d*f_in) == pytest.approx(1)
    
f_out0 = factor_out(r, d, eta=0.99)
f_out = f_out0
assert np.sum(np.pi/6 * (d*f_out)**3) == pytest.approx(0.99)

tau = 2000

rho = 0.1


for iteration in range(10000):
    f_in = factor_in(r, d) 

    if f_in < 0:
        raise RuntimeError(f"{f_in = } !?")
    
    visu.update_data(r, d * f_in)
    # if iteration % 100 == 0:
    # visu.show()

    V_real = np.pi / 6 * np.sum((d*f_in) ** 3)
    V_virt = np.pi / 6 * np.sum((d*f_out)**3)

    nu = np.ceil(-np.log10(V_virt - V_real))

    print(iteration, V_real, V_virt, nu)

    f_out -= 0.5 ** nu * f_out0 / (2 * tau)

    if f_in > f_out:
        break

    with zaya.TTimer("Force spheres"):
        dr_sphere = np.zeros_like(r)  # delta r resulting from an "overlap force"

        for i in range(N):
            for j in range(N):
                if i == j:
                    continue
                d_out_i, d_out_j = f_out * d[i], f_out * d[j]

                rji = r[j] - r[i]
                abs_rji = np.linalg.norm(rji)

                overlap = (d_out_i + d_out_j)/2 - abs_rji
                if overlap < 0:
                    # accounts for 1_ij in eq.(2)
                    continue

                def V_sphere_intersection(R, r, d):
                    return np.pi * (R+r-d)**2 * (d**2 + 2*d*r -3 *r**2 + 2*d*R + 6*r*R-3*R**2) / (12 * d)

                p_ij = d_out_i * d_out_j * (abs_rji**2 / (0.25 * (d_out_i + d_out_j)**2) - 1)
                # sigma = 0.5 * (d_out_i + d_out_j)
                # p_ij = -(sigma - abs_rji) / sigma
                # p_ij = -V_sphere_intersection(d_out_i, d_out_j, abs_rji)
                # p_ij = abs_rji**2 - d_out_i * d_out_j
#
                F_ij = rho * p_ij * rji / abs_rji
                dr_sphere[i] += F_ij / d[i]

    with zaya.TTimer("Force wall"):
        dr_wall = np.zeros_like(r)  
        distance_wall = distance_to_unit_cube(r)

        current_radius = d * f_out / 2
        allowed_distance_wall = current_radius
        
        overlaps = [np.maximum(allowed_distance_wall- distance, 0) for distance in distance_wall]

        # overlap_V = [1/3 * np.pi * overlap**2 * (3 * current_radius - overlap) for overlap in overlaps]

        scaled_overlaps = [V/current_radius*2 / d for V in overlaps]
       
        dr_wall[:, 0] += rho * scaled_overlaps[0] 
        dr_wall[:, 0] -= rho * scaled_overlaps[1]

        dr_wall[:, 1] += rho * scaled_overlaps[2]
        dr_wall[:, 1] -= rho * scaled_overlaps[3]

        dr_wall[:, 2] += rho * scaled_overlaps[4]
        dr_wall[:, 2] -= rho * scaled_overlaps[5]

    # print(dr_wall)
    # print(dr_sphere)


    with zaya.TTimer("Update positions"):
        r += dr_wall
        r += dr_sphere

    # kappa = min(kappa * 1.01, 0.1)
    
visu.update_data(r, d * f_in)
visu.show()
zaya.list_timings()

# p = AdaptivePlot(N)
# p(r, r)
# p.keep()
