import matplotlib.pyplot as plt
import numpy as np
import zaya
import pytest
from rsa import rsa
from scipy.spatial import KDTree
import numba
from scipy.optimize import root_scalar, newton

from sphere_visu import SphereVisualizer
import edmd

L = 50

# @numba.jit(nopython=True)
def wrapped_euclidean(p, x):
    diff = np.abs(p-x)
    periodic_diff = np.minimum(diff, L-diff)
    return np.sum(periodic_diff**2, axis=1)**0.5

# @numba.jit(nopython=True)
def factor_in(r, d):
    """
    Calculate the factor_in in 
        d_i_in = d_i  + min_offset
    such that the spheres (centers `r` and diameters `d`) do not overlap with 
    themselves or the walls.
    """
    min_offset = np.inf
    for i in range(len(r)-1):
        current_distance = wrapped_euclidean(r[i+1:len(r)], r[i])
        radii = (d[i+1:len(r)] + d[i]) / 2
        
        offset = (current_distance - radii) 
        min_offset = min(min_offset, min(offset))

    return min_offset

def factor_out(d, eta=1.0, guess=None):
    """
    Calculate the factor_out in 
        d_i_out = factor_out * d_i 
    such that the spheres (centers `r` and diameters `d`) have a volume 
    fraction of `eta`, where the volume fraction is defined as

        eta = V_spheres / V_box = sum(pi/6 d³) / V_box
            = pi/6 factor_out³ * sum(d_i³)

    """
    V_box = L**3

    def V(offset):
        eta_num = np.pi / 6 *np.sum((d + offset)**3) / V_box 
        return eta_num - eta

    if guess is None:
        guess = 10.
    r = newton(V, x0=guess)
    return r
    # return r

    # return (V_box * eta * 6 / np.pi / np.sum(d**3))**(1/3)

# @numba.jit(nopython=True)
def get_dr_sphere(r, d, f_out, rho):
    dr_sphere = np.zeros_like(r)  # delta r resulting from an "overlap force"
    for i in range(len(r)):

        d_out_i, d_out_j = f_out + d[i], f_out + d
        sigma = (d_out_i + d_out_j) / 2

        rji = r - r[i] # unperiodic
        swap = (np.abs(rji) > L - np.abs(rji)) * np.sign(rji)
        # print(swap)
        rji -= L * swap

        abs_rji = np.sqrt(rji[:,0]**2 + rji[:, 1]**2+ rji[:, 2]**2)

        # abs_rji = np.linalg.norm(rji, axis=1)

        abs_rji[i] = np.inf

        p_ij = (sigma - abs_rji) / sigma
        # p_ij = d_out_i * d_out_j * (1 - abs_rji**2 / sigma**2 ) / sigma**2
        p_ij = np.maximum(p_ij, 0)
        assert np.all(p_ij < 1)
    
        scalar_factors = rho / d[i] * p_ij / abs_rji
            
        dr_sphere[i] -= np.sum(rji * scalar_factors.reshape(-1, 1), axis=0)

    return dr_sphere

np.random.seed(6174)

# gc = edmd.GradingCurve()
# gc.add_grading_class(2, 4, 0.12)
# gc.add_grading_class(4, 8, 0.24)
# gc.add_grading_class(8, 16, 0.32)

# box = edmd.Cube(50, 50, 50)
# phi = box.volume() * (0.22 + 0.32 + 0.12)
# radii = gc.sample(phi, 1)

N = 100
# d = 1 *
d = np.ones(N)
# r = np.random.uniform(0.2, 0.8, (3**3,3))
# d = np.random.uniform(0.2, 1.0, N)**1/3
# d = np.random.choice([1, 0.6], len(r))

np.random.seed(6174)
with zaya.TTimer("RSA"):
    r, _ = rsa(d/2, L=L)
    assert r is not None

d_in = factor_in(r, d)
print(d_in)

print(factor_in(r, d+d_in))
#
d_out = factor_out(d)

# exit()

N = len(r)

visu = SphereVisualizer(len(r))
visu.add_box(L, L, L)
# visu.update_data(r, d)
# visu.show()

# f_in = factor_in(r, d)
# assert factor_in(r, d*f_in) == pytest.approx(1)
    

def fba(r, d, rho=0.01, tau=1e4):

    f_out0 = factor_out(d, eta=0.99)
    f_out = f_out0
    f_in_old = 0
    # assert np.sum(np.pi/6 * (d*f_out)**3) == pytest.approx(0.99)

    for iteration in range(100000):

        f_in = factor_in(r, d)

        if f_in < 0:
            raise RuntimeError(f"{f_in = } !?")
        
        # visu.update_data(r, d * f_in)
        # visu
        # visu.show()
        # if iteration % 10 == 0:

        V_real = np.pi / 6 * np.sum((d+f_in) ** 3) / L**3
        V_virt = np.pi / 6 * np.sum((d+f_out)**3) / L**3

        nu = np.ceil(-np.log10(V_virt - V_real))

        print(iteration, V_real, V_virt, nu)

        # f_out -= 0.5 ** nu * f_out0 / (2 * tau)
        f_out -= 0.5 ** nu * f_out / (2 * tau)

        if f_in > f_out:
            break

        with zaya.TTimer("Force spheres"):
            dr_sphere = get_dr_sphere(r, d, f_out, rho)    

        with zaya.TTimer("Update positions"):
            # print(f"{np.max(dr_wall)  = }")
            # print(f"{np.max(dr_sphere)= }")
            # r += dr_wall
            r += dr_sphere
            r = np.mod(r, L)

            f_in_old = f_in
            # assert np.all(r >= 0)
            # assert np.all(r <= 1)

        # break
        # rho = min(rho * 1.01, 1)
    return r, f_in_old

r, f_in = fba(r, d, rho=0.1, tau=1e3)
# r, f_in = fba(r, d, tau=1e4)
# r, f_in = fba(r, d, tau=1e5)
# r, f_in = fba(r, d, tau=1e3)
print(f_in)
print(f_in)

visu.update_data(r, d * f_in)
visu.show()
zaya.list_timings()

# p = AdaptivePlot(N)
# p(r, r)
# p.keep()
