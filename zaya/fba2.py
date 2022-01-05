import matplotlib.pyplot as plt
import numpy as np
import zaya
import pytest
from rsa import rsa
from scipy.spatial import KDTree
import numba

from sphere_visu import SphereVisualizer



@numba.jit(nopython=True)
def wrapped_euclidean(p, x):
    diff = np.abs(p-x)
    periodic_diff = np.minimum(diff, 1-diff)
    return np.sum(periodic_diff**2, axis=1)**0.5

@numba.jit(nopython=True)
def factor_in(r, d):
    """
    Calculate the factor_in in 
        d_i_in = factor_in * d_i 
    such that the spheres (centers `r` and diameters `d`) do not overlap with 
    themselves or the walls.
    """
    d_spheres = np.inf
    for i in range(len(r)-1):
        current_distance = wrapped_euclidean(r[i+1:len(r)], r[i])
        radii = (d[i+1:len(r)] + d[i]) / 2
        
        factor = current_distance / radii
        d_spheres = min(d_spheres, min(factor))

    return d_spheres

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

@numba.jit(nopython=True)
def get_dr_sphere(r, d, f_out, rho):
    dr_sphere = np.zeros_like(r)  # delta r resulting from an "overlap force"
    for i in range(len(r)):

        d_out_i, d_out_j = f_out * d[i], f_out * d
        sigma = (d_out_i + d_out_j) / 2

        rji = r - r[i] # unperiodic
        swap = (np.abs(rji) > 1 - np.abs(rji)) * np.sign(rji)
        # print(swap)
        rji -= swap

        abs_rji = np.sqrt(rji[:,0]**2 + rji[:, 1]**2+ rji[:, 2]**2)

        # abs_rji = np.linalg.norm(rji, axis=1)

        abs_rji[i] = np.inf

        # p_ij = (sigma - abs_rji) / sigma
        p_ij = d_out_i * d_out_j * (1 - abs_rji**2 / sigma**2 ) / sigma**2
        p_ij = np.maximum(p_ij, 0)
        assert np.all(p_ij < 1)
    
        scalar_factors = rho / d[i] * p_ij / abs_rji
            
        dr_sphere[i] -= np.sum(rji * scalar_factors.reshape(-1, 1), axis=0)

    return dr_sphere



np.random.seed(6174)
N = 500
d = np.ones(N)
# r = np.random.uniform(0.2, 0.8, (3**3,3))
# d = np.random.uniform(0.2, 1.0, len(r))
# d = np.random.choice([1, 0.6], len(r))

np.random.seed(6174)
with zaya.TTimer("RSA"):
    r, _ = rsa(d/2 * 0.01)
    assert r is not None


N = len(r)

visu = SphereVisualizer(len(r))
visu.add_box(1, 1, 1)
# visu.update_data(r, d)
# visu.show()

# f_in = factor_in(r, d)
# assert factor_in(r, d*f_in) == pytest.approx(1)
    

@profile
def fba(r, d, rho=0.01, tau=1e4):

    f_out0 = factor_out(r, d, eta=0.99)
    f_out = f_out0
    assert np.sum(np.pi/6 * (d*f_out)**3) == pytest.approx(0.99)

    for iteration in range(1000):

        f_in = factor_in(r, d)

        if f_in < 0:
            raise RuntimeError(f"{f_in = } !?")
        
        # visu.update_data(r, d * f_in)
        # if iteration % 100 == 0:
            # visu.show()

        V_real = np.pi / 6 * np.sum((d*f_in) ** 3)
        V_virt = np.pi / 6 * np.sum((d*f_out)**3)

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
            r = np.mod(r, 1)
            # assert np.all(r >= 0)
            # assert np.all(r <= 1)

        # break
        # rho = min(rho * 1.01, 1)
    return r, f_in

r, f_in = fba(r, d, tau=1e3)
print(f_in)

visu.update_data(r, d * f_in)
visu.show()
zaya.list_timings()

# p = AdaptivePlot(N)
# p(r, r)
# p.keep()
