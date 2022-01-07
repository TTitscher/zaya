from . import _particles as cpp

import numpy as np
from itertools import cycle
from loguru import logger
from scipy.optimize import newton



def sample_grading_class(dmin, dmax, V, chunk=100):
    radii = []
    sampled_volume = 0.0

    while True:
        xs = np.random.random(chunk)
        rs = 0.5 * dmin * dmax / ((1 - xs) * dmax ** 3 + xs * dmin ** 3) ** (1 / 3)
        vs = 4.0 / 3.0 * np.pi * np.sum(rs ** 3)

        if sampled_volume + vs < V:
            # Adding the particles will not yet reach the target volume, so we
            # can add them all at once
            sampled_volume += vs
            radii += list(rs)
        else:
            # Almost reached the target volume. Adding particles one by one.
            remaining_radii = cycle(rs)
            while sampled_volume < V:
                r = next(remaining_radii)
                v = 4.0 / 3.0 * np.pi * r ** 3
                sampled_volume += v
                radii.append(r)
            break

    return np.sort(radii)[::-1]

def V(r, dr):
    return 4.0 / 3.0 * np.pi * np.sum((r + dr) ** 3)

def fba(x, r, rho=0.001, tau=1e3, iter_max=int(1e6), info_inverval=100):
    
    def error(dr):
        return V(r, dr) - 1.
    dr_out0 = newton(error, x0=100)
    print(dr_out0)

    dr_out = dr_out0

    dx = np.empty_like(x)

    for iteration in range(iter_max):
        dr_in = cpp.dx_sphere(x, r, dr_out, rho, dx)
        print(f"{dr_in = }")

        V_real, V_virt = V(r, dr_in), V(r, dr_out)
        
        info = f"{iteration:5d}: {V_real:6.3f} | {V_virt:6.3f}"

        if dr_in > dr_out:
            logger.info(info)
            break

        if iteration % info_inverval == 0:
            logger.info(info)


        # print(info)


        x += dx
        x = np.mod(x, 1.0)

        nu = np.ceil(-np.log10(V_virt - V_real))
        dr_out -= 0.5 ** nu * dr_out / (2.0 * tau)

    return x, dr_in
