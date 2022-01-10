import zaya.particles
import zaya.sphere_visu
import numpy as np


def main():
    np.random.seed(0)
    N = 4
    r = np.random.choice([0.1, 0.3], N)
    FACTOR = 1.0
    print(zaya.particles.V(r * FACTOR, 0.0))
    # r*= 0.02
    # return

    visu = zaya.sphere_visu.SphereVisualizer(N)
    visu.add_box(1, 1, 1)

    # return
    x = zaya.particles.cpp.rsa(r * FACTOR + 0.01, 1000)

    for i in range(N):
        for j in range(i + 1, N):
            d = np.linalg.norm(x[i] - x[j])
            rs = r[i] + r[j]
            print(d, rs, d > rs)

    return

    visu.update_data(x, r * FACTOR)
    visu.show()

    with zaya.TTimer("FBA"):
        x, dr = zaya.particles.fba(x, r, rho=0.001, info_inverval=1, iter_max=100)

    # print(zaya.particles.V(r*FACTOR, dr))
    #
    # print(x, dr)
    # visu.update_data(x, r*FACTOR)
    # visu.show()
    # visu.update_data(x, r*FACTOR+dr)
    # visu.show()
    #

    zaya.list_timings()


if __name__ == "__main__":
    main()
