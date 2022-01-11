import zaya.particles
import zaya.sphere_visu
import numpy as np


def main():
    np.random.seed(0)
    N = 500
    r = np.random.uniform(2, 16, N)
    r = np.ones(N)
    r = np.sort(r)[::-1] 
    r *= 0.5 / r[0]
    # print(r)
    while(True):
        try:
            print("Running RSA", flush=True)
            x = zaya.particles.cpp.rsa(r, 1000)
            break
        except:
            r *= 0.95

    print(zaya.particles.V(r, 0))
    
    # r*= 0.02
    # return

    visu = zaya.sphere_visu.SphereVisualizer(N)
    visu.add_box(1, 1, 1)

    # return
    # print(x)

    # print(r)
    if False:
        for i in range(N):
            for j in range(i + 1, N):
                rji = x[i] - x[j]
                correction = np.round(rji)

                d = np.linalg.norm(rji - correction)
                rs = r[i] + r[j]
                if d < rs:
                    print(i, j, d, rs, d > rs)
                
        return

    visu.update_data(x, r)
    # visu.show()

    with zaya.TTimer("FBA"):
        x, dr = zaya.particles.fba(x, r, rho=0.001, info_inverval=100, iter_max=1e7, tau=1e3)

    # print(zaya.particles.V(r*FACTOR, dr))
    #
    print("dr = ", dr)
    visu.update_data(x, r+dr)
    visu.show()
    visu.update_data(x, r)
    visu.show()
    #

    zaya.list_timings()


if __name__ == "__main__":
    main()
