import zaya.particles
import zaya.sphere_visu
import numpy as np

def main():
    np.random.seed(0)
    N = 100
    r = np.random.random(N) + 0.5
    r*= 0.02
    
    visu = zaya.sphere_visu.SphereVisualizer(N)
    visu.add_box(1,1,1)

    print(zaya.particles.V(r, 0.))

    

    # return
    x = zaya.particles.cpp.rsa(r, 1000)
    
    visu.update_data(x, r)
    visu.show()


    with zaya.TTimer("FBA"):
        x, dr = zaya.particles.fba(x, r, rho=0.001, info_inverval=100)

    print(zaya.particles.V(r, dr))

    print(x, dr)
    visu.update_data(x, r)
    visu.show()
    visu.update_data(x, r+dr)
    visu.show()


    zaya.list_timings()
    

if __name__ == "__main__":
    main()
