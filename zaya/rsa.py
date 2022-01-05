import numpy as np
from tqdm import tqdm

def rsa(radii, L = 1., n_tries=1000, seed=6174):
    """
    Places N spheres of various `radii` in a box with length "L".
    Stops, a non-overlapping position of a sphere could not be found within
    `n_tries`.
    """
    assert np.all(2*radii < L) # must at least theoretically fit in...

    np.random.seed(seed)
    N = len(radii)
    spheres = -np.ones((N, 3))  # initial positions outside of box
    tries = np.ones(N) * np.inf

    for i, r in enumerate(tqdm(radii)):
        n = 0
        while True:
            x = np.random.uniform(r, L-r, size=3)

            d = np.sum((spheres[:i] - x) ** 2, axis=1)
            allowed_d = (radii[:i] + r)**2

            if np.all(d > allowed_d):
                # print(f"place sphere {i} at {x}")
                spheres[i] = x
                tries[i] = n
                break

            n += 1
            if n > n_tries:
                return None, tries
    return spheres, tries


if __name__ == "__main__":

    phi = 0.3
    N = 10000
    r = np.ones(N) * (phi * 3 / (4 * N * np.pi))**(1/3)

    print(4/3 * np.pi * np.sum(r**3))
    # print(r)
    spheres, tries = rsa(r, n_tries=10000)
    assert spheres is not None

    print("Average number of tries:", np.sum(tries)/len(r))
