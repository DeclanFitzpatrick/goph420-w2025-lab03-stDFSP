import numpy as np
import matplotlib.pyplot as plt

from src.goph420_lab03.open_methods import (
    root_newton_raphson,
)


def main():
    # b_1 and B_2 are in units of m/s
    b_1 = 1900
    b_2 = 3200
    # density are in units of kg/m^3
    p_1 = 1800
    p_2 = 2500
    # thickness is in units of m
    H = 4000

    # frequency is in Hz
    freq = [0.1,0.5,1,1.5,2]
    nfreq = len(freq)

    def asym_finder(f, H, b_1, b_2):
        asyms = [0.0]
        a, k = 0, 0

        while a <= zeta_max:
            a = 0.25 * (2 * k + 1) / f
            if a < zeta_max:
                asyms.append(a)
            k += 1
        asyms.append(zeta_max)
        n = len(asyms)
        return asyms

    zeta_max = H * np.sqrt(b_1 ** (-2) + b_2 ** (-2))  # equation 2, zeta in terms of love wave velocity
    mode_list = [[],[],[]]

    for j, f in enumerate(freq):
        asyms = asym_finder(f, H, b_1, b_2)

        def dispersion(x):
            return (p_1 / p_2) * (np.sqrt(zeta_max**2 - x**2) / x) - np.tan(2 * np.pi * f * x)  # equation 1

        def dispersion_deriv(x):
            return ((p_2 / p_1) * (np.sqrt(zeta_max) / x ** 2) -
                            ((p_2 / p_1) * (1 / np.sqrt(H ** 2 * (b_1 ** (-2) - b_2 ** (-2))))) -
                            (2 * np.pi * (1 / np.cos(2 * np.pi * f * x) ** 2)))

        guesses = []
        for j, a in enumerate(asyms):
            if a == 0 or (a == zeta_max and dispersion(a) > 0):
                continue
            x0 = asyms[j] - 1e-3
            guesses.append(x0)

        roots = []
        for x0 in guesses:
            root_val = root_newton_raphson(x0, f, dispersion_deriv)[0]
            root_val.append(root_val)

        for k, mode in enumerate(mode_list):
            if k < len(roots):
                mode.append(roots[k])

            print(mode_list)




# plotting...

    # zeta and c_L over a range of freq (modes)
        # each mode will have a curve of c_L vs f    and lambda vs. f





if __name__ == "__main__":
    main()