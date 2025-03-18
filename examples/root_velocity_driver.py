import numpy as np
import matplotlib.pyplot as plt

from src.goph420_lab03.open_methods import (
    root_newton_raphson,
)


def asym_finder(f, H, b_1, b_2):
    """
    To find asymptotes.
    """
    asyms = [0.0]
    a, k = 0, 0
    zeta_max = H * np.sqrt(b_1 ** (-2) + b_2 ** (-2))  # equation 2, zeta in terms of love wave velocity
    while a <= zeta_max:
        a = 0.25 * (2 * k + 1) / f
        if a < zeta_max:
            asyms.append(a)
        k += 1
    asyms.append(zeta_max)
    return asyms


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
    freq = [0.1, 0.5, 1, 1.5, 2]

    zeta_max = H * np.sqrt(b_1 ** (-2) + b_2 ** (-2))
    mode_list = [[], [], []]

    for f in freq:
        asyms = asym_finder(f, H, b_1, b_2)  # call asyms func

        # we need to run function and derivative for each freq
        def dispersion(x):
            return (p_1 / p_2) * (np.sqrt(zeta_max**2 - x**2) / x) - np.tan(2 * np.pi * f * x)  # equation 1

        def dispersion_deriv(x):
            return ((p_2 / p_1) * (np.sqrt(zeta_max) / x ** 2) -
                            ((p_2 / p_1) * (1 / np.sqrt(H ** 2 * (b_1 ** (-2) - b_2 ** (-2))))) -
                            (2 * np.pi * (1 / np.cos(2 * np.pi * f * x) ** 2)))

        guesses = []  # list for initial guesses
        for j, a in enumerate(asyms):
            if a == 0 or (a == zeta_max and dispersion(a) > 0):
                continue
            x0 = asyms[j] - 1e-3
            guesses.append(x0)

        root_list = []  # storing roots
        for x0 in guesses:
            root_val = root_newton_raphson(x0, f, dispersion_deriv)[0]  # call root finder to solve dispersion
            root_list.append(root_val)

        for k, mode in enumerate(mode_list):  # store each root found in the list at each mode
            if k < len(root_list):
                mode.append(root_list[k])

        print(f'root list: {root_list}')


if __name__ == "__main__":
    main()