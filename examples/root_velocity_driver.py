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
        a = (0.25 * 1 / f) * (2 * k + 1)
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
            func = (H ** 2) * (b_1 ** (-2) - b_2 ** (-2))
            if func - x**2 < 0:
                return np.nan
            return (((p_2 / p_1) * np.sqrt(func - x ** 2) / x)
                    - np.tan(2 * np.pi * f * x))

        def dispersion_deriv(x):
            func = (H ** 2) * (b_1 ** (-2) - b_2 ** (-2))
            if func - x ** 2 < 0:
                return np.nan
            return ((-(p_2 / p_1) * func / (x ** 2 * np.sqrt(func - x ** 2)))
                            - 2 * np.pi * f * (1 / np.cos(2 * np.pi * f * x)) ** 2)

        guesses = []  # list for initial guesses
        for j, a in enumerate(asyms):
            if a == 0 or (a == zeta_max and dispersion(a) > 0):
                continue
            x0 = asyms[j] - 1e-3
            print(f'Initial guess: {x0}')
            if (H ** 2) * (b_1 ** (-2) - b_2 ** (-2)) - x0 ** 2 >= 0:  # Only add valid guesses
                guesses.append(x0)
            else:
                pass

        root_list = []  # storing roots
        for x0 in guesses:
            # call root finder to solve dispersion
            root_val = root_newton_raphson(x0, dispersion, dispersion_deriv)[0]
            root_list.append(root_val)

        for k, mode in enumerate(mode_list):  # store each root found in the list at each mode
            if k < len(root_list):
                mode.append(root_list[k])

    mode_list = np.array(mode_list, dtype=object)

    c_L_0 = [np.sqrt(1 / (b_1**(-2) - (r/H)**2)) for r in mode_list[0]]
    c_L_1 = [np.sqrt(1 / (b_1**(-2) - (r/H)**2)) for r in mode_list[1]]
    c_L_2 = [np.sqrt(1 / (b_1**(-2) - (r/H)**2)) for r in mode_list[2]]
    print(f'mode 0: {c_L_0}')
    print(f'mode 1: {c_L_1}')
    print(f'mode 2: {c_L_2}')

    plt.plot(freq[-len(c_L_0):], c_L_0, label="mode 0")
    plt.plot(freq[-len(c_L_1):], c_L_1, label="mode 1")
    plt.plot(freq[-len(c_L_2):], c_L_2, label="mode 2")

    plt.xlabel('Frequency (Hz)')
    plt.ylabel('c_L Values')
    plt.title('c_L vs frequency for mode 0, 1, and 2')
    plt.grid()
    plt.legend()
    plt.savefig('C:/Users/sydne/git/goph420/goph420-w2025-lab03-stDFSP/figures/root_velocities.png')
    plt.show()


if __name__ == "__main__":
    main()