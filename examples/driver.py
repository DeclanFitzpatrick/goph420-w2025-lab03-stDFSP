import numpy as np
import matplotlib.pyplot as plt

from src.goph420_lab03.open_methods import (
    root_newton_raphson,
)


# creating function to make the plot using part 2
def dispersion(x, p_1, p_2, sqrt_f):
    """
    Uses equation 1.
    """
    f_zeta_max = (p_1 / p_2) * (np.sqrt(sqrt_f)/x) - np.tan(2*np.pi*x)  # equation 1
    print(f'Dispersion at x: {x}, where f(x) = {f_zeta_max} ')
    return f_zeta_max


def dispersion_deriv(x, b_1, b_2, p_1, p_2, H, sqrt_f):
    """
    Derived equation 1.
    """
    dF_dzeta_max= (p_2 / p_1) * (np.sqrt(sqrt_f) / x ** 2) - ((p_2/p_1) * (1/np.sqrt(H ** 2 * np.sqrt(
        b_1 ** (-2) + b_2 ** (-2)))) - x ** 2) - (2*np.pi * (1 / np.cos(2*np.pi*x)**2))
    print(f'Dispersion at x: {x}, where dF/dx(x) = {dF_dzeta_max} ')
    return dF_dzeta_max


def zeta_func(f, b_1, b_2, p_1, p_2, H):
    """
    Using the Newton_Raphson root-finding method to solve for zeta.
    """
    initial_guess = [0.1, 0.5, 1.0, 1.5]
    zeta_max = H * np.sqrt(b_1 ** (-2) + b_2 ** (-2))   # equation 2, zeta in terms of love wave velocity
    sqrt_f = np.maximum(H**2 * (b_1 ** (-2) - b_2 ** (-2)) - f**2, 0)  # np.maximum is to ensure no neg values

    # checking zeta_max range!
    if zeta_max < 0 or np.isnan(zeta_max):
        return None

    print(f'sqrt_f: {sqrt_f}, zeta_max: {zeta_max}')

    for x0 in initial_guess:
        # note that 'lambda x' is used to pass the necessary arguments to the functions (so ignores x)
        zeta = root_newton_raphson(x0,
                               lambda x: dispersion(x, b_1, b_2, p_1),
                               lambda x: dispersion_deriv(x, b_1, b_2, p_1, p_2, H, sqrt_f))
        print(f'Initial x0: {x0}, Zeta: {zeta}')
        if zeta is not None and zeta >= 0:
            return zeta

    return None


def main():
    # b_1 and B_2 are in units of m/s
    b_1 = 1900
    b_2 = 3200

    # density are in units of kg/m^3
    p_1 = 1800
    p_2 = 2500

    # thickness is in units of m
    H = 4000

    #frequency is in Hz
    freq = np.linspace(0.1, 5, 10)
    c_L_list = []
    lambda_L_list = []

    # to plot asymptotes
    asyms_list = []

    for f in freq:
        if f != 0:
            asyms = (1 / (4 * f)) + 1 + (0.1 / f**2)  # creating asymptotes
            asyms_list.append(asyms)
            zeta = zeta_func(f, b_1, b_2, p_1, p_2, H)
            print(f'Zeta: {zeta}')

            if zeta is not None and zeta >=0:
                # solving for various c_L depending on zetas
                c_L_denom = (b_1**-2 - (zeta / H)**2)
                if c_L_denom > 0:
                    c_L = (1 / c_L_denom)**0.5  # equation 2 rearranged
                    lambda_L = c_L/f  # equation 3
                    c_L_list.append(c_L)
                    lambda_L_list.append(lambda_L)
                else:
                    c_L_list.append(np.nan)
                    lambda_L_list.append(np.nan)
            else:
                c_L_list.append(np.nan)
                lambda_L_list.append(np.nan)
        else:
            c_L_list.append(np.nan)
            lambda_L_list.append(np.nan)
            asyms_list.append(np.nan)

        # printing lists to debug UGH
        print(f'c_L values: {c_L_list}')
        print(f'wavelength values: {lambda_L_list}')
        print(f'asyms: {asyms_list}')

    # plotting...

    # zeta and c_L over a range of freq (modes)
        # each mode will have a curve of c_L vs f    and lambda vs. f

    for asyms in asyms_list:
        if not np.isnan(asyms):
            plt.axvline(x=asyms, color='r', linestyle='--')


if __name__ == "__main__":
    main()
