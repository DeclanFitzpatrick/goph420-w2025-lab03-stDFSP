import numpy as np
import matplotlib.pyplot as plt

from src.goph420_lab03.open_methods import (
    root_newton_raphson,
)


# creating function to make the plot using part 2
def dispersion(x, b_1, b_2, p_1, p_2, H, zeta_max, sqrt_f):
    """
    Uses equation 1.
    """
    f_zeta_max = (p_1 / p_2) * (np.sqrt(sqrt_f)/zeta_max) - np.tan(2*np.pi*zeta_max)  # equation 1
    return f_zeta_max


def dispersion_deriv(x, b_1, b_2, p_1, p_2, H, zeta_max, sqrt_f):
    """
    Derived equation 1.
    """
    dF_dzeta_max= (p_2 / p_1) * (np.sqrt(sqrt_f) / zeta_max ** 2) - ((p_2/p_1) * (1/np.sqrt(H ** 2 * np.sqrt(b_1 ** (-2) + b_2 ** (-2)))) - zeta_max ** 2) - (2*np.pi * (1 / np.cos(2*np.pi*zeta_max)**2))
    return dF_dzeta_max

def zeta_func(f, b_1, b_2, p_1, p_2, H):
    """
    Using the Newton_Raphson root-finding method to solve for zeta.
    """
    x0 = 0.1  # initial guess
    zeta_max = H * np.sqrt(b_1 ** (-2) + b_2 ** (-2))   # equation 2, zeta in terms of love wave velocity
    sqrt_f = np.maximum(H**2 * (b_1 ** (-2) - b_2 ** (-2)) - zeta_max**2, 0)  # np.maximum is to ensure no neg values

    return root_newton_raphson(x0,
                               lambda x: dispersion(x, b_1, b_2, p_1, p_2, H, zeta_max, sqrt_f),
                               lambda x: dispersion_deriv(x, b_1, b_2, p_1, p_2, H, zeta_max, sqrt_f))
    # 'lambda x' is used to pass the necessary arguments to the functions (so ignores x)


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
    freq = np.linspace(-5, 5, 5)
    c_L_list = []
    lambda_L_list = []

    for f in freq:
        if f != 0:  # to avoid dividing by zero
            zeta = zeta_func(f, b_1, b_2, p_1, p_2, H)
            if zeta is not None:
                # solving for various c_L depending on zetas
                c_L = (1 / (b_1**-2 - (zeta / H)**2))**0.5  # equation 2 rearranged
                lambda_L = c_L/f  # equation 3
                c_L_list.append(c_L)
                lambda_L_list.append(lambda_L)
            else:
                c_L_list.append(np.nan)
                lambda_L_list.append(np.nan)
        else:
            c_L_list.append(np.nan)
            lambda_L_list.append(np.nan)

    # to plot asymptotes
    asyms_list = []
    for k, f in enumerate(freq):  # enumerate requires index, value
        asyms = (1/4) * f * 2 * k + 1 # creating asymptotes
        asyms_list.append(asyms)
    print(asyms_list)

    # plotting...

    # plot 1
    plt.figure()
    plt.plot(freq, c_L_list, label="Love Wave Velocity")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Love Wave Velocity (m/s)")
    plt.legend()
    plt.show()



if __name__ == "__main__":
    main()
