import numpy as np
import matplotlib.pyplot as plt

from src.goph420_lab03.open_methods import (
    root_newton_raphson,
)


# creating function to make the plot using part 2
def dispersion(b_1, b_2, p_1, p_2, c_L, H, zeta_max):
    f_zeta_max = (p_1/p_2) * (np.sqrt(H**2 * (b_1 ** (-2) - b_2 ** (-2)) - zeta_max**2)/zeta_max) - np.tan(2*np.pi*zeta_max)  # equation 1
    return f_zeta_max


def disperson_deriv(b_1, b_2, p_1, p_2, c_L, H, zeta_max):
    dF_dzeta_max= (p_2 / p_1) * (np.sqrt(H ** 2 * (b_1 ** (-2) - b_2 ** (-2)) - zeta_max ** 2) / zeta_max ** 2) - ((p_2/p_1) * (1/np.sqrt(H ** 2 * np.sqrt(b_1 ** (-2) + b_2 ** (-2)))) - zeta_max ** 2) - (2*np.pi * (1 / np.cos(2*np.pi*zeta_max)**2))
    return dF_dzeta_max


def main():
    # b_1 and B_2 are in units of m/s
    b_1 = 1900
    b_2 = 3200

    # density are in units of kg/m^3
    p_1 = 1800
    p_2 = 2500

    # thickness is in units of m
    H = 4000

    c_L = 1

    zeta_max = H * np.sqrt(b_1 ** (-2) + c_L ** (-2))   # equation 2

    #frequency is in Hz
    freq = np.linspace(-5, 5, 5)
    asyms = [0,0]
    asyms_list = []
    for k in enumerate(freq):
        asyms = ((1/4)*freq * 2*k[0] + 1)  # creating asymptotes
        asyms_list.append(asyms)
    print(asyms_list)

    plt.plot(asyms_list, freq)  # x, y
    plt.xlabel("Zeta")
    plt.ylabel("Frequency (Hz)")
    plt.show()


if __name__ == "__main__":
    main()
