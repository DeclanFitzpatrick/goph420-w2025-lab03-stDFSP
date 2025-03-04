import numpy as np
import matplotlib.pyplot as plt


# creating function to make the plot using part 2
def dispersion(b_1, b_2, p_1, p_2, c_L, H):
    Smax = H * np.sqrt(b_1 ** (-2) + c_L ** (-2))  # equation 2
    g_S = (p_1/p_2) * (np.sqrt(H**2 * (b_1 ** (-2) - b_2 ** (-2)) - Smax**2)/Smax) - np.tan(2*np.pi*Smax)  # equation 1
    return g_S

def disperson_deriv(b_1, b_2, p_1, p_2, c_L, H):


def main():
    #b_1 and B_2 are in units of m/s
    b_1 = 1900
    b_2 = 3200

    #density are in units of kg/m^3
    p_1 = 1800
    p_2 = 2500

    #thickness is in units of m
    H = 4000

    print(dispersion(b_1, b_2, p_1, p_2, c_L, H))

    #frequency is in Hz
    f = np.linspace(0, 0.1, 10)
    for k in enumerate(f):
        Stots = (1/4*f[0])*(2*k + 1)
        print(Stots)

if __name__ == "__main__":
    main()
