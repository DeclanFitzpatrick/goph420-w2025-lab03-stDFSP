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
    return f_zeta_max


def dispersion_deriv(x, b_1, b_2, p_1, p_2, H, sqrt_f):
    """
    Derived equation 1.
    """
    dF_dzeta_max = ((p_2 / p_1) * (np.sqrt(sqrt_f) / x ** 2) -
                   ((p_2/p_1) * (1/np.sqrt(H ** 2 * (b_1 ** (-2) - b_2 ** (-2))))) -
                   (2*np.pi * (1 / np.cos(2*np.pi*x)**2)))
    return dF_dzeta_max


"""def zeta_func(b_1, b_2, p_1, p_2, H):

    Using the Newton_Raphson root-finding method to solve for zeta.
    **** might not use
    
    zeta_max = H * np.sqrt(b_1 ** (-2) + b_2 ** (-2))   # equation 2, zeta in terms of love wave velocity
    sqrt_f = H * np.sqrt((b_1 ** -2) - (b_2 ** -2))

    # checking zeta_max range!
    if zeta_max < 0 or np.isnan(zeta_max):
        return None

    initial_guess = np.linspace(0.1, zeta_max, 10)
    zeta_list = []  # for the modes

    for x0 in initial_guess:
        # note that 'lambda x' is used to pass the necessary arguments to the functions (so ignores x)
        # actually solving for roots of dispersion
        zeta = root_newton_raphson(x0,
                               lambda x: dispersion(x, p_1, p_2, sqrt_f),
                               lambda x: dispersion_deriv(x, b_1, b_2, p_1, p_2, H, sqrt_f))
        if zeta is not None and zeta >= 0:
            zeta_list.append(zeta)  # only keeping positive roots
    # returning a list of valid zeta values
    return zeta_list if zeta_list else None"""


def main():
    # b_1 and B_2 are in units of m/s
    b_1 = 1900
    b_2 = 3200
    # density are in units of kg/m^3
    p_1 = 1800
    p_2 = 2500
    # thickness is in units of m
    H = 4000

    zeta_max = H * np.sqrt(b_1 ** (-2) + b_2 ** (-2))  # equation 2, zeta in terms of love wave velocity
    sqrt_f = H * np.sqrt((b_1 ** -2) - (b_2 ** -2))

    # frequency is in Hz
    freq = np.linspace(0.1, 2, 4)
    nfreq = len(freq)

    plt.figure(figsize=(10, 8))

    for j, f in enumerate(freq):
        asyms = [0.0]
        t = 0
        k = 0

        while t < zeta_max:
            t = 0.25 * (2 * k + 1) / f
            if t < zeta_max:
                asyms.append(t)
            k += 1
        asyms.append(zeta_max)
        n = len(asyms)

        plt.subplot(nfreq, 1, j + 1)  # creating subplots

        for k, tote in enumerate(asyms):  # looping over asyms list
            if k and k < n - 1:
                plt.plot([tote, tote], [-5, 5], "--r")  # plotting a red dash line at each asym
            if k < n - 1:
                zetaval = np.linspace(tote + 1e-3, asyms[k + 1] - 1e-3)  # creating zeta array between two asyms
                print(f'Zeta values: {zetaval}')
                disperse = dispersion(zetaval, p_1, p_2, H)
                print(f'Dispersion values: {disperse}')
                plt.plot(zetaval, disperse, "-b", label='Dispersion Values vs. Zeta')  # x,y
            plt.grid()
            plt.xlabel("Love Wave Velocity")
            plt.ylabel("Dispersion")
            plt.xlim((0.0, zeta_max))
            plt.ylim((-5.0, 5.0))
    plt.title('Dispersion vs Love Wave Velocity')
    plt.legend()
    plt.show()



# plotting...

    # zeta and c_L over a range of freq (modes)
        # each mode will have a curve of c_L vs f    and lambda vs. f





if __name__ == "__main__":
    main()
