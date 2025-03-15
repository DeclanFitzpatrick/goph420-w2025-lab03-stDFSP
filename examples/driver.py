import numpy as np
import matplotlib.pyplot as plt

from src.goph420_lab03.open_methods import (
    root_newton_raphson,
)


# creating function to make the plot using part 2
#def dispersion_deriv(x, b_1, b_2, p_1, p_2, H, sqrt_f):
# derived equ
    #dF_dzeta_max = ((p_2 / p_1) * (np.sqrt(sqrt_f) / x ** 2) -
                   #((p_2/p_1) * (1/np.sqrt(H ** 2 * (b_1 ** (-2) - b_2 ** (-2))))) -
                   #(2*np.pi * (1 / np.cos(2*np.pi*x)**2)))
    #return dF_dzeta_max


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

    fig, axes = plt.subplots(nrows=nfreq, sharex=True, figsize=(12, 12))

    for j, (ax, f) in enumerate(zip(axes, freq)):
        zeta_max = H * np.sqrt(b_1 ** (-2) + b_2 ** (-2))  # equation 2, zeta in terms of love wave velocity
        sqrt_f = H * np.sqrt((b_1 ** -2) - (b_2 ** -2))

        def dispersion(x):
            return (p_1 / p_2) * (np.sqrt(sqrt_f) / x) - np.tan(2 * np.pi * x)  # equation 1

        asyms = [0.0]
        a, k = 0, 0

        while a <= zeta_max:
            a = 0.25 * (2 * k + 1) / f
            if a < zeta_max:
                asyms.append(a)
            k += 1
        asyms.append(zeta_max)
        asyms = sorted(set(asyms))
        asyms = [float(a) for a in asyms]
        print(f'Asyms: {asyms}')
        n = len(asyms)

        #plt.subplot(nfreq, 1, j + 1)  # creating subplots
        for k, ak in enumerate(asyms[:-1]):
            if k > 0:
                ax.plot([ak, ak], [-5, 5], '--b', label='Asymptotes')
            if asyms[k + 1] - ak > 1e-3:
                points = np.linspace(ak + 1e-3, asyms[k + 1] - 1e-3, num=100)
                func = dispersion(points)
                ax.plot(points, func, '-r', label='Love Wave Velocity')
        ax.grid()
        ax.set_xlabel("Love Wave Velocity")
        ax.set_ylabel("Dispersion")
        ax.set_xlim((0.0, zeta_max))
        ax.set_ylim((-5.0, 5.0))
        ax.set_title(f'Dispersion vs Love Wave Velocity (f = {f} Hz)')
        ax.legend()
    plt.subplots_adjust(hspace=0.4)  # Prevents overlap
    plt.savefig("C:/Users/sydne/git/goph420/goph420-w2025-lab03-stDFSP/figures/disperse_love_wave.png")
    plt.show()



# plotting...

    # zeta and c_L over a range of freq (modes)
        # each mode will have a curve of c_L vs f    and lambda vs. f





if __name__ == "__main__":
    main()
