import numpy as np
import matplotlib.pyplot as plt


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

    for j, f in enumerate(freq):
        zeta_max = H * np.sqrt(b_1 ** (-2) + b_2 ** (-2))  # equation 2, zeta in terms of love wave velocity

        asyms = [0.0]
        a, k = 0, 0

        while a <= zeta_max:
            a = 0.25 * (2 * k + 1) / f
            if a < zeta_max:
                asyms.append(a)
            k += 1
        asyms.append(zeta_max)
        n = len(asyms)

        def dispersion(x):
            return (p_1 / p_2) * (np.sqrt(zeta_max**2 - x**2) / x) - np.tan(2 * np.pi * f * x)  # equation 1

        plt.subplot(nfreq, 1, j + 1)  # creating subplots
        for k, ak in enumerate(asyms):
            if k and k < n - 1:
                plt.plot([ak, ak], [-5, 5], '--b')
            if k < n - 1:
                pnts = np.linspace(ak + 1e-3, asyms[k + 1] - 1e-3)
                func = dispersion(pnts)
                # print(f'x: {pnts}')
                # print(f'y: {func}')
                plt.plot(pnts, func, '-r')  # x y
        plt.grid()
        plt.xlabel("Love Wave Velocity")
        plt.ylabel("Dispersion")
        plt.xlim(0, zeta_max)
        plt.ylim(-5, 5)

    plt.subplots_adjust(hspace=0.4)
    plt.suptitle("Dispersion as a Function of Love Wave Velocity")
    plt.savefig("C:/Users/sydne/git/goph420/goph420-w2025-lab03-stDFSP/figures/disperse_love_wave.png")
    plt.show()


if __name__ == "__main__":
    main()
