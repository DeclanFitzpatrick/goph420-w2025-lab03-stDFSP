import numpy as np
import matplotlib.pyplot as plt


#creating function to make the plot
def dispersion(b_1, b_2, p_1, p_2, H):
    Smax = H**2 * np.sqrt(b_1 **-2 + b_2 **-2)
    g_S = (((p_1)/(p_2)) * ((Smax - Smax**2)/(Smax))) - np.tan(2*np.pi*Smax)
    return g_S

def main():
    #b_1 and B_2 are in units of m/s
    b_1 = 1900
    b_2 = 3200

    #density are in units of kg/m^3
    p_1 = 1800
    p_2 = 2500

    #thickness is in units of m
    H = 4000


    print(dispersion(b_1, b_2, p_1, p_2, H))

    #frequency is in Hz
    f = np.linspace(0, 0.1, 10)
    for k in enumerate(f):
        Stots = (1/4*f[0])*(2*k + 1)
        print(Stots)

if __name__ == "__main__":
    main()
