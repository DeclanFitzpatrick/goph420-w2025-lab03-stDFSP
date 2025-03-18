import numpy as np
import matplotlib.pyplot as plt

from src.goph420_lab03.open_methods import (
    root_newton_raphson,
)

def dispersion_deriv(x, b_1, b_2, p_1, p_2, H, sqrt_f):
"""Derived Equation"""
    dF_dzeta_max = ((p_2 / p_1) * (np.sqrt(sqrt_f) / x ** 2) -
                    ((p_2/p_1) * (1/np.sqrt(H ** 2 * (b_1 ** (-2) - b_2 ** (-2))))) -
                    (2*np.pi * (1 / np.cos(2*np.pi*x)**2)))

# b_1 and B_2 are in units of m/s
b_1 = 1900
b_2 = 3200
# density are in units of kg/m^3
p_1 = 1800
p_2 = 2500
# thickness is in units of m
H = 4000
