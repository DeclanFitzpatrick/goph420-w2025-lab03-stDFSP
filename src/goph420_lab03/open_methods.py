# for step 4 part 1
# "Write Python function(s) that find the root of a given function using the Newton-Raphson Method

import numpy as np

def root_newton_raphson(x0, f, dfdx):
    """
    Find the root of a given function using the Newton-Raphson Method.

    Parameters:
    -----
    x0: initial guess.
    f: the function to find the root of.
    dfdx: derivative of the function, where dx is the step size.

    Returns:
    -----
    final_est: float.
        Final root estimate.
    iter: int.
        Iterations to convergence.
    rel_error: numpy.ndarray
        1D vector of the approx. relative error at each iteration.
    """
    x = x0
    while True:
        x = x - f(x) / dfdx(x)
        return x