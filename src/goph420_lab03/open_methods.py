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
    eps_a = 1
    max_iter = 100
    iter_count = 1
    tol = 5e-6
    error_list = []

    while eps_a > tol and iter_count < max_iter:
        divide = -f(x) / dfdx(x)
        x += divide

        eps_a = np.abs(divide / x)
        error_list.append(eps_a)
        iter_count += 1

        if iter_count >= max_iter and eps_a > tol:
            print(f'{iter_count} iterations have been completed with a relative error of {eps_a}')

    return x, iter_count, np.array(error_list)
