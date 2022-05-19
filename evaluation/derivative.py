import numpy as np

def finite_first_order(y, dx, method):
    if method == 'central':
        dy_central = np.array([y[i+1] - y[i-1] for i in range(1, len(y)-1)])
        return dy_central/(2*dx)
    elif method == 'forward':
        dy_forward = np.array([y[i+1] - y[i] for i in range(0, len(y)-1)])
        return  dy_forward/dx
    elif method == 'backward':
        dy_backward = np.array([y[i] - y[i-1] for i in range(1, len(y))])
        return dy_backward/dx
    else:
        raise ValueError("Method must be 'central', 'forward' or 'backward'.")

def finite_second_order(y, dx):
    ddy = np.array([y[i+1] - 2*y[i] + y[i-1] for i in range(1, len(y)-1)])
    return ddy/(dx*dx)