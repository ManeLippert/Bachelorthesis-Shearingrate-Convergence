import numpy as np

def finite_first_order(y, dx, method, PERIODIC = False):
    
    #print(np.insert(np.insert(y[:,0], 0, y[:,0][-1]), len(y[:,0]) + 1, y[:,0][0]))
    
    if PERIODIC:
        
        y_periodic = []
        
        for i in range(0, y.shape[1]):
            
            y_periodic.append(np.insert(np.insert(y[:,i], 0, y[:,i][-1]), len(y[:,i]) + 1, y[:,i][0]))
        
        y = np.array(y_periodic).transpose()
    
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

def finite_second_order(y, dx, method):
    if method == 'central':
        ddy = np.array([y[i+1] - 2*y[i] + y[i-1] for i in range(1, len(y)-1)])
    elif method == 'period':
        ddy_list = []
        for i in range(0, len(y)):
            if i == 0:
                ddy_list.append(y[i+1] - 2*y[i] + y[len(y)-1])
            elif i == len(y)-1:
                ddy_list.append(y[0] - 2*y[i] + y[i-1])
            else:
                ddy_list.append(y[i+1] - 2*y[i] + y[i-1])
        
        ddy = np.array(ddy_list)
    
    return ddy/(dx*dx)