#!/usr/bin/env python3


import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


if __name__ == "__main__":

    print("""This script can help with gkwdata.h5 of restarted runs,
    where the data was properly appended to, but e.g. due to a
    missing FDS.dat the time grid has jumps.
    This script tries to eliminate these jumps in the time grid.
    
    This may not treat each and every restart case correctly, in
    particular not if parameters have been changed.
    
    """)
    
    filename = 'gkwdata.h5'
    f = h5py.File(filename, "r+")
    if(f['/'].attrs['control.io_legacy'][0] == b'F'):

        #fig, ax = plt.subplots()
        plt.plot(f['/grid/time'].value.transpose(), linestyle='solid', marker='o')

        time_grid = f['/grid/time'].value.copy()
        #                         x
        #                      x
        #          x        x
        #       x        x 
        #    x        x
        # x           
        #
        
        # make a grid without negative jumps
        new_time_grid = time_grid + np.cumsum(np.roll(-1*np.append(np.minimum(np.diff(time_grid),0), [0]),1))
        
        # at every negative jump, this leaves a one-step wide plateau:
        #                    x
        #                x
        #        x   x
        #     x
        #  x
        #
        d = np.append(np.diff(time_grid),[0])
        # a mask for alle the indices which just follow on a negative jump
        mask = np.roll(d < 0, 1)
        print(mask)
        # get the
        d2 = d.copy()
        d2[np.logical_not(mask)] = 0.0
        smooth_plateaus = np.cumsum(d2)
        print(smooth_plateaus.shape)
        print(new_time_grid.shape)
        
        new_time_grid += smooth_plateaus
        # now it should be a line without plateaus.
        
        plt.plot(new_time_grid.transpose(), linestyle='solid', marker='x')
        
        plt.show()

        answer = input("Do you want to replace the old /grid/time with the new one [y/n] ? ")
        if(answer == 'y'):
            f['/grid/time_backup'] = f['/grid/time']
            del f['/grid/time']
            f['/grid/time'] = new_time_grid
            for k in f['/grid/time_backup'].attrs:
                f['/grid/time'].attrs[k] = f['/grid/time_backup'].attrs[k]
        
    else:
        print(f['/'].attrs['control.io_legacy'])
        print('sorry, this only deals with io_legacy=F data')
    f.close()
