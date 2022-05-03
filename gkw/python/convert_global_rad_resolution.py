#!/usr/bin/env python3
import numpy as np
import sys
import h5py
import argparse
import shutil
import os
import scipy.interpolate

try:
    import namelist_python
except:
    # if the GKW_HOME/python folder is not contained in the PYTHONPATH
    import os
    sys.path.append(os.path.join(os.environ['GKW_HOME'],'python'))
    import namelist_python

import gkwlib

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Change the radial resolution of a GKW global input file, and its profiles in the associated input.prof file.")
    parser.add_argument('--input-file',
                        help="GKW local input file (input file)", 
                        type=str,
                        default='input.dat'
    )
    parser.add_argument('--inputprof-file',
                        help="GKW input.prof file (input file)", 
                        type=str,
                        default='input.prof'
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--n_x_grid',
                        help="new number of radial grid points", 
                        type=int,
                        default=None
    )
    group.add_argument('--factor',
                        help="factor, to augment number of radial grid points", 
                        type=float,
                        default=None
    )

    parser.add_argument('--output-folder',
                        help="output folder (will be created if not existent). This is where the new input.dat and input.prof will be put. You may use {n_x_grid} and {factor} to format n_x_grid or the factor into the folder name.", 
                        type=str,
                        default=None
    )
    
    args = parser.parse_args()

    if(args.n_x_grid is not None):
        if(args.output_folder is None):
            args.output_folder = 'n_x_grid_{n_x_grid}'
        foldername = args.output_folder.format(n_x_grid=args.n_x_grid, factor=args.factor)
    else:
        if(args.output_folder is None):
            args.output_folder = 'n_x_grid_factor_{factor}'
        foldername = args.output_folder.format(n_x_grid=args.n_x_grid, factor=args.factor)

    # parse the file into a Namelist object, and take the groups member of that.
    nlist = namelist_python.read_namelist_file(args.input_file)
    n = nlist.groups
    # that object is an ordered case-insensitive dictionary

    complete_input_prof_data = gkwlib.read_input_prof_data(args.inputprof_file)

    old_n_x_grid = n['gridsize']['nx']
    if(args.n_x_grid is None):
        n['gridsize']['nx'] = int(np.round(n['gridsize']['nx'] * args.factor))
    else:
        args.factor = args.n_x_grid / old_n_x_grid
        n['gridsize']['nx'] = args.n_x_grid

    # now we try to interpolate each and every profile
    iblock = 0

    print("Interpolation from n_x_grid %d to %d." % (old_n_x_grid, n['gridsize']['nx']))
    
    while iblock < len(complete_input_prof_data):
        if(complete_input_prof_data[iblock][0].startswith("Species")):
            # due to the not-so-great file format of
            # input.prof there is no header line specified
            # to address the block containing the
            # background profiles. Instead, it is the next block!
            iblock += 1
        # now, whether this is a block with profiles for a species, or
        # a block containing geometry profiles, we intepolate all columns.
        elif(complete_input_prof_data[iblock][0].startswith("Geometry")):
            pass
        else:
            iblock += 1
            continue

        ncolumns = complete_input_prof_data[iblock][1].shape[1]
        print("Block %d (%s) has %d columns." % (iblock, complete_input_prof_data[iblock][0], ncolumns))
        new_block = np.zeros((n['gridsize']['nx'], ncolumns))

        xgrid = complete_input_prof_data[iblock][1][:,0]
        new_xgrid = n['gridsize']['psil'] + (n['gridsize']['psih']-n['gridsize']['psil'])*(np.arange(1,n['gridsize']['nx']+1)-0.5)/n['gridsize']['nx']
        new_block[:,0] = new_xgrid

        for icol in range(1,ncolumns):
            smoothing_cond = 1.0e0
            tck = scipy.interpolate.splrep(xgrid,complete_input_prof_data[iblock][1][:,icol], k=3, s=smoothing_cond)
            # evaluate the profile
            new_block[:,icol] = scipy.interpolate.splev(new_xgrid, tck, der=0)

        complete_input_prof_data[iblock] = (complete_input_prof_data[iblock][0],new_block)
        iblock += 1


    # further parametrs to adjust:
    n['gyroaverage']['n_points_ga'] = n['gyroaverage']['n_points_ga'] * args.factor
        
    gkwlib.make_run_folder(foldername, input_file_namelist=nlist, input_prof_data=complete_input_prof_data)
    print("Data has been written to %s/" % foldername)
    

    
