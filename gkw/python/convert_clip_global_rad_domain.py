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
        description="Clip the radial domain of a GKW global input file and all input profiles in the associated input.prof file.")
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
    parser.add_argument('--psil',
                        help="lower boundary of the radial grid", 
                        type=float,
                        default=None
    )
    parser.add_argument('--psih',
                        help="upper boundary of the radial grid", 
                        type=float,
                        default=None
    )
    parser.add_argument('--keep-n-gridpoints',
                        help="Set True, to keep the original number of radial gridpoints. Set False, to keep the original resolution, i.e. decrease the number of gridpoints proportionally. ", 
                        type=bool,
                        default=False
    )

    parser.add_argument('--avoid-primes',
                        help="Set True, to adjust psil or psih automatically by a small amount, to avoid a prime number of gridpoints. Note that the given (not the adjusted!) values will be used to name the output folder", 
                        type=bool,
                        default=True
    )

    parser.add_argument('--output-folder',
                        help="output folder (will be created if not existent). This is where the new input.dat and input.prof will be put. You may use {n_x_grid} and {factor} to format n_x_grid or the factor into the folder name.", 
                        type=str,
                        default=None
    )
    
    args = parser.parse_args()

    if(args.output_folder is None):
        args.output_folder = ""
        if(args.psil is not None):
            args.output_folder = 'psil_{psil}'
        if(args.psih is not None):
            args.output_folder = '_'.join([args.output_folder,'psih_{psih}'])
    foldername = args.output_folder.format(psil=args.psil, psih=args.psih)

    # parse the file into a Namelist object, and take the groups member of that.
    nlist = namelist_python.read_namelist_file(args.input_file)
    n = nlist.groups
    # that object is an ordered case-insensitive dictionary
    
    if(args.psil is None and args.psih is None):
        raise ValueError("At least one of the arguments psil or psih must be given.")
    if(args.psil is None):
        args.psil = n['gridsize']['psil']
    if(args.psih is None):
        args.psih = n['gridsize']['psih']


    old_n_x_grid = n['gridsize']['nx']
    old_n_procs_x = n['gridsize']['n_procs_x']
    old_nx = old_n_x_grid / old_n_procs_x
    old_psil = n['gridsize']['psil']
    old_psih = n['gridsize']['psih']

    while True:
        # basic sanity tests
        if(n['control']['flux_tube']):
            raise ValueError("The file specified as global input file is actually a local input file.")
        if(args.psil > args.psih):
            raise ValueError("The parameter psil must be smaller than psih.")
        if(args.psil < old_psil or args.psil > old_psih):
            raise ValueError("The parameter psil must be inside the original domain: %f to %f" % (old_psil, old_psih))
        if(args.psih < old_psil or args.psih > old_psih):
            raise ValueError("The parameter psih must be inside the original domain: %f to %f" % (old_psil, old_psih))

        # now start to do the job

        complete_input_prof_data = gkwlib.read_input_prof_data(args.inputprof_file)
        if(args.keep_n_gridpoints):
            ix_l = 1
            ix_h = old_n_x_grid
            raise NotImplementedError("The switch --keep-n-gridpoints has not been implemented for the True value. Please use the separate convert_global_rad_resolution.py script.")
        else:
            # these indices are 0-based:
            ix_l = int((args.psil-old_psil)/(old_psih - old_psil) * old_n_x_grid)
            ix_h = int((args.psih-old_psil)/(old_psih - old_psil) * old_n_x_grid)
        n['gridsize']['nx'] = ix_h - ix_l

        # now we clip each and every profile
        iblock = 0

        print("Clipping...")

        grid_bounds_flag = False
        dx = (old_psih - old_psil)/old_n_x_grid
        #dx = new_block[1,0] - new_block[0,0]
        while iblock < len(complete_input_prof_data):
            if(complete_input_prof_data[iblock][0].startswith("Species")):
                # due to the not-so-great file format of
                # input.prof there is no header line specified
                # to address the block containing the
                # background profiles. Instead, it is the next block!
                iblock += 1
            # now, whether this is a block with profiles for a species, or
            # a block containing geometry profiles, we clip all columns.
            elif(complete_input_prof_data[iblock][0].startswith("Geometry")):
                pass
            else:
                iblock += 1
                continue

            ncolumns = complete_input_prof_data[iblock][1].shape[1]
            print("Block %d (%s) has %d columns." % (iblock, complete_input_prof_data[iblock][0], ncolumns))
            new_block = complete_input_prof_data[iblock][1][ix_l:ix_h,:]
            if(not grid_bounds_flag):
                n['gridsize']['psil'] = new_block[0,0] - dx/2
                n['gridsize']['psih'] = new_block[-1,0] + dx/2
                grid_bounds_flag = True
            complete_input_prof_data[iblock] = (complete_input_prof_data[iblock][0],new_block)
            iblock += 1


        # further parameters to adjust:

        # this is difficult to predict from outside GKW:
        #n['gyroaverage']['n_points_ga'] = int(n['gyroaverage']['n_points_ga'] * (n['gridsize']['nx'] * 1.0 / old_n_x_grid))

        divisors = [ d for d in range(2,int(n['gridsize']['nx']/2)) if n['gridsize']['nx'] % d == 0 ]
        if(len(divisors) == 0):
            print("WARNING: new grid size n_x_grid = %d is a prime number" % n['gridsize']['nx'])
            if(args.avoid_primes):
                print("WARNING: prime number n_x_grid is being avoided. Will now modify the range.")
                if(args.psil - dx > old_psil):
                    args.psil -= dx
                elif(args.psih + dx < old_psih):
                    args.psih += dx
                elif(args.psil + dx > old_psih):
                    args.psil += dx
                elif(args.psih - dx < old_psil):
                    args.psih -= dx
                else:
                    print("ERROR: There is no way to modify the domain to prevent a prime number. This should never happen.")
                    exit(1)
                continue

        #start with the original number of gridpoints per process..
        for i in range(0, int(n['gridsize']['nx']/2)):
            # check, if there is a nx nearby which is a divisor of the new gridsize
            new_nx = old_nx + i
            if n['gridsize']['nx'] % new_nx == 0:
                # if it is possible to take it, then take it
                break
            new_nx = old_nx - i
            # but don't take too few points per process either:
            if new_nx >= 4 and n['gridsize']['nx'] % new_nx == 0:
                break
        n['gridsize']['n_procs_x'] = int(n['gridsize']['nx'] / new_nx)

        break
        
    gkwlib.make_run_folder(foldername, input_file_namelist=nlist, input_prof_data=complete_input_prof_data)
    print("Data has been written to %s/" % foldername)
    

    
