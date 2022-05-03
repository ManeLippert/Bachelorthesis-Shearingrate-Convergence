#!/usr/bin/env python3

import numpy as np
import sys
import h5py
import argparse
import shutil
import os

try:
    import namelist_python
except:
    # if the GKW_HOME/python folder is not contained in the PYTHONPATH
    import os
    sys.path.append(os.path.join(os.environ['GKW_HOME'],'python'))
    import namelist_python

import gkwlib

TINY_NUMBER = 1.0e-10
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a GKW global input file to a local input file.")
    parser.add_argument('--local-input-file',
                        help="GKW local input file (output file)", 
                        type=str,
                        default='-'
    )
    parser.add_argument('--global-input-file',
                        help="GKW local input file (input file)", 
                        type=str,
                        default='input.dat'
    )
    parser.add_argument('--gkwdatah5-file',
                        help="a GKW HDF5 data file (input file)", 
                        type=str,
                        default='gkwdata.h5'
    )
    parser.add_argument('--geomdat-file',
                        help="GKW geom.dat file (input file)", 
                        type=str,
                        default='geom.dat'
    )
    parser.add_argument('--inputprof-file',
                        help="GKW input.prof file (input file)", 
                        type=str,
                        default='input.prof'
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--radial-ix',
                        help="radial index (0-based)", 
                        type=int,
                        default=None
    )
    group.add_argument('--radial-x',
                        help="radial x position", 
                        type=float,
                        default=None
    )
    args = parser.parse_args()

    if(args.local_input_file == '-'):
        output = sys.stdout
    else:
        output = open(args.local_input_file, 'w')

    # parse the file into a Namelist object, and take the groups member of that.
    nlist = namelist_python.read_namelist_file(args.global_input_file)
    n = nlist.groups
    # that object is an ordered case-insensitive dictionary

    # sanity check
    if(n['control']['flux_tube']):
        raise ValueError("The file specified as global input file is actually a local input file.")

    n['control']['flux_tube'] = True

    psil = n['gridsize']['psil']
    psih = n['gridsize']['psih']
    if(args.radial_ix is not None):
        if (args.radial_ix < 0 or args.radial_ix >= n['gridsize']['nx']):
            raise ValueError("The value radial-ix %d is invalid. (n_x_grid = %d)" % (args.radial_ix, n['control']['nx']))
        else:
            args.radial_x = psil + (psih-psil)*(args.radial_ix-0.5)/n['gridsize']['nx']
    else:
        if (args.radial_x < psil or args.radial_x > psih):
            raise ValueError("The value radial-x %f is invalid. (x grid is on interval %f to %f)" % (args.radial_x, psil, psih))
        else:
            args.radial_ix = int(np.round((args.radial_x - psil) * n['gridsize']['nx']/(psih-psil)+0.5))

    print("! This local run was generated from radial position", file=output)
    print("!  ix = %d" % args.radial_ix, file=output)
    print("!  x = %f" % args.radial_x, file=output)
    print("! of", file=output)
    print("!  global_parent = %s" % os.path.abspath(args.global_input_file), file=output)
    print("! ", file=output)

    # pick background parameters at this radial position

    for isp in range(n['gridsize']['number_of_species']):
        if(n['species'][isp]['dens_prof_type'] == 'file'):
            # I can read input.prof without running GKW
            complete_input_prof_data = gkwlib.read_input_prof_data(args.inputprof_file)
            for iblock in range(len(complete_input_prof_data)):
                if(complete_input_prof_data[iblock][0].startswith("Species")):
                    # check if this is indeed the input.prof entry for this species:
                    if(abs(n['species'][isp]['mass'] - complete_input_prof_data[iblock][1][0,0]) < TINY_NUMBER
                       and (n['species'][isp]['z'] - complete_input_prof_data[iblock][1][0,1]) < TINY_NUMBER):
                        # due to the not-so-great file format of
                        # input.prof there is no header line specified
                        # to address the block containing the
                        # background profiles. Instead, it is the next block!
                        n['species'][isp]['rln'] = complete_input_prof_data[iblock+1][1][args.radial_ix,gkwlib.input_prof_background_prof_columns.index("rln")]
                        n['species'][isp]['dens'] = complete_input_prof_data[iblock+1][1][args.radial_ix,gkwlib.input_prof_background_prof_columns.index("dens")]
                
        else:
            # I do not want to compute any complicated formula
            # here. So let's use profiles output by GKW, based on the
            # parameters. This means you must run GKW (just one step,
            # if you would like) in order to be able to convert to a
            # local run.
            f = h5py.File(args.gkwdatah5_file, "r+")
            n['species'][isp]['rln'] = f['/diagnostic/diagnos_rad/background_profiles/'+'background_rln'].value[isp,args.radial_ix]
            n['species'][isp]['dens'] = f['/diagnostic/diagnos_rad/background_profiles/'+'background_dens'].value[isp,args.radial_ix]
            n['species'][isp]['uprim'] = f['/diagnostic/diagnos_rad/background_profiles/'+'background_uprim'].value[isp,args.radial_ix]

    for isp in range(n['gridsize']['number_of_species']):
        if(n['species'][isp]['temp_prof_type'] == 'file'):
            # I can read input.prof without running GKW
            complete_input_prof_data = gkwlib.read_input_prof_data(args.inputprof_file)
            for iblock in range(len(complete_input_prof_data)):
                if(complete_input_prof_data[iblock][0].startswith("Species")):
                    # check if this is indeed the input.prof entry for this species:
                    if(abs(n['species'][isp]['mass'] - complete_input_prof_data[iblock][1][0,0]) < TINY_NUMBER
                       and (n['species'][isp]['z'] - complete_input_prof_data[iblock][1][0,1]) < TINY_NUMBER):
                        # due to the not-so-great file format of
                        # input.prof there is no header line specified
                        # to address the block containing the
                        # background profiles. Instead, it is the next block!
                        n['species'][isp]['rlt'] = complete_input_prof_data[iblock+1][1][args.radial_ix,gkwlib.input_prof_background_prof_columns.index("rlt")]
                        n['species'][isp]['temp'] = complete_input_prof_data[iblock+1][1][args.radial_ix,gkwlib.input_prof_background_prof_columns.index("temp")]
                
        else:
            for isp in range(n['gridsize']['number_of_species']):
                # I do not want to compute any complicated formula
                # here. So let's use profiles output by GKW, based on the
                # parameters. This means you must run GKW (just one step,
                # if you would like) in order to be able to convert to a
                # local run.
                f = h5py.File(args.gkwdatah5_file, "r+")
                n['species'][isp]['rlt'] = f['/diagnostic/diagnos_rad/background_profiles/'+'background_rlt'].value[isp,args.radial_ix]
                n['species'][isp]['temp'] = f['/diagnostic/diagnos_rad/background_profiles/'+'background_temp'].value[isp,args.radial_ix]
                n['species'][isp]['uprim'] = f['/diagnostic/diagnos_rad/background_profiles/'+'background_uprim'].value[isp,args.radial_ix]


    # pick geometry parameters at this radial position
    if(n['geom']['prof_type'] == 'file'):
        # The profiles of miller parameters are currently not output into the gkwdata.h5 file.
        # This is why they are read from input.prof here.
        complete_input_prof_data = gkwlib.read_input_prof_data(args.inputprof_file)
        iblock_geometry = gkwlib.get_input_prof_block_index(complete_input_prof_data, "Geometry")
        ncols = complete_input_prof_data[iblock_geometry][1].shape[1]

        n['geom']['q'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("q")]
        n['geom']['shat'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("shat")]


        if(n['geom']['geom_type'] == 'miller'):
            n['geom']['kappa'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("kappax")]
            n['geom']['skappa'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("skappax")]
            n['geom']['delta'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("deltax")]
            n['geom']['sdelta'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("sdeltax")]
            n['geom']['square'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("squarex")]
            n['geom']['ssquare'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("ssquarex")]
            n['geom']['zmil'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("Zmilx")]
            n['geom']['drmil'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("dRmilx")]
            n['geom']['dzmil'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("dZmilx")]
            n['geom']['gradp'] = complete_input_prof_data[iblock_geometry][1][args.radial_ix,gkwlib.input_prof_miller_columns.index("gradpx")]
    else:
        f = h5py.File(args.gkwdatah5_file, "r+")
        n['geom']['q'] = f['/diagnostic/diagnos_rad/background_profiles/'+'qx'].value[args.radial_ix]
        n['geom']['shat'] = f['/diagnostic/diagnos_rad/background_profiles/'+'shatx'].value[args.radial_ix]

    if(n['geom']['prof_type'] == 'file'):
        n['geom']['prof_type'] = 'none'

    # calculate the corresponding kthrho for the local run.
    #kzeta = 2*np.pi * n['spcgeneral']['rhostar'] * n

    # one probably wants to reduce gridsize, too. Put some values.
    if(n['control']['non_linear']):
        n['gridsize']['nx'] = 167
    else:
        n['gridsize']['nx'] = 5
    n['gridsize']['n_s_grid'] = 20

    # adopt the size of the velocity space, with same resolution. This
    # might be more than is actually necessary, but I see no other way,
    # since it is determined by physics.
    # n['gridsize']['n_mu_grid'] = #do not change
    # n['gridsize']['n_vpar_grid'] = #do not change
    # n['gridsize']['vpmax'] = #do not change
    # n['gridsize']['mumax'] = #do not change

    #for nonspectral fluxtube runs, one must specify
    n['gridsize']['lx'] = 140 #FIXME do this properly

    # one can probably use a much larger timestep
    n['control']['dtim'] = 1.5e-3
    n['control']['nl_dtim_est'] = True

    # it is recommended, to compute the free energy balance diagnostic for local runs
    n['diagnostic']['lcalc_energetics'] = True

    # if this is local, reduce naverage, to avoid floating overflow
    if (not n['control']['non_linear']):
        n['control']['naverage'] = 100

    print("!", file=output)
    print("! You may want to change the dissipation coefficients manually.", file=output)
    print("!", file=output)
    
    print(nlist.dump(), file=output)
