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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Change the electron/ion mass ratio of a GKW global input file, and the associated input.prof file.")
    parser.add_argument('--input-file',
                        help="GKW local input file (default: input.dat)", 
                        type=str,
                        default='input.dat'
    )
    parser.add_argument('--inputprof-file',
                        help="GKW input.prof file (default: input.prof)", 
                        type=str,
                        default='input.prof'
    )
    ELECTRON_SPECIES = "electron"
    ION_SPECIES = "ion"
    c = [ELECTRON_SPECIES,ION_SPECIES]
    parser.add_argument('--species',
                       help="which species to change. default 'electron' changes electron mass, and leaves ion mass as it is", 
                       type=str,
                       choices=c,
                       default=ELECTRON_SPECIES
    )
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--mass',
                        help="new mass", 
                        type=float,
                        default=None
    )
    group.add_argument('--mass-factor',
                        help="factor, to change mass", 
                        type=float,
                        default=None
    )
    group.add_argument('--mass-ratio',
                        help="change mass, to meet given mass ratio of first species and electron species", 
                        type=float,
                        default=None
    )
    group.add_argument('--larmor',
                        help="change mass, to meet given larmor radius", 
                        type=float,
                        default=None
    )
    group.add_argument('--larmor-factor',
                        help="change mass, to meet larmor radius times this factor.", 
                        type=float,
                        default=None
    )
    group.add_argument('--larmor-ratio',
                        help="change mass, to meet given larmor radius ratio of first species and electron species", 
                        type=float,
                        default=None
    )

    parser.add_argument('--output-folder',
                        help="output folder (will be created if not existent). This is where the new input.dat and input.prof will be put. You may use {species}, {mass}, {larmor}, {factor} and {ratio} to format parameters into the folder name. (default folder name depends on other parameters)", 
                        type=str,
                        default=None
    )
    
    args = parser.parse_args()

    if(args.mass_factor is not None):
        factor = args.mass_factor
    else:
        factor = args.larmor_factor

    if(args.mass_ratio is not None):
        ratio = args.mass_ratio
    else:
        ratio = args.larmor_ratio
    
    if(args.output_folder is None):
        if(args.mass is not None):
            args.output_folder = '{species}_mass_{mass}'
        elif(args.mass_factor is not None):
            args.output_folder = '{species}_mass_factor_{factor}'
        elif(args.mass_ratio is not None):
            args.output_folder = '{species}_mass_ratio_{ratio}'
        elif(args.larmor is not None):
            args.output_folder = '{species}_larmor_{larmor}'
        elif(args.larmor_factor is not None):
            args.output_folder = '{species}_larmor_factor_{factor}'
        elif(args.larmor_ratio is not None):
            args.output_folder = '{species}_larmor_ratio_{ratio}'
    foldername = args.output_folder.format(species=args.species, mass=args.mass, larmor=args.larmor,
                                           factor=factor, ratio=ratio)

    # parse the file into a Namelist object, and take the groups member of that.
    nlist = namelist_python.read_namelist_file(args.input_file)
    n = nlist.groups
    # that object is an ordered case-insensitive dictionary

    complete_input_prof_data = gkwlib.read_input_prof_data(args.inputprof_file)

    # find the target species in the input file namelist
    ion_isp = None
    electron_isp = None
    for isp in range(n['gridsize']['number_of_species']):
        if(n['species'][isp]['Z'] > 0):
            if(ion_isp is None):
                ion_isp = isp
                if(args.species == ION_SPECIES):
                    # simply consider the first ion species found
                    target_isp = isp
        else:
            if(electron_isp is None):
                electron_isp = isp
                if(args.species == ELECTRON_SPECIES):
                    target_isp = isp
    mass = n['species'][target_isp]['mass']
    Z = n['species'][target_isp]['Z']
    # try:
    #     dens = n['species'][target_isp]['dens']
    # except:
    #     dens = 1.0
    try:
        T = n['species'][target_isp]['temp']
    except:
        T = 1.0

    print("mass: %f" % mass)
    print("charge Z: %f" % Z)
    print("temperature T: %f" % T)
    
    # find the target species in the input.prof file
    for iblock in range(len(complete_input_prof_data)):
        if(complete_input_prof_data[iblock][0].startswith("Species")):
            #Species: m , Z, ngrid, Tgrid
            #print(complete_input_prof_data[iblock][1][0,:])
            tiny = 1.0e-3
            if(np.abs(mass - complete_input_prof_data[iblock][1][0,0]) < tiny
               and np.abs(Z - complete_input_prof_data[iblock][1][0,1]) < tiny):
                target_iblock = iblock
                print("species has been found in input.prof")

    print("target species: %d" % target_isp)
    print("target species input.prof block: %d" % target_iblock)
    print("electron species: %d" % electron_isp)
    print("ion species: %d" % ion_isp)
    print()
    #T = 1/2 m v_th^2
    #T_ref T_N = 1/2 m_R m_ref (v_thref v_R)^2
    # in GKW, we define T_ref = 1/2 m_ref v_thref^2, hence
    #T_N = m_R v_R^2
    #v_R = sqrt(T_N/m_R)
    
    # rho = m v_th/(eB)
    # rho_N rho_ref = m_R m_ref v_R v_thref/(e B_N B_ref)
    # in GKW, we define rho_ref = (m ref v thref)/(eB ref), hence
    # rho_N = m_R v_R/B_N
    # and use the temperature and say that B_N \approx 1
    # rho_N = m_R v_R/B_N
    # rho_N = m_R sqrt(T_N/m_R) = sqrt(m_R T_N)

    rho = np.sqrt(mass * T)
    print("rho: %f" % rho)

    electron_mass = n['species'][electron_isp]['mass']
    ion_mass = n['species'][ion_isp]['mass']
    print("electron mass: %f" % electron_mass)
    print("ion mass: %f" % ion_mass)

    electron_rho = np.sqrt(electron_mass * n['species'][electron_isp]['temp'])
    ion_rho = np.sqrt(ion_mass * n['species'][ion_isp]['temp'])

    print("mass ratio: %f" % (ion_mass/electron_mass))
    print("rho ratio: %f" % (ion_rho/electron_rho))
    print()

    if(args.mass is not None):
        new_mass = args.mass
    elif(args.mass_factor is not None):
        new_mass = mass * args.mass_factor
    elif(args.mass_ratio is not None):
        if(args.species == ION_SPECIES):
            new_mass = electron_mass * args.mass_ratio
        else:
            new_mass = ion_mass / args.mass_ratio
    elif(args.larmor is not None):
        new_mass = args.larmor**2/T
    elif(args.larmor_factor is not None):
        new_mass = (rho * args.larmor_factor)**2/T
    elif(args.larmor_ratio is not None):
        if(args.species == ION_SPECIES):
            new_mass = (electron_rho * args.larmor_factor)**2/T
        else:
            new_mass = (ion_rho / args.larmor_factor)**2/T

    print("new mass: %f" % new_mass)
    new_rho = np.sqrt(new_mass * T)
    print("new rho: %f" % new_rho)
    
    new_ion_mass = ion_mass
    new_electron_mass = electron_mass
    new_ion_rho = ion_rho
    new_electron_rho = electron_rho
    if(args.species == ION_SPECIES):
        new_ion_mass = new_mass
        new_ion_rho = new_rho
    else:
        new_electron_mass = new_mass
        new_electron_rho = new_rho
        
    print("new mass ratio: %f" % (new_ion_mass/new_electron_mass))
    print("new larmor ratio: %f" % (new_ion_rho/new_electron_rho))

    # put the desired mass value into the datastructures for input.prof and input.dat
    new_block = complete_input_prof_data[target_iblock][1]
    new_block[0,0] = new_mass
    complete_input_prof_data[target_iblock] = (complete_input_prof_data[target_iblock][0],new_block)
    n['species'][target_isp]['mass'] = new_mass

    print()
    
    # now write the new mass back to a file
    gkwlib.make_run_folder(foldername, input_file_namelist=nlist, input_prof_data=complete_input_prof_data)
    print("Data has been written to %s/" % foldername)
    

