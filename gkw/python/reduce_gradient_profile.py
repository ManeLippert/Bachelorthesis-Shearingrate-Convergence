#!/usr/bin/env python3

import sys
try:
    import gkwlib
except:
    # if the GKW_HOME/python folder is not contained in the PYTHONPATH
    import os
    sys.path.append(os.path.join(os.environ['GKW_HOME'],'python'))
    import gkwlib

import os
import shutil
import numpy as np
import scipy.interpolate
import scipy.integrate
import copy

import matplotlib.pyplot as plt

def clip_electron_background(complete_input_prof_data, background_dset, max_value):
    reduction_ratio_profile = None
    reduction_ratio_grad_profile = None
    for iblock in range(len(complete_input_prof_data)):
        if(complete_input_prof_data[iblock][0].startswith("Species")):
            # check if this is indeed the input.prof entry for this species:
            mass = complete_input_prof_data[iblock][1][0,0]
            Z = complete_input_prof_data[iblock][1][0,1]
            if(Z < 0):
                # due to the not-so-great file format of
                # input.prof there is no header line specified
                # to address the block containing the
                # background profiles. Instead, it is the next block!
                xgrid = complete_input_prof_data[iblock+1][1][:,gkwlib.input_prof_background_prof_columns.index("xgr")]
                icol = gkwlib.input_prof_background_prof_columns.index(background_dset)
                profile = complete_input_prof_data[iblock+1][1][:,icol]
                profile = np.clip(profile, -max_value, +max_value)
                # obtain a smoothed cubic spline representation of the curve
                smoothing_cond = 1.0e0
                tck = scipy.interpolate.splrep(xgrid,profile, k=3, s=smoothing_cond)
                # evaluate the profile
                new_profile = scipy.interpolate.splev(xgrid, tck, der=0)
        
                if(background_dset.startswith("rl")):
                    # having modified the gradient, now reintegrate the profile
                    
                    scalar_profile = complete_input_prof_data[iblock+1][1][:,icol-2]
                    leftmost_value = scalar_profile[0]
                    print("Re-integrate the profile starting from leftmost value %f..." % leftmost_value)
                    
                    
                    # Right now, there seems to be a bug in cumtrapz
                    # regarding the initial value. This is why one
                    # must take initial=0 and add the constant.
                    
                    # Note that what we specify as paramater is an inverse normalised
                    # gradient length 1/L_{T,N} = R/L_T = -1/T dT/dPsi .
                    # Thus, to obtain T(Psi), integrate this:
                    # dT/dPsi = - T R/L_T
                    new_scalar_profile = scipy.integrate.cumtrapz(- scalar_profile * new_profile, x=xgrid, initial=0) + leftmost_value
                    reduction_ratio_grad_profile = (new_scalar_profile * new_profile)/(scalar_profile * complete_input_prof_data[iblock+1][1][:,icol])
                    reduction_ratio_profile = new_scalar_profile/complete_input_prof_data[iblock+1][1][:,icol-2]
                    
                    complete_input_prof_data[iblock+1][1][:,icol-2] = new_scalar_profile
                    

                else:
                    # having modified the profile, now rederive the gradient
                    
                    new_grad_profile = -scipy.interpolate.splev(xgrid, tck, der=1)/new_profile

                    reduction_ratio_grad_profile = (new_scalar_profile * new_profile)/(complete_input_prof_data[iblock+1][1][:,icol+2] * complete_input_prof_data[iblock+1][1][:,icol])
                    reduction_ratio_profile = new_profile/complete_input_prof_data[iblock+1][1][:,icol]
                    
                    complete_input_prof_data[iblock+1][1][:,icol+2] = new_grad_profile
                complete_input_prof_data[iblock+1][1][:,icol] = new_profile

    # print("reduction_ratio_profile:")
    # print(reduction_ratio_profile)
    # print("reduction_ratio_grad_profile:")
    # print(reduction_ratio_grad_profile)
    # print("reduction_ratio diff:")
    # print(reduction_ratio_profile - reduction_ratio_grad_profile)

    if(background_dset in ("dens", "rln")):
        for iblock in range(len(complete_input_prof_data)):
            if(complete_input_prof_data[iblock][0].startswith("Species")):
                # check if this is indeed the input.prof entry for this species:
                mass = complete_input_prof_data[iblock][1][0,0]
                Z = complete_input_prof_data[iblock][1][0,1]
                dens_icol = gkwlib.input_prof_background_prof_columns.index("dens")
                rln_icol = gkwlib.input_prof_background_prof_columns.index("rln")
                if(Z > 0):
                    # this is an ion species. Modify the density profile to ensure quasineutrality.
                    # \sum_is signz_G(is)*de_G(ix,is)*fp_G(ix,is) = 0
                    # \sum_is signz_G(is)*de_G(ix,is) = 0
                    # We just have changed both electron density and electron density gradient.
                    # Therefore, change the ion densities, too:
                    complete_input_prof_data[iblock+1][1][:,dens_icol] *= reduction_ratio_profile
                    complete_input_prof_data[iblock+1][1][:,rln_icol] *= reduction_ratio_grad_profile/reduction_ratio_profile
                    
    return complete_input_prof_data


if __name__ == "__main__":
    inputprof_file = "input.prof"
    complete_input_prof_data = gkwlib.read_input_prof_data(inputprof_file)

    # 1) reduce electron temp gradient profile:

    for max_value in [12.0,11.0,10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0]:
        work_data = copy.deepcopy(complete_input_prof_data)
        work_data = clip_electron_background(work_data, "rlt", max_value)

        dirname = "clipped_rlt_%f" % max_value
        gkwlib.make_run_folder(dirname, orig_input_file="input.dat", input_prof_data=work_data)
        print()

    # 2) reduce electron density gradient profile:

    for max_value in [12.0,11.0,10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0]:
        work_data = copy.deepcopy(complete_input_prof_data)
        work_data = clip_electron_background(work_data, "rln", max_value)

        # ensure quasi-neutrality
        #\sum_is signz_G(is)*de_G(ix,is)*fp_G(ix,is) = 0

        dirname = "clipped_rln_%f" % max_value
        gkwlib.make_run_folder(dirname, orig_input_file="input.dat", input_prof_data=work_data)
        print()

    # 3) reduce electron density AND temp gradient profile:

    for max_value in [12.0,11.0,10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0]:
        work_data = copy.deepcopy(complete_input_prof_data)
        work_data = clip_electron_background(work_data, "rln", max_value)
        work_data = clip_electron_background(work_data, "rlt", max_value)
        
        dirname = "clipped_rln_rlt_%f" % max_value
        #gkwlib.make_run_folder(dirname, orig_input_file="input.dat", input_prof_data=work_data)
        gkwlib.make_run_folder(dirname, orig_input_file="input.dat", input_prof_data=work_data)
        print()

        
