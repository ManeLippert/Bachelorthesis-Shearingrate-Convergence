#!/usr/bin/env python3

# the namelist_python module used here is available at
# https://github.com/stegro/namelist_python
# or a copy in the $GKW_HOME/python folder.

import sys
try:
    import namelist_python
except:
    # if the GKW_HOME/python folder is not contained in the PYTHONPATH
    import os
    sys.path.append(os.path.join(os.environ['GKW_HOME'],'python'))
    import namelist_python
    
from math import pi
import numpy as np

if __name__=='__main__':

    inputfile_name = 'input.dat'
    outputfile_name = 'input.dat.converted'
    gkwdatah5_name = 'gkwdata.h5'
    geomdat_name = 'geom.dat'
    
    print("""
 Usage:
    %s [inputfile [outputfile [gkwdata.h5|geom.dat]]]

 This script reads a GKW inputfile (spectral, flux_tube) and
 converts it from modebox=F to T.

 If the inputfile specifies
    kr_type='chin'
 then this script needs to read some geometrical quantities from
 a gkwdata.h5 or geom.dat file.
    
 For the more complicated geometries 'miller' and 'chease' specifying
 a gkwdata.h5 or geom.dat file is obligatory.

 Default inputfile name: %s
 Default outputfile name: %s
 Default geometry data file name: %s or %s

    """ % (sys.argv[0], inputfile_name, outputfile_name, gkwdatah5_name, geomdat_name))

    try:
        inputfile_name = sys.argv[1]
        outputfile_name = sys.argv[2]
        if(outputfile_name.endswith("h5") or outputfile_name.startswith('geom')):
            raise SystemExit('Output file %s resembles an HDF5 file or a geom.dat file - probably a mistake?' % outputfile_name)
        gkwdatah5_name = sys.argv[3]
        
    except SystemExit as err:
        print(err)
        exit(1)
    except:
        # no argument given
        pass

    # parse the file into a Namelist object, and take the groups member of that.
    nlist = namelist_python.read_namelist_file(inputfile_name)
    n = nlist.groups
    # that object is an ordered case-insensitive dictionary

    # set some default values if they are not there
    if('geom_type' not in n['geom']):
        n['geom']['geom_type'] = 's-alpha'
    if('q' not in n['geom']):
        n['geom']['q'] = 1.0
    if('eps' not in n['geom']):
        n['geom']['eps'] = 1.0
    if('spectral_radius' not in n['control']):
        n['control']['spectral_radius'] = True
    if('flux_tube' not in n['control']):
        n['control']['flux_tube'] = True
    if('mode_box' not in n['mode']):
        n['mode']['mode_box'] = False
    if('kr_type' not in n['mode']):
        n['mode']['kr_type'] = 'chin'
        if('chin' not in n['mode']):
            n['mode']['chin'] = '0.0'

    if(n['mode']['kr_type'] == 'kr'):
        if('krrho' not in n['mode']):
            n['mode']['krrho'] = '0.0'

        
    if(not n['control']['spectral_radius']):
        raise SystemExit('Error in inputfile: must be spectral_radius = true')
    
    if(not n['control']['flux_tube']):
        raise SystemExit('Error in inputfile: must be flux_tube = true')

    print("The inputfile has modebox = %s" % str(n['mode']['mode_box']))
    print(" -> Convert this to modebox = %s" % str(not str(n['mode']['mode_box'])))
    print()
    
    if(n['mode']['mode_box']):
        # generate an output file with:
        n['mode']['mode_box'] = False

        if(n['grid']['nperiod'] != 1):
            raise ValueError('error in inputfile: if mode_box=T then must be nperiod = 1')
        
        raise NotImplementedError('This conversion is not yet implemented')
        
    else:
        # generate an output file with:
        n['mode']['mode_box'] = True
        
        if(n['gridsize']['nx'] != 1):
            raise ValueError('Error in inputfile: must have nx = 1')

        nturns = (2*n['gridsize']['nperiod'] - 1)
        points_per_turn = n['gridsize']['n_s_grid'] / nturns
        # put the new n_s_grid into the datastructure:
        n['gridsize']['n_s_grid'] = points_per_turn

        print('The given file specifies nperiod=%d, i.e. %d turns of the fluxtube.' % (n['gridsize']['nperiod'], nturns))
        all_or_some_nturns = None
        while(all_or_some_nturns not in ["all","some"]):
            all_or_some_nturns = input("""Do you want at least %d turns for *every* ky mode (ENTER "all") in the computational box,
 or is it enough if only the largescale modes do %d turns (ENTER "some")? """ % (nturns, nturns))

        #for the ky=0 mode, every kx connects to itself at the parallel boundary condition.
        #for other ky modes, the kx connect to kx+delta_kx.
        # delta_kx depends on ky.

        eps = n['geom']['eps']
        
        if(n['geom']['geom_type'] == 's-alpha'):
            q = n['geom']['q']
            shat = n['geom']['shat']
            kthnorm = q  / ( 2 * pi * eps)
        elif(n['geom']['geom_type'] == 'circ'):
            q = n['geom']['q']
            shat = n['geom']['shat']
            kthnorm = sqrt(1.0 + (q/eps)**2*(1-eps**2)) /(2 * pi * (1.0+eps))
        elif(n['geom']['geom_type'] in ['miller', 'chease']):
            try:
                import h5py
                f = h5py.File(gkwdatah5_name, "r+")
                q = f["/geom/q"][0]
                print("Obtained q = %e from the file %s" % (q, gkwdatah5_name))
                shat = f["/geom/shat"][0]
                print("Obtained shat = %e from the file %s" % (shat, gkwdatah5_name))
                kthnorm = f["/geom/kthnorm"][0]
                print("Obtained kthnorm = %e from the file %s" % (kthnorm, gkwdatah5_name))
            except Exception as err:
                print("HDF5 file %s could not be read." % gkwdatah5_name)
                print(err)
                try:
                    raise NotImplementedError('reading geom.dat is not yet implemented')
                except:
                    raise SystemExit("The file %s could not be read." % geomdat_name)
                
            
        nmod = n['gridsize']['nmod']
        answer = input("Enter a new nmod if you want (ENTER to keep nmod = %d): " % nmod)
        try:
            nmod = int(answer)
        except:
            print(" -> nmod = %d is used for the new box" % nmod)

        try:
            krhomax = max(n['mode']['kthrho'])
        except:
            # just one kthrho value given
            krhomax = n['mode']['kthrho']
        print("Taking the largest of the kthrho values would mean krhomax = %f" % krhomax)
        try:
            krhomax = float(input("Enter a new krhomax if you want (ENTER to keep krhomax = %f): " % krhomax))
            print("WARNING, FIXME: With arbitrary nmod and krhomax you may not have the orig. krho in the new box.")
            print("WARNING, FIXME: This script does not check this error/help you making a good choice, at the moment.")
            print()
        except:
            pass
        n['mode']['krhomax'] = krhomax
        print(" -> krhomax = %f" % krhomax)

        print("Summary of the ky modes:")
        print("original box krho:")
        print(n['mode']['kthrho'])
        print("new box krho:")
        if(nmod > 1):
            print(np.linspace(0,krhomax,nmod))
        else:
            print(krhomax)

        # distinguish the cases nmod = 1 and nmod > 1
        krho_min = krhomax / max(nmod-1, 1)
        ky_min = krho_min / kthnorm


        kxrh = None
        if(n['mode']['kr_type'] == 'chin'):
            if(n['gridsize']['nmod']  > 1):
                raise ValueError("chin can only be used for nmod=1")
            chin = n['mode']['chin']
            print(" Given chin = %f" % chin)

            try:
                import h5py
                f = h5py.File(gkwdatah5_name, "r+")
                pol_angle = f["/geom/poloidal_angle"]
                g_eps_zeta = f["/geom/g_eps_zeta"]
                g_eps_eps = f["/geom/g_eps_eps"]
                print("Obtained pol_angle, g_eps_zeta and g_eps_eps from the file %s" % (gkwdatah5_name))
            except Exception as err:
                print("HDF5 file %s could not be read." % gkwdatah5_name)
                print(err)
                try:
                    raise NotImplementedError('reading geom.dat is not yet implemented')
                except:
                    raise SystemExit("The file %s could not be read." % geomdat_name)

            # for each kthrho:
            # ky = kthrho / kthnorm
            #kx = - ky * g_psizeta / g_psipsi

            # Note that we have
            # krho_min = krhomax = krho(1)
            # since we consider only the case nmod=1

            # compute the value of kx to minimise k_perp @ pol_angle=CHIN
            r_tiny = 1e-10
            if(abs(chin) < r_tiny):
                kxrh = 0
            else:
                # n_s_grid = len(pol_angle)
                # dum_angle = n_s_grid * [None]
                # kx_on_sgrid = n_s_grid * [None]
                # for i in range(n_s_grid):
                #     dum_angle[i] = pol_angle[i] - chin
                #     kx_on_sgrid[i] = - krho_min * g_eps_zeta[i] / g_eps_eps[i]
                dum_angle = pol_angle.value - chin
                kx_on_sgrid = - krho_min * g_eps_zeta.value / g_eps_eps.value

                #print("dum_angle: " + str(np.abs(dum_angle)))
                i_chin = np.argmin(np.abs(dum_angle))
                if (abs(pol_angle[i_chin] - chin) < r_tiny):
                    # if the given angle chin coincides with the poloidal angle corresponding to an s gridpoint
                    kxrh = kx_on_sgrid[i_chin]
                else:
                    # if the chin is somewhere inbetween
                    # then linear interpolation:
                    if(pol_angle[i_chin] > chin):
                        j_chin = i_chin - 1
                    else:
                        j_chin = i_chin + 1
                    kxrh = (kx_on_sgrid[i_chin] - kx_on_sgrid[j_chin]) / (pol_angle[i_chin] - pol_angle[j_chin]) * (chin - pol_angle[i_chin]) + kx_on_sgrid[i_chin]
                # In s-alpha the found value corresponds to
                # kxrh(1) = - chin* abs(q*shat*krho(imod=1)) / (2*pi*eps)
                if(n['geom']['geom_type'] == 's-alpha'):
                    kxrh_comparison = - chin* q*abs(shat)*krho_min / (2*pi*eps)
                    print('Comparison of kxrh (computed from chin) with an analytical s-alpha expression:')
                    print(' kxrh s-alpha analytical: %f' % kxrh_comparison)
                    print(' difference computed - analytical: %f -> %s' %
                          (kxrh - kxrh_comparison,
                           abs(kxrh - kxrh_comparison) < r_tiny and "like it should be" or "code mistake?" )
                    )
            print(' kxrh computed: %f' % kxrh)
            
            # if(chin != 0.0):
            #     raise NotImplementedError('calculating kx from chin is not yet implemented for chin!=0')
        elif(n['mode']['kr_type'] == 'kr'):
            kxrh = n['mode']['krrho']
        else:
            raise ValueError("kr_type has unknown value: %s" % n['mode']['kr_type'])

        kx = kxrh / kthnorm
        print(' -> striving to have kxrh = %f (kx = %f) on the grid' % (kxrh,kx))
        if(kx < 0):
            kxrh = abs(kxrh)
            kx = abs(kx)
            print(' -> striving to have kxrh = %f (kx = %f) on the grid' % (abs(kxrh),abs(kx)))
        
        # the kx value after one poloidal turn depends on the ky of the mode
        #kxplus = kx + abs(q*shat*ky/eps)
        # to be able to connect the kx for each ky mode, construct the grid based on the smallest nonzero ky
        kxplus_min = abs(q*shat*ky_min / eps)
        
        # ikxspace then allows to increase the radial resolution further:
        #kxspace = kxplus_min / real(ikxspace)
        # So choose ikxspace:

        # If kx == 0 then we can choose any ikxspace and the mode from the orig box will always be in the system
        # (it is the/a mode with ixzero). We can freely choose ikxspace to get a certain aspect ratio of
        # the computational box.
        # However, if kx!=0 then we may want to choos ikxspace in such a way that the new kx grid contains
        # this kx of the original mode. If we try to get aspect ratio=1 then we may not exactly have the
        # original kx in the system.
        # It may also need very high ikxspace to have accurately certain kx.

        #lyn = 2*pi/(krho_min);
        #kx_min = kxplus_min/ikxspace
        
        #lx = 2*pi/kx_min
        # we want to:
        # lyn = lxn
        # krho_min = kx_min
        #          = kxplus_min/ikxspace;
        # -> ikxspace = kxplus_min/krho_min
        ikxspace_for_aspratio1 = round(kxplus_min / krho_min)

        # to make sure the kx of the orig box is contained in the new box as
        # best as possible, we want to
        #  1) try feasible ikxspace values
        #  2) for each ikxspace value, compute the minimum of the differences of all
        #     kx gridpoints with the orig. kx
        ikxspace_uplim = 20
        difference = np.ones(ikxspace_uplim)*np.inf
        print(" ikxspace | asp. ratio | kxspace       | kx error")
        print("----------+------------+---------------+---------------")
        for ikxspc in range(1,ikxspace_uplim):
            # the kx gridspacing with that value would be:
            kxspace = kxplus_min / (1.0*ikxspc)
            modulo = kx % kxspace
            difference[ikxspc] = min(modulo, kxspace-modulo)
            kx_min = (kxplus_min/ikxspc)
            asp_ratio = kx_min/krho_min
            print("   %2d    | %11f | %13.5f | %13.5f" % (ikxspc, asp_ratio, kxspace, difference[ikxspc]))            
            
        #  3) choose ikxspace which  minimizes that minimum and is as small as
        #     possible (i.e. dont take ikxspace=8 if 4 does it, too).
        print()
        print(" - Just trying to approximate aspect ratio 1 leads to : ikxspace = %d" % ikxspace_for_aspratio1)
        print()
        
        answer = input("Enter ikxspace now manually (just ENTER to try to choose ikxspace=%d based on aspect ratio only) : "%ikxspace_for_aspratio1)
        try:
            n['mode']['ikxspace'] = int(answer)
        except ValueError as err:
            n['mode']['ikxspace'] = ikxspace_for_aspratio1

        print(" -> ikxspace = %d is proposed" % n['mode']['ikxspace'])
        kx_min = kxplus_min/n['mode']['ikxspace']
        print(" -> leads to aspect ratio = lyn/lx = kx_min/krho_min = %f " % (kx_min/krho_min))

        # now choose the number of radial modes according to the turns
        if(all_or_some_nturns == "some"):
            # this is probably what you expect:
            n['gridsize']['nx'] = 1 + (nturns-1)*n['mode']['ikxspace']
            # I think then the kx!=0 modes may end with one turn
            # less than than the kx=0 mode. Therefore:
            n['gridsize']['nx'] = nturns*n['mode']['ikxspace']
        elif(all_or_some_nturns == "all"):
            # kxplus is proportional to ky - and indeed equal to kxplus_min * imod.
            # this means that at the parallel boundary the radial modes connect like this
            #   (ky=0, kx) -> (ky=0, kx)
            #   (ky!=0,kx) -> (ky!=0,kx+kxplus_min*imod)
            #       => the largest kx needed is nperiod * kxplus_min*nmod
            #       => the number of kx gridpoints is then nturns * nmod * ikxspace
            n['gridsize']['nx'] = nturns*n['mode']['ikxspace']*nmod

        # make sure that nx is an odd number:
        n['gridsize']['nx'] += n['gridsize']['nx'] % 2 == 0 and 1 or 0
        
        print(" -> radial gridpoints nx=%d" % n['gridsize']['nx'])

        # put the new nmod into the config
        n['gridsize']['nmod'] = nmod

        #remove settings without meaning now:
        n['gridsize'].pop('nperiod')
        n['mode'].pop('kthrho')
        n['mode'].pop('kr_type')
        if('chin' in n['mode']):
            n['mode'].pop('chin')
        if('krrho' in n['mode']):
            n['mode'].pop('krrho')
        
        with open(outputfile_name, 'w') as f:
            f.write(nlist.dump())
            print()
            print("Output has been written to file %s" % outputfile_name)

