### to remove later ###
import sys
sys.path.append('/home/yann/projects/gkdb/gkdb_gitlab/gkdb/tools/')
import matplotlib.pyplot as plt
#####
import f90nml
import re
import numpy as np
import os
from scipy.interpolate import interp1d
import FS_param as FS # in gkdb/gkdb/tools


class gkwrun:
  """ Class containing data from GKW simulations 
      gkwrun.storage          -> path of GKW input and outputs files
                                 filled by set_storage
            .input            -> content of GKW input files input.dat and input_out.dat
            .input_out           filled by get_inputs
            .geom             -> content of geom.dat output file
                                 filled by get_geom
            .grids            -> following grids: time, s, kthrho, krrho, keps, kzeta, file_count
                                 filled by get_grids
            .info             -> content of GKW 'out' file (as a string)
                                 filled by get_info
            .eigenvalues      -> mode growth rate and real frequency
                                 filled by get_eigenvalues
            .collisions       -> collision frequencies for all species
                                 filled by get_collisions
            .cfdens           -> background density and density gradient when poloidal asymmetry
                                 filled by get_cfdens
            .fluxes           -> fluxes in co-moving and laboratory frames
            .fluxes_lab          filled by get_fluxes
            .parallel         -> parallel structure of fields and moments from parallel.dat
                                 filled by get_parallel
            .kykxs            -> fields and moments from kykxs diagnostics
                                 filled by get_kykxs
     If data is missing, empty dict, strings or arrays are returned
     1) To add
            .spectra          -> fields, fluxes and vorticity spectra          
     2) Implement verbose / silent modes + store log of warning error messages in the data structure (list of strings?)

  """


  def set_storage(self,flpth,file_style='multiple_folders',flnm=None):
    """ Get the path of GKW input and output files 
    Inputs
      file_style   either 'single_folder' or 'multiple_folders' file storage
      flpth        path of the project main folder (multi-folders file storage) 
                   or path to the run folder (single folder file storage))
      flnm         name of the run (used only for multi-folders file storage)
    Outputs
      gkwrun.storage['path']    dict containing the path to the various GKW inputs and outputs
    """
    # checks
    if flpth[-1]!='/':
      flpth=flpth+'/'    
    assert (file_style in ('single_folder','multiple_folders')),\
           "Accepted values for 'file_style' are 'single_folder' or 'multiple_folders'"
    assert not((file_style=='multiple_folders') & (flnm==None)),\
           "flnm required for file_style='multiple_folders'"
    # store input values
    self.storage={}
    self.storage['style']=file_style
    self.storage['flpth']=flpth
    self.storage['flnm']=flnm
    # store path to GKW inputs and outputs
    self.storage['path']={}
    if self.storage['style']=='single_folder':
      self.storage['path']['input']=flpth+'input.dat'
      # to be completed
    elif self.storage['style']=='multiple_folders':      
      self.storage['path']['input']=flpth+'input/'+flnm
      self.storage['path']['input_out']=flpth+'input_out/'+flnm
      self.storage['path']['file_count']=flpth+'grids/file_count/'+flnm
      self.storage['path']['kthrho']=flpth+'grids/krho/'+flnm
      self.storage['path']['keps']=flpth+'grids/kxrh/'+flnm
      self.storage['path']['geom']=flpth+'geom/'+flnm
      self.storage['path']['time']=flpth+'time/'+flnm
      self.storage['path']['out']=flpth+'out/'+flnm
      self.storage['path']['CollFreqs']=flpth+'other/'+flnm +'/CollFreqs'
      self.storage['path']['cfdens']=flpth+'cfdens/'+flnm
      self.storage['path']['parallel']=flpth+'parallel/'+flnm
      self.storage['path']['fluxes_phi']=flpth+'fluxes/'+flnm
      self.storage['path']['fluxes_apar']=flpth+'fluxes/em/'+flnm
      self.storage['path']['fluxes_bpar']=flpth+'fluxes/bpar/'+flnm
      self.storage['path']['fluxes_lab_phi']=flpth+'fluxes_lab/'+flnm
      self.storage['path']['fluxes_lab_apar']=flpth+'fluxes_lab/em/'+flnm
      self.storage['path']['fluxes_lab_bpar']=flpth+'fluxes_lab/bpar/'+flnm
      self.storage['path']['fields']=flpth+'spectra_3D/fields/'+flnm+'/'
      self.storage['path']['moments']=flpth+'spectra_3D/moments/'+flnm+'/'
      self.storage['path']['j0_moments']=flpth+'spectra_3D/j0_moments/'+flnm+'/'
      self.storage['path']['j1_moments']=flpth+'spectra_3D/j1_moments/'+flnm+'/'


  def check_attributes(attr_list): #to be called only from the class
    """
    """

  def check_keys(dicts_keys_list): #to be called only from the class
    """
    """

  def load_file(flname,fltype,endianness=""): #to be called only from the class
    """ Generic routine to load GKW I/O files depending of their type 
    Inputs    
      flname      path+name of the file to open
      fltype      one among nml       -> fortran namelist, loaded with f90nml.read(f)
                            ftext     -> formated text (e.g. geom.dat), loaded with f.read
                            array     -> array of floats (e.g. fluxes.dat), loaded with np.loadtxt(f)
                            binary    -> array of floats in binary format (double precision)
                            binary_sp -> array of floats in binary format (single precision)
      endianness  can be '>'(big) or '<' (little). If empty attempts automatic detection.
                  only used for binary format

    Outputs
      Depending on fltype: nml   ->  dict with the content of the namelist
                           ftext ->  string
                           array ->  numpy array 
                           binary -> numpy array
    """
    assert (fltype in ('nml','ftext','array','binary','binary_sp')), \
           "Accepted values for 'fltype' are 'nml', 'ftext', 'array', 'binary' or 'binary_sp'"
    if fltype=='binary':
      prec='d'
    elif fltype=='binary_sp':
      prec='f'
    try:
      with open(flname,'r') as f:
        if fltype=='nml': # can be replaced by match/case, once python 3.10 is out
          return f90nml.read(f)
        if fltype=='ftext':
          return f.read()
        if fltype=='array':
          return np.loadtxt(f,dtype=float)
        if fltype=='binary' or fltype=='binary_sp':
          if endianness!='<' and endianness!='>': # if not known, attempts automatic detection
            endianness='<' # test little endian first
            dum=np.fromfile(f,dtype=endianness+prec)
            f.seek(0) #rewind file 
            if np.any(dum>1e50):
              endianness='>'
          return np.fromfile(f,dtype=endianness+prec), endianness
    except FileNotFoundError as err:
      print("Warning: file "+err.filename+" not found.")
      if fltype=='nml': # can be replaced by match/case, once python 3.10 is out
        return dict()
      if fltype=='ftext':
        return ""
      if fltype=='array':
        return np.array([])
      if fltype=='binary' or fltype=='binary_sp':
        return np.array([]),endianness


  def get_inputs(self):
    """ Load GKW input files (both input.dat and input_out.dat) in the following dicts:
          gkwrun.input
          gkwrun.input_out
    """
    self.input=gkwrun.load_file(self.storage['path']['input'],'nml')
    self.input_out=gkwrun.load_file(self.storage['path']['input_out'],'nml')
    self.nsp=self.input_out['gridsize']['number_of_species']
    self.ns=self.input_out['gridsize']['n_s_grid']
    self.nkx=self.input_out['gridsize']['nx']
    self.nky=self.input_out['gridsize']['nmod']
   

  def get_geom(self):
    """ Load GKW geom file in the following dict:
          gkwrun.geom
    """
    str=gkwrun.load_file(self.storage['path']['geom'],'ftext')
    p=re.compile('^(\S+)\n',flags = re.MULTILINE) # finds text starting from the beginning of the line
    m=p.split(str)
    self.geom={}
    for i in range(1,len(m),2):
      try:
        self.geom[m[i]]=int(m[i+1])
      except ValueError:
        try:
          self.geom[m[i]]=float(m[i+1])
        except ValueError:
          self.geom[m[i]]=np.fromstring(m[i+1],dtype=float,sep=' ')
    # compute krnorm manually if missing from geom.dat file
    if 'krnorm' not in self.geom:
      print("'krnorm' not available in geom.dat, computed from 'g_eps_eps'")
      f=interp1d(self.geom['s_grid'],np.sqrt(self.geom['g_eps_eps']),kind='cubic')
      self.geom['krnorm']=f(0)


  def get_eigenvalues(self):
    """ Load mode growth rate and real frequency in
          gkwrun.eigenvalues['gamma']
          gkwrun.eigenvalues['freq']
    """
    dum=gkwrun.load_file(self.storage['path']['time'],'array')
    self.eigenvalues={}
    # this works for initial value runs, need to treat eigenvalue runs differently
    if dum.size!=0 and dum.ndim==2:
      self.eigenvalues['gamma']=dum[:,1]
      self.eigenvalues['freq']=dum[:,2]
    else:
      self.eigenvalues['gamma']=np.array([])
      self.eigenvalues['freq']=np.array([])

 
  def get_grids(self):
    """ Load the various grids: time, s, kthrho, krrho, keps, kzeta, file_count in the dict
          gkwrun.grids
    """
    self.grids={}
    # time grid
    # this works for initial value runs, need to treat eigenvalue runs differently?
    # if not, each "time slice" corresponds to a different eigenvalue, maybe ok but not very straigtforward
    dum=gkwrun.load_file(self.storage['path']['time'],'array')
    if dum.ndim<2:
      self.grids['time']=dum
    else:
      self.grids['time']=dum[:,0]
    self.nt=dum.shape[0]

    # s grid (from geom outputs)
    self.grids['s']=self.geom['s_grid']

    # kx and ky grids
    dum=gkwrun.load_file(self.storage['path']['keps'],'array')
    if dum.ndim==2:
      self.grids['keps']=dum[0,:]
    else:
      self.grids['keps']=dum
    self.grids['krrho']=np.array(self.grids['keps']*self.geom['krnorm'])

    dum=gkwrun.load_file(self.storage['path']['kthrho'],'array')
    if dum.ndim==2:
      self.grids['kthrho']=dum[:,0]
    else:
      self.grids['kthrho']=dum
    self.grids['kzeta']=np.array(self.grids['kthrho']/self.geom['kthnorm'])

    # file counter for kykxs diagnostics
    dum=gkwrun.load_file(self.storage['path']['file_count'],'array')
    if dum.size==0 and not self.input_out['CONTROL']['non_linear']: # missing file_count, attempts to built it for linear runs
      self.grids['file_counter']=np.arange(1,len(self.grids['time'])+1) 
      print("Warning, file_count missing: build it from the time grid (may fail)")
    else:
      self.grids['file_counter']=dum.astype('int')


  def get_info(self):
    """ Load the content of the GKW 'out.dat' file as a string and parse some if its content
          gkwrun.info
    """
    self.info={}
    self.info['out']=gkwrun.load_file(self.storage['path']['out'],'ftext')


  def get_collisions(self):
    """ Load the collision frequency for all species (i.e. gammab in GKW)
        nu_i_j [vth_i/Rref] with i lines, j columns in the array
          gkwrun.collisions
    """
    dum=gkwrun.load_file(self.storage['path']['CollFreqs'],'array')
    if dum.size!=0:
      dum=dum[:,2]
      self.collisions=np.reshape(dum,(self.nsp,self.nsp))
    else:
      self.collisions=np.array([])


  def get_cfdens(self):
    """ Load the poloidally varying background density and density gradient on the s grid
          gkwrun.cfdens['n_pol']  [nspecies] with n_pol=n(s)/n(s=0)
          gkwrun.cfdens['rln_pol']    [nspecies]
    """
    self.cfdens={}
    dum=gkwrun.load_file(self.storage['path']['cfdens'],'array')    
    if dum.size!=0:
      self.cfdens['rln_pol']=dum[:,1:1+self.nsp]
      self.cfdens['n_pol']=dum[:,1+self.nsp:-1]

  def get_fluxes(self):
    """ Load the fluxes (GKW units) in 
          gkwrun.fluxes['phi']   [ntime x nspecies]
          gkwrun.fluxes['apar']  [ntime x nspecies]
          gkwrun.fluxes['bpar']  [ntime x nspecies]
          gkwrun.fluxes_lab['phi']   [ntime x nspecies]
          gkwrun.fluxes_lab['apar']  [ntime x nspecies]
          gkwrun.fluxes_lab['bpar']  [ntime x nspecies]
    """
    fl_list1=['fluxes_phi','fluxes_apar','fluxes_bpar']
    fl_list2=['fluxes_lab_phi','fluxes_lab_apar','fluxes_lab_bpar']
    key_list=['phi','apar','bpar']
    self.fluxes={}
    self.fluxes_lab={}
    for fl1,fl2,k in zip(fl_list1,fl_list2,key_list):
      dum=gkwrun.load_file(self.storage['path'][fl1],'array')
      self.fluxes[k]=dum
      if dum.size!=0:
        dum=np.reshape(dum,(self.nt,self.nsp*3))
        self.fluxes[k]={}
        self.fluxes[k]['particles']=dum[:,0::3]
        self.fluxes[k]['heat']=dum[:,1::3]
        self.fluxes[k]['momentum_parallel']=dum[:,2::3]
      dum=gkwrun.load_file(self.storage['path'][fl2],'array')
      self.fluxes_lab[k]=dum
      if dum.size!=0:
        dum=np.reshape(dum,(self.nt,self.nsp*3))
        self.fluxes_lab[k]={}
        self.fluxes_lab[k]['particles']=dum[:,0::3]
        self.fluxes_lab[k]['heat']=dum[:,1::3]
        self.fluxes_lab[k]['momentum_parallel']=dum[:,2::3]


  def get_parallel(self):
    """ Load the parallel structure of fields and moments from parallel.dat in
          gkwrun.parallel['phi']      [ns]
          gkwrun.parallel['apar']     [ns]
          gkwrun.parallel['bpar']     [ns]
          gkwrun.parallel['dens']        [ns x nspecies]
          gkwrun.parallel['vpar']     [ns x nspecies]
          gkwrun.parallel['Tpar']     [ns x nspecies]
          gkwrun.parallel['Tperp']    [ns x nspecies]
    """
    dum=gkwrun.load_file(self.storage['path']['parallel'],'array')
    self.parallel={}
    if dum.size!=0:
      self.parallel['phi']=dum[0:self.ns,1]+1j*dum[0:self.ns,2]
      self.parallel['apar']=dum[0:self.ns,3]+1j*dum[0:self.ns,4]
      self.parallel['bpar']=dum[0:self.ns,13]+1j*dum[0:self.ns,14]
      self.parallel['dens']=dum[0:self.ns,5]+1j*dum[0:self.ns,6]
      self.parallel['vpar']=dum[0:self.ns,11]+1j*dum[0:self.ns,12]
      self.parallel['Tpar']=dum[0:self.ns,7]+1j*dum[0:self.ns,8]
      self.parallel['Tperp']=dum[0:self.ns,9]+1j*dum[0:self.ns,10]

  def get_kykxs(self,to_load,time_interval=np.array([]),prec='sp',kx_interval=np.array([])):
    """ Load the kykxs fields/moments specified in the to_load list for a given time interval and kx interval
        Available fields/moments are:
          phi, apar, bpar
          dens, vpar, Tpar, Tperp, Qpar, M12, M24
          dens_J0, vpar_J0, Tpar_J0, Tperp_J0, Qpar_J0, M12_J0, M24_J0
          Tperp_J1 

        Inputs:
          to_load          list of fields/moments to load (string or set of strings, e.g. {'phi','apar'})
          time_interval    time interval for which the files will be loaded [t_start t_end] (numpy.ndarray)
                           if empty, loads all available files
                           not used if xy_estep=F (fields/moments output at last time step only)
          prec             precision used for the binary files, available 'dp' (double) and 'sp' (single)
          kx_interval      not implemented yet
       
        Outputs
         gkwrun.kykxs['time_interval']   effective time interval for which the data is loaded
                     ['It']              corresponding indexes in the time grid
                     ['phi']             loaded data, complex array of size nky,nkx,ns,nt,(nsp) 
                     ['apar']
                     etc.             
    """
    kykxs_list={'phi': ('Phi','fields'), 'apar': ('Apa','fields'), 'bpar': ('Bpa','fields'), 
                'dens': ('dens','moments'), 'vpar': ('vpar','moments'), 'Tpar': ('Tpar','moments'),
                'Tperp': ('Tperp','moments'), 'Qpar': ('Qpar','moments'), 'M12': ('M12','moments'),
                'M24': ('M24','moments'),
                'dens_J0': ('dens_ga','j0_moments'), 'vpar_J0': ('vpar_ga','j0_moments'),
                'Tpar_J0': ('Tpar_ga','j0_moments'), 'Tperp_J0': ('Tperp_ga','j0_moments'),
                'Qpar_J0': ('Qpar_ga','j0_moments'), 'M12_J0': ('M12_ga','j0_moments'), 
                'M24_J0': ('M24_ga','j0_moments'), 'Tperp_J1': ('Tperp_J1','j1_moments')
               }
    # checks
    if type(to_load)==str:
      to_load={to_load}
    time_interval=np.asarray(time_interval)
    kx_interval=np.asarray(kx_interval)
    # initialise kykxs dict if needed
    if not hasattr(self,'kykxs'):
      self.kykxs={}
    # find time indexes (check if output was requested at each time step)
    if self.input_out['diagnostic']['xy_estep']:
      if time_interval.size==0:
        dum=self.grids['time'][[0,-1]]
        It=range(self.nt)
      else:
        I1=np.abs(self.grids['time']-time_interval[0]).argmin()
        I2=np.abs(self.grids['time']-time_interval[1]).argmin()
        dum=self.grids['time'][[I1, I2]]
        It=range(I1,min(I2+1,self.nt))
    else:
      dum=self.grids['time'][[-1,-1]]
      It=range(self.nt-1,self.nt)
      if time_interval.size!=0:
        print("Warning, only last time step available for kykxs diagnostic")
    # store time interval and local time grid
    if ('time_interval' in self.kykxs) and np.any(self.kykxs['time_interval']!=dum):
      print("Warning: 'time_interval' already defined for kykxs fields/moments and different from requested. No data loaded.")
      return
    else:
      self.kykxs['time_interval']=dum
    self.kykxs['It']=It
    # find kx indexes (to implement later if needed)
    # either by specifying the kx interval or with a switch so that kx_max=3*ky_max
    Ikx=[]
    self.kykxs['Ikx']=Ikx

    endianness=''
    for kk in to_load:
      if kk in kykxs_list:
        flpth=self.storage['path'][kykxs_list[kk][1]]
        if kykxs_list[kk][1]=='fields':  # fields
          self.kykxs[kk]=np.full((self.nky,self.nkx,self.ns,len(It)),np.nan+1j*np.nan)
          for ii in range(len(It)):
            flroot=kykxs_list[kk][0]+'_kykxs'+str(self.grids['file_counter'][It[ii]]).zfill(8)
            if prec=='sp': # file stored in single precision
              dum_r,endianness=gkwrun.load_file(flpth+flroot+'_real_sp','binary_sp',endianness)
              dum_i,endianness=gkwrun.load_file(flpth+flroot+'_imag_sp','binary_sp',endianness)
            else:          # assume double precision
              dum_r,endianness=gkwrun.load_file(flpth+flroot+'_real','binary',endianness)
              dum_i,endianness=gkwrun.load_file(flpth+flroot+'_imag','binary',endianness)
            if dum_r.size==self.nky*self.nkx*self.ns and dum_i.size==self.nky*self.nkx*self.ns:
              self.kykxs[kk][:,:,:,ii]=np.reshape(dum_r+1j*dum_i,(self.nky,self.nkx,self.ns),order='F')
            elif dum_r.size>0 and dum_i.size>0:
              print("Warning: wrong file size for "+flroot+". No data loaded")
        else:                            # moments
          self.kykxs[kk]=np.full((self.nky,self.nkx,self.ns,len(It),self.nsp),np.nan+1j*np.nan)
          for ii in range(len(It)):
            for jj in range(0,self.nsp):
              flroot=kykxs_list[kk][0]+'_kykxs'+str(jj+1).zfill(2)+'_'+str(self.grids['file_counter'][It[ii]]).zfill(6)
              if prec=='sp': # file stored in single precision
                dum_r,endianness=gkwrun.load_file(flpth+flroot+'_real_sp','binary_sp',endianness)
                dum_i,endianness=gkwrun.load_file(flpth+flroot+'_imag_sp','binary_sp',endianness)
              else:          # assume double precision
                dum_r,endianness=gkwrun.load_file(flpth+flroot+'_real','binary',endianness)
                dum_i,endianness=gkwrun.load_file(flpth+flroot+'_imag','binary',endianness)
              if dum_r.size==self.nky*self.nkx*self.ns and dum_i.size==self.nky*self.nkx*self.ns:
                self.kykxs[kk][:,:,:,ii,jj]=np.reshape(dum_r+1j*dum_i,(self.nky,self.nkx,self.ns),order='F')
              elif dum_r.size>0 and dum_i.size>0:
                print("Warning: wrong file size for "+flroot+". No data loaded")
        if np.all(np.isnan(self.kykxs[kk])):
          self.kykxs[kk]=np.array([])
      else:
        print("Warning: Field/moments "+kk+" not defined. Skipped.")



  def gkw2ids(self,time_interval=np.array([]),Nsh=10):
    """ Convert GKW inputs and outputs to an IMAS 'gyrokinetics' IDS
        Inputs:
          time interval used for the outputs (fluxes + eigenmodes for NL, eigenmodes only for L)
            If empty:
              NL run: all available times
              linear run: last time step only
            If dim 2:
              NL run: time selection (time averaged fluxes + time dependent eigenmodes)
              linear run, initial value: time selection (time dependent eigenmodes only) 
            Not used for linear eigenvalue runs (or use it to select the eigenvalue?)         
          Nsh  number of moments for flux surface Fourier parametrisation (default Nsh=10)
    """
    self.ids={}
    time_interval=np.asarray(time_interval)

    ## constants (from IMAS 3.23.2)
    me=9.1093837015e-31 # electron mass [kg] 
    mD=3.3435837724e-27 # deuterium mass [kg]
    eV=3.3435837724e-19 # elementary charge [C] / electron volt [J]

    ## code parameters
    self.ids['code']={}
    self.ids['code']['name']="GKW"
    self.ids['code']['commit']=""
    self.ids['code']['parameters']="" # to deal with later, can use dict2xml and back conversion with untangle.parse())
                                      # fill the general code table + the specific code table in eigenmode for linear runs only
    ## model
    self.ids['model']={}
    if self.input_out['CONTROL']['method']=='EXP':
      self.ids['model']['initial_value_run']=True
    elif self.input_out['CONTROL']['method']=='EIV':
      self.ids['model']['initial_value_run']=False
    else:
      print("Only explicit time integration or eigenvalue runs are dealt with")
      return

    if self.input_out['CONTROL']['non_linear']:
      self.ids['model']['non_linear_run']=True
      if time_interval.size==0: # take all available times 
        self.ids['model']['time_interval_norm']=self.grids['time'][[0,-1]]
        It_fluxes=range(self.nt)
        It_moments=range(self.nt)
      else: 
        I1=np.abs(self.grids['time']-time_interval[0]).argmin()
        I2=np.abs(self.grids['time']-time_interval[1]).argmin()
        self.ids['model']['time_interval_norm']=self.grids['time'][[I1, I2]]
        It_fluxes=range(I1,min(I2+1,self.nt))
        It_moments=range(I1,min(I2+1,self.nt))
    else: 
      self.ids['model']['non_linear_run']=False 
      self.ids['model']['time_interval_norm']=self.grids['time'][[-1,-1]] # fluxes given at last time step only
      It_fluxes=range(self.nt-1,self.nt)
      It_moments=range(self.nt-1,self.nt)
      if self.ids['model']['initial_value_run']==True & time_interval.size!=0: # time interval for moments
        I1=np.abs(self.grids['time']-time_interval[0]).argmin()
        I2=np.abs(self.grids['time']-time_interval[1]).argmin()
        It_moments=range(I1,min(I2+1,self.nt))

    self.ids['model']['include_a_field_parallel']=self.input_out['CONTROL']['nlapar']
    self.ids['model']['include_b_field_parallel']=self.input_out['CONTROL']['nlbpar']

    if self.input_out['ROTATION']['cf_drift'] or self.input_out['ROTATION']['cf_trap']:
      if (self.input_out['ROTATION']['cf_drift'] and self.input_out['ROTATION']['cf_trap']
          and self.input_out['ROTATION']['cf_upphi'] and self.input_out['ROTATION']['cf_upsrc']):
        self.ids['model']['include_centrifugal_effects']=True
      else:
        print("Inconsistent CF switches: all or none need to be ON for the GK IDS")
        return
    else:
      self.ids['model']['include_centrifugal_effects']=False

    self.ids['model']['collisions_pitch_only'] = (self.input_out['COLLISIONS']['pitch_angle'] 
      and not self.input_out['COLLISIONS']['en_scatter'] and not self.input_out['COLLISIONS']['friction_coll'])
    self.ids['model']['collisions_momentum_conservation'] = self.input_out['COLLISIONS']['mom_conservation']
    self.ids['model']['collisions_energy_conservation'] = self.input_out['COLLISIONS']['ene_conservation']
    self.ids['model']['collisions_finite_larmor_radius']=False


    ## flux_surface
    self.ids['flux_surface']={}
    assert (self.input_out['GEOM']['geom_type'] in ('miller','fourier','chease')), \
           " 'miller', 'fourier', 'chease' are the only supported equilibrium descriptions for the GK IDS"

    if self.input_out['GEOM']['geom_type']=='miller' or self.input_out['GEOM']['geom_type']=='fourier':

      if self.input_out['GEOM']['gradp_type']=='beta_prime':
        if self.input_out['SPCGENERAL']['betaprime_type']=='ref':
          self.ids['model']['include_full_curvature_drift']=True
          betapr_gkw=self.input_out['SPCGENERAL']['betaprime_ref']
        else:
          print("With 'miller' or 'fourier' parametrisation, betaprime_type='ref' is the only supported option (could be extended)")
          # could be extended to support betaprime_type='sp' if needed 
      else:
        print("With 'miller' or 'fourier' parametrisation, gradp_type=beta_prime is the only supported option (could be extended)")
        return

      if self.input_out['GEOM']['geom_type']=='fourier':
        # First compute R/Rref_GKW and Z/Rref_GKW from the Fourier parametrisation
        Nsh=self.input_out['GEOM']['n_shape']
        c=self.input_out['GEOM']['c'][0:Nsh]
        s=self.input_out['GEOM']['s'][0:Nsh]
        dcdr=self.input_out['GEOM']['c_prime'][0:Nsh]
        dsdr=self.input_out['GEOM']['s_prime'][0:Nsh]
        rdum,Rdum,Zdum=FS.fourier2rz(c,s,dcdr,dsdr,'gkw')
        # and compute back the fourier parametrisation (in case the initial parametrisation was not with the right center R0,Z0)
        cc,ss,dccdr,dssdr,R0,Z0,R_out,Z_out,err_out=FS.rz2fourier(Rdum,Zdum,rdum[1],'imas',Nsh,doplots=False)
 
      if self.input_out['GEOM']['geom_type']=='miller':
        # First compute R/Rref_GKW and Z/Rref_GKW from the Miller parametrisation
        r0=self.input_out['GEOM']['eps']
        Rmil=1
        Zmil=self.input_out['GEOM']['zmil']
        dRmildr=self.input_out['GEOM']['drmil']
        dZmildr=self.input_out['GEOM']['dzmil']
        k=self.input_out['GEOM']['kappa']
        d=self.input_out['GEOM']['delta']
        z=self.input_out['GEOM']['square']
        sk=self.input_out['GEOM']['skappa']
        sd=self.input_out['GEOM']['sdelta']
        sz=self.input_out['GEOM']['ssquare']
        rdum,Rdum,Zdum=FS.miller2rz(r0,Rmil,Zmil,k,d,z,dRmildr,dZmildr,sk,sd,sz)
        # and compute the Fourier parametrisation for the IDS
        cc,ss,dccdr,dssdr,R0,Z0,R_out,Z_out,err_out=FS.rz2fourier(Rdum,Zdum,rdum[1],'imas',Nsh,doplots=False)

      R_rat=1/R0   # Rref_GKW/Rref_GKDD
      B_rat=1   # Bref_GKW/Bref_GKDD

    elif self.input_out['GEOM']['geom_type']=='chease':
      betapr_gkw=self.geom['betaprime_eq']
      if self.input_out['SPCGENERAL']['betaprime_type']=='eq':
        self.ids['model']['include_full_curvature_drift']=True
      elif (self.input_out['SPCGENERAL']['betaprime_type']=='sp' and self.input_out['SPCGENERAL']['beta']==0) or \
           (self.input_out['SPCGENERAL']['betaprime_type']=='ref' and self.input_out['SPCGENERAL']['betaprime_ref']==0):
        self.ids['model']['include_full_curvature_drift']=False
      else:
        print("Warning: No way to check that the betaprime used in the curvature drift is consistent with the magnetic equilibrium") 
        return
      # need to read hamada file, and get R,Z from it
      print("Treatment of chease equilibrium not yet implemented")
      return

    sb_gkw=self.input_out['GEOM']['signb']
    sj_gkw=self.input_out['GEOM']['signj']
    self.ids['flux_surface']['r_minor_norm']=self.geom['eps']*R_rat
    self.ids['flux_surface']['q']=(-sb_gkw)*(-sj_gkw)*self.geom['q']
    self.ids['flux_surface']['magnetic_shear_r_minor']=self.geom['shat']
    self.ids['flux_surface']['pressure_gradient_norm']=-betapr_gkw*B_rat**2/R_rat
    self.ids['flux_surface']['b_field_tor_sign']=-sb_gkw
    self.ids['flux_surface']['ip_sign']=-sj_gkw
    self.ids['flux_surface']['shape_coefficients_c']=cc*R_rat
    self.ids['flux_surface']['shape_coefficients_s']=ss*R_rat
    self.ids['flux_surface']['dc_dr_minor_norm']=dccdr
    self.ids['flux_surface']['ds_dr_minor_norm']=dssdr

    # compute theta_imas=f(s_gkw)
    R_dum = self.geom['R']
    Z_dum = self.geom['Z']
    th = np.arctan2(-(Z_dum-Z0),R_dum-R0)
    dum=np.cumsum(np.concatenate(([False],np.diff(th)>np.pi)))
    th_imas_of_s_grid=th-2*np.pi*dum+2*np.pi*dum[int(self.ns/2)]

    ## compute n_cf=n_s(theta=0)/n_s(s=0) and rln_th0=rln_s(theta=0)
    n_cf=np.ones(self.nsp)
    rln_th0=np.array([x['rln'] for x in self.input_out['SPECIES']])
    if self.ids['model']['include_centrifugal_effects']==True:
      assert self.input_out['GEOM']['r0_loc']=='LFS', \
             "On axis specification of density and density gradient not implemented (R0_LOC='LFS' required)."
      for ii in range(self.nsp):
        f=interp1d(th_imas_of_s_grid,self.cfdens['n_pol'][:,ii],kind='cubic')
        n_cf[ii]=f(0)
        f=interp1d(th_imas_of_s_grid,self.cfdens['rln_pol'][:,ii],kind='cubic')
        rln_th0[ii]=f(0)

    ## find electron species
    assert not(self.input_out['SPCGENERAL']['adiabatic_electrons']), \
           "Adiabatic electrons not allowed in a GK IDS"
    dum=[sp['z']<0 for sp in self.input_out['SPECIES']]
    assert (dum.count(True)==1), \
           "Runs with multiple or no electron species not handled yet"
    Iele=dum.index(True) # index of the electron species

    ## compute ratio of reference quantities for GKW to IMAS conversion
    # reference charge ratio, qref_GKW/qref_IMAS    
    q_rat = -1/self.input_out['SPECIES'][Iele]['z'] 
    # reference mass ratio, mref_GKW/mref_IMAS 
    m_rat = me/mD*1/self.input_out['SPECIES'][Iele]['mass']
    # reference temperature ratio, Tref_GKW/Tref_IMAS
    T_rat = 1/self.input_out['SPECIES'][Iele]['temp']
    # reference density ratio, nref_GKW/nref_IMAS
    n_rat = 1/self.input_out['SPECIES'][Iele]['dens']/n_cf[Iele]
    # additional reference ratio
    vth_rat = np.sqrt(T_rat/m_rat)
    rho_rat = q_rat*B_rat/(m_rat*vth_rat)

    ## species
    self.ids['species']=[]
    for ii in range(self.nsp):
      dum={}
      dum['charge_norm']=self.input_out['SPECIES'][ii]['z']*q_rat
      dum['mass_norm']=self.input_out['SPECIES'][ii]['mass']*m_rat
      dum['density_norm']=self.input_out['SPECIES'][ii]['dens']*n_rat*n_cf[ii]
      dum['density_log_gradient_norm']=rln_th0[ii]/R_rat
      dum['temperature_norm']=self.input_out['SPECIES'][ii]['temp']*T_rat
      dum['temperature_log_gradient_norm']=self.input_out['SPECIES'][ii]['rlt']/R_rat
      dum['velocity_tor_gradient_norm']=sb_gkw*self.input_out['SPECIES'][ii]['uprim']*vth_rat/R_rat**2
      self.ids['species'].append(dum)


    ## species_all
    self.ids['species_all']={}
    self.ids['species_all']['debye_length_reference']=0.0 # not taken into account in GKW

    self.ids['species_all']['beta_reference']=0.0
    if self.input_out['CONTROL']['nlapar']==True | self.input_out['CONTROL']['nlbpar']==True: # beta value
      if self.input_out['SPCGENERAL']['beta_type']=='ref':
        self.ids['species_all']['beta_reference']=self.input_out['SPCGENERAL']['beta_ref']*B_rat**2/(n_rat*T_rat)
      elif self.input_out['SPCGENERAL']['beta_type']=='eq':
        dum=np.sum([x['dens']*x['temp'] for x in self.input_out['SPECIES']])
        self.ids['species_all']['beta_reference']=self.geom['beta_eq']/dum*B_rat**2/(n_rat*T_rat)
      else:
        print("Unknown beta_type option in SPCGENERAL namelist")
        return
    self.ids['species_all']['velocity_tor_norm']=sb_gkw*self.input_out['ROTATION']['vcor']*vth_rat/R_rat

    if self.ids['model']['initial_value_run']==False:
      self.ids['species_all']['shearing_rate_norm']=0.0
    else:
      if self.input_out['rotation']['shear_profile']=='none':
        self.ids['species_all']['shearing_rate_norm']=0.0
      elif self.input_out['rotation']['shear_profile']=='wavevector_remap':
        if self.input_out['rotation']['t_shear_begin']<self.grids['time'][It_moments[0]]:
          if self.input_out['rotation']['toroidal_shear']=='none':
            self.ids['species_all']['shearing_rate_norm']=-self.input_out['rotation']['shear_rate']*B_rat*vth_rat/R_rat
          elif (self.input_out['rotation']['toroidal_shear']=='use_uprim' or 
               self.input_out['rotation']['toroidal_shear']=='add_uprim'):
            # extract actual shear_rate value from out file
            dum=np.float(re.search('shear\_rate used \:(.*)',self.info['out']).group(1))
            self.ids['species_all']['shearing_rate_norm']=-dum*B_rat*vth_rat/R_rat
          else:
             print("Unknown or unsupported toroidal_shear option in ROTATION namelist")
        elif self.input_out['rotation']['t_shear_begin']>self.grids['time'][It_moments[-1]]:
          self.ids['species_all']['shearing_rate_norm']=0.0 # shearing applied after the specified interval
        else:
          print("Varying shearing rate within the specified time interval not allowed")
          return
      else:
        print("Unknown or unsupported shear_profile option in ROTATION namelist") 
        return

    ## collisions
    self.ids['collisions']={}
    self.ids['collisions']['collisionality_norm']=np.full((self.nsp,self.nsp),0.0)
    if self.collisions.size>0:
      for ii in range(self.nsp):
        self.ids['collisions']['collisionality_norm'][ii,:]=(self.collisions[ii,:]
                            *vth_rat/R_rat
                            *np.sqrt(self.input_out['species'][ii]['temp']/self.input_out['species'][ii]['mass']))

    ## wavevector
    self.ids['wavevector']=[]
    f=interp1d(th_imas_of_s_grid,self.geom['g_eps_eps'],kind='cubic')
    gxx_th0=f(0)
    f=interp1d(th_imas_of_s_grid,self.geom['g_zeta_zeta'],kind='cubic')
    gyy_th0=f(0)
    if self.nkx*self.nky==1:
        dum={}
        dum['radial_component_norm']=-sb_gkw*sj_gkw*self.grids['keps']*gxx_th0/rho_rat
        dum['binormal_component_norm']=self.grids['kzeta']*gyy_th0/rho_rat
        self.ids['wavevector'].append(dum)
    else:
      for ii_kx in range(self.nkx):
        for ii_ky in range(self.nky):
          dum={}
          dum['radial_component_norm']=-sb_gkw*sj_gkw*self.grids['keps'][ii_kx]*gxx_th0/rho_rat
          dum['binormal_component_norm']=self.grids['kzeta'][ii_ky]*gyy_th0/rho_rat
          self.ids['wavevector'].append(dum)


    ## eigenmode
    Ith_sorted=np.argsort(th_imas_of_s_grid)
    th_sorted=th_imas_of_s_grid[Ith_sorted]

  
  def write_ids(self,file_format):
    """ Write the 'gyrokinetics' IDS on file (JSON or HDF5)
    """

  def __init__(self): # initialise attributes here
    # possibly load all inputs here 
    return

  def __str__(self): # informal string representation of the object, for the user
    return ''

