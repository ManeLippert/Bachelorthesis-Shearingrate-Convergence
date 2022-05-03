cd /home/yann/codes/gkw/python/yc
import main 
import importlib #to reload modules by doing importlib.reload(main)

# create an instance of the run class
tt=main.gkwrun()

# specify the type of file storage and get the corresponding paths
gkw_home='/home/yann/runs1/gkw/'
proj='PFS_Samuele'
flnm='Waltz_disp_1_1_1'  # linear run

gkw_home='/home/yann/runs1/gkw/'
proj='TCVREV_3'
flnm='p2LINhigh_2_2'    # linear run

gkw_home='/home/yann/runs1/gkw/'
proj='JET_KBMbis'
flnm='j15_NLfiW_rlte_exb_1'     # non-linear run


gkw_home='/home/yann/runs1/gkw/'
proj='GKDB_TEST'
flnm='ref_fourier'     # linear run, fourier

tt.set_storage(gkw_home+proj,'multiple_folders',flnm)

# load input files (input.dat and input_out.dat)
tt.get_inputs()

# load geom
tt.get_geom()

# load grids 
tt.get_grids()

# load output file
tt.get_info()

# load eigenvalues
tt.get_eigenvalues()

# load collisions
tt.get_collisions()

# load background density
tt.get_cfdens()

# load fluxes
tt.get_fluxes()

# load fields and moments from parallel.dat
tt.get_parallel()

# load kykxs fields and moments 
tt.get_kykxs('phi',[103.5,104.2])
tt.get_kykxs({'apar','Tperp','Tperp_J0'},[103.5,104.2])

# make a GK IDS from GKW inputs and outputs
tt.gkw2ids([])



# to do: put all IDS related routines in a separate class


# temporary, look how a JSON file is loaded in python
import json
with open('/home/yann/projects/gkdb/GENE_GKW/GKW-case-1-1.json') as f:
 ids=json.load(f)


# structure represented as a combination of dicts and lists of dicts
type(ids['wavevector']) # list
type(ids['wavevector'][0]) # list


    # ids properties (to fill later)
    self.ids['ids_properties']['comment']=""
    self.ids['ids_properties']['provider']=""
    self.ids['ids_properties']['creation_date']=""
