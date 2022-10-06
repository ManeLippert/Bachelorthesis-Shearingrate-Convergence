import h5py, os, subprocess, pandas as pd, numpy as np

filename = [f for f in os.listdir() if f.endswith('.h5')]

dataset = {'evaluation':['second_derivative_phi', 'zonalflow_potential', 'shearing_rate', 'shearing_rate_maximum']}

for dirpath, dirnames, filenames in os.walk("."):
    for filename in [f for f in filenames if f.endswith(".h5")]:
        
        dirpath = dirpath.replace('./', '')
        
        file = dirpath + '/' + filename        

        print(file)
        
        #subprocess.run(['rm', '-rf', file])
        
        #f = h5py.File(file, 'r')

        #for node in dataset:
        #    for data in dataset[node]:
        #        df = pd.DataFrame(f[node + '/' + data][()].T)
        #        df.to_csv(dirpath + '/' + 'data-evaluation-' + data + '.csv', index=False, header=False)

        #time = f['diagnostic/diagnos_growth_freq/time'][()]
        #nt = time.shape[1]     

        #df = pd.DataFrame(time.T)
        #df.to_csv(dirpath + '/' + 'data-diagnostic-diagnos_growth_freq-' + 'time.csv', index=False, header=False)

        #eflux = np.reshape(f['diagnostic/diagnos_fluxes/eflux_species01'][()].T,(nt,2))[:,0]

        #df = pd.DataFrame(eflux)
        #df.to_csv(dirpath + '/' + 'data-diagnostic-diagnos_fluxes-' + 'eflux_species01.csv', index=False, header=False)