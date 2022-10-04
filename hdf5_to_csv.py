import h5py, os, pandas as pd

filename = [f for f in os.listdir() if f.endswith('.h5')]

dataset = {'evaluation':['second_derivative_phi', 'zonalflow_potential', 'shearing_rate', 'shearing_rate_maximum'],
           'diagnostics/diagnos_growth_freq': ['time'],
           'diagnostics/diagnos_fluxes': ['eflux_species01']}

for file in filename:
    f = h5py.File(file, 'r')
    
    for node in dataset:
        for data in dataset[node]:
            df = pd.DataFrame(f[node + '/' + data][()].T)
            df.to_csv(data + '.csv', index=False, header=False)