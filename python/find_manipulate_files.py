import h5py, os, subprocess
import pandas as pd
import numpy as np

# VARIABLES ===========================================================================================================

pattern = "data.h5"
dataset = {"evaluation": ["derivative_stepsize"]}

# FUNKTIONS ===========================================================================================================

def file_loop(pattern):
    
    data = []
    
    for dirpath, dirnames, filenames in os.walk("."):
        for filename in [f for f in filenames if f.endswith(pattern)]:
            data.append(dirpath + "/" + filename)
    
    return data

def hdf5_extract(input_file, output_file, value):
    subprocess.run(["h5copy", "-p", "-i", input_file, "-o", output_file, "-s", value, "-d", value])

def hdf5_to_csv(dirpath, filename, dataset):

    dirpath = dirpath.replace("./", "")
    file = dirpath + "/" + filename

    f = h5py.File(file, "r")

    for node in dataset:
        for data in dataset[node]:
            # print(f[node + '/' + data][()].T)
            df = pd.DataFrame(f[node + "/" + data][()].T)
            df.to_csv(dirpath + "/" + data + ".csv", index=False, header=False)

    # time = f['diagnostic/diagnos_growth_freq/time'][()]
    # nt = time.shape[1]

    # df = pd.DataFrame(time.T)
    # df.to_csv(dirpath + '/' + 'data-diagnostic-diagnos_growth_freq-' + 'time.csv', index=False, header=False)

    # eflux = np.reshape(f['diagnostic/diagnos_fluxes/eflux_species01'][()].T,(nt,2))[:,0]

    # df = pd.DataFrame(eflux)
    # df.to_csv(dirpath + '/' + 'data-diagnostic-diagnos_fluxes-' + 'eflux_species01.csv', index=False, header=False)

def hdf5_combine(input_file, combined_file, groupname):
    
    def hdf5_write_data(f, data, groupname = 'added_data'):
    
        try:
            if type(data) == list:
                if type(groupname) == list:

                    for d, n in zip(data, groupname):
                        f.create_dataset(n, data = d)

                else:
                    print('Every data entry should have own groupname')        
            else:
                f.create_dataset(groupname, data = data)
        except (ValueError,OSError):
            pass

        f.close()
        
    f_input = h5py.File(input_file, "r+")
    f_combined = h5py.File(combined_file, "r+")
    
    data_input = f_input[groupname][()]
    
    hdf5_write_data(f_combined, data_input, groupname)
    
def delete_file(file):
    subprocess.run(['rm', '-rf', file])
    
def hdf5_close(file):
    
    try:
        f = h5py.File(file, 'r')
        f.close()
        print(file + ' successfully closed!')
    
    except OSError:
        print('! ' + file + ' might be broken !')

# MANIPULATION ================================================================================================================

data = file_loop('data.h5')
stepsize = file_loop('stepsize.h5')

data_data = []

for i in data:
    if "/data.h5" in i:
        data_data.append(i)

#print(data[0], stepsize[0])

#hdf5_extract(data[0], "./test.h5", "evaluation/shearing_rate_maximum")

for i in stepsize:
    delete_file(i)