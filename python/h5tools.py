import h5py
import os

# Returns a list with three levels, with all keys of gkw data
def get_keys(f):

    # Maximum level until datasets => 3
    data = []

    # First level
    i = 0
    for keys0 in list(f.keys()):
        data.append([])
        # Second level
        j = 0
        for keys1 in list(f[keys0].keys()):
            data[i].append([])
            # Third level
            if type(f[keys0 + '/' + keys1]) == h5py._hl.dataset.Dataset:
                data[i][j] = keys0 + '/' + keys1
                print('"' + keys0 + '/' + keys1 + '"')
            elif type(f[keys0 + '/' + keys1]) == h5py._hl.group.Group:
                for keys2 in list(f[keys0 + '/' + keys1].keys()):
                    data[i][j].append(keys0 + '/' + keys1 + '/' + keys2)
                    print('"' + keys0 + '/' + keys1 + '/' + keys2 + '"')
            j += 1
        i += 1
        
    return data

# Returns key if the end of the level is equal the vairable search
def find_key(f, search):

    # i, j, k are indexes for the function get_data_keys()
    
    # First level
    #i = 0
    for keys0 in list(f.keys()):
        # Second level
        #j = 0
        for keys1 in list(f[keys0].keys()):
            # Third level
            #k = 0
            if type(f[keys0 + '/' + keys1]) == h5py._hl.dataset.Dataset:
                key = keys0 + '/' + keys1
                if search == keys1:
                    return key
            elif type(f[keys0 + '/' + keys1]) == h5py._hl.group.Group:
                for keys2 in list(f[keys0 + '/' + keys1].keys()):
                    key = keys0 + '/' + keys1 + '/' + keys2
                    if search == keys2:
                        return key
                    #k += 1
            #j += 1
        #i += 1


# Prints all possible keys with the variable search in it
def find_keys(f, search):
    
    key_list = []

    # i, j, k are indexes for the function get_data_keys()
    
    # First level
    #i = 0
    for keys0 in list(f.keys()):
        # Second level
        #j = 0
        for keys1 in list(f[keys0].keys()):
            # Third level
            #k = 0
            if type(f[keys0 + '/' + keys1]) == h5py._hl.dataset.Dataset:
                key = keys0 + '/' + keys1
                if search in key:
                    key_list.append(key)
            elif type(f[keys0 + '/' + keys1]) == h5py._hl.group.Group:
                for keys2 in list(f[keys0 + '/' + keys1].keys()):
                    key = keys0 + '/' + keys1 + '/' + keys2
                    if search in key:
                        key_list.append(key)
                    #k += 1
            #j += 1
        #i += 1
     
    for i in key_list:
        print(i)
        
def hdf5_close():
    filename = [f for f in os.listdir() if f.endswith('.h5')]
    for file in filename:
        try:
            f = h5py.File(file, 'r')
            f.close()
            print(file + ' successfully closed!')
        except OSError:
            print('! ' + file + ' might be broken !')
            
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
    except ValueError:
        pass
        
    f.close()