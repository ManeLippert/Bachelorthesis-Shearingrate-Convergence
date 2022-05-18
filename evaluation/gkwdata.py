import h5py

def gkw_data(f):

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
            elif type(f[keys0 + '/' + keys1]) == h5py._hl.group.Group:
                for keys2 in list(f[keys0 + '/' + keys1].keys()):
                    data[i][j].append(keys0 + '/' + keys1 + '/' + keys2)
            j += 1
        i += 1
        
    # !Important! close h5 file after usage 
    f.close()
    
    return data