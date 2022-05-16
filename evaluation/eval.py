import h5py
import numpy

filename = "evaluation/data/S6_rtl6.3/gkwdata.h5"

file = h5py.File(filename)

keys = file.keys()
keys_list = list(file.keys())

group = file[keys_list[3]]
group_list = list(file[keys_list[3]])

dataset = group[group_list[0]]
dataset_list = list(group[group_list[0]])

data = dataset[dataset_list[0]]

print(list(data))

file.close()