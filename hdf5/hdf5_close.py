#!/usr/bin/env python3

import h5py
import os

filename = [f for f in os.listdir() if f.endswith('.h5')]

for file in filename:

	try:
		f = h5py.File(file, 'r')
		f.close()
		print(file + ' successfully closed!')
	except OSError:
		print('! ' + file + ' might be broken !')