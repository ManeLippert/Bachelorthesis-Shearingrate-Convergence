#!/usr/bin/env python3

import h5py

try:
	f = h5py.File('gkwdata.h5', 'r')
	f.close()
	print('File successfully closed!')
except OsError:
	print('! File might be broken !')

