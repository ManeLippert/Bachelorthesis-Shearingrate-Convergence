import sys
import subprocess
import pkg_resources

required  = {'numpy', 'pandas', 'h5py'} 
installed = {pkg.key for pkg in pkg_resources.working_set}
missing   = required - installed

if missing:
    # implement pip as a subprocess:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', *missing], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
import numpy
import pandas
import h5py

print("Modules imported")