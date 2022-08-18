import subprocess
import os

backupLocation = '/scratch/bt712347/backup/'

user = os.getlogin()
path = os.path.dirname(os.path.abspath(__file__)).split(user + '/')[1]

backupPath = backupLocation + path

# Restore
subprocess.run(['rsync', '-a', backupPath + '/', ''])
