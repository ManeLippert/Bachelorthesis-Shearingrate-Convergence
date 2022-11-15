import argparse

description_text = '''

=========================== DESCRIPTION ===========================

FEATURES:
o Creats jobscript file with defined content
o Restarts job until criteria is suffused
o Makes backup after each run before Restart
o Sends mail at the beginning, end and restart with status file 
  as attachment (mailx has to be defined)
o Creates status file with current status and progress bar

START SCRIPT IN BACKGROUND:
>>> nohup python3 -u slurm_monitor.py &> status.txt &

With arguments (example):
>>> nohup python3 -u slurm_monitor.py -n 30000 &> status.txt &

Output:
>>> [1] 10537
This will write every output in the status file 'status.txt'

LIST PROCESS:
>>> ps ax | grep slurm_monitor.py

Output:
>>> 10537 pts/1    S      0:00 python3 -u slurm_monitor.py
>>> 23426 pts/1    S+     0:00 grep --color=auto slurm_monitor.py

KILL PROCESS:
>>> kill 10537

OUTPUT STATUS:
>>> cd $DATA && find . -name status.txt -exec tail --lines=+2 {} \;

============================ ARGUMENTS ============================
'''

parser = argparse.ArgumentParser(description=description_text, formatter_class=argparse.RawTextHelpFormatter)
#parser._action_groups.pop()

required = parser.add_argument_group('required arguments')
additional = parser.add_argument_group('additional arguments')

required.add_argument('-n', dest='timesteps', nargs='?', type=int, required=True,
                    help='required timesteps                   (REQUIRED)')

additional.add_argument('--mail', dest='mail', nargs='?', type=str,
                    help='mail address (mail@server.de)        (default=None)')

additional.add_argument('--restart-mail', dest='bool', nargs='?', type=bool, default=True,
                    help='mail after every restart             (default=False)')

additional.add_argument('-b', '--backup', dest='backup', nargs='?', type=str,
                    help='backup location for files            (default=None)')

additional.add_argument('-m', '--monitorfile', dest='monitorFile', nargs='?', type=str, default='status.txt',
                    help='file with output from nohup command  (default=status.txt)')

additional.add_argument('-r', '--restartfile', dest='restartFile', nargs='?', type=str, default='FDS.dat',
                    help='restart file with data               (default=FDS.dat)')

additional.add_argument('-j', '--job', dest='jobscriptFile', nargs='?', type=str, default='jobscript',
                    help='jobscript to run SLURM job           (default=jobscript)')

additional.add_argument('--job-name', dest='jobname', nargs='?', type=str, default='Name',
                    help='job name not longer than 8 character (default=Name)')

additional.add_argument('--ntask-per-node', dest='tasks', nargs='?', type=str, default='32',
                    help='MPI task per node                    (default=32)')

additional.add_argument('--nodes', dest='nodes', nargs='?', type=str, default='3',
                    help='number of nodes                      (default=3)')

additional.add_argument('--time', dest='walltime', nargs='?', type=str, default='0-24:00:00',
                    help='walltime of server (d-hh:mm:ss)      (default=0-24:00:00)')







#parser.add_argument('text', action='store', type=str, help='The text to parse.')

args = parser.parse_args()
print(args.jobscriptFile)
print(args.jobname)
print(args.remail)