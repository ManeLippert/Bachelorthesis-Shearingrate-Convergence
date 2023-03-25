import os

#data = int(os.path.getmtime("gkwdata.h5"))
#restart_dat = int(os.path.getmtime("FDS.dat"))
#restart_file = int(os.path.getmtime("FDS"))
#
#print(data)
#print(restart_dat)
#print(restart_file)

#print(data - restart_file)

walltime = "0-24:00:00"

def get_time_in_seconds(time):
    
    # Format D-HH:MM:SS or HH:MM:SS or MM:SS
    
    time = time.replace("-", ":")
    time_split = time.split(":")

    seconds = [7*24*60*60,24*60*60, 60*60, 60, 1]
    seconds = seconds[-len(time_split):]

    time_sec = sum([a*b for a,b in zip(seconds, map(int,time_split))])

    return time_sec

print(get_time_in_seconds(walltime))


Test1 = False
Test2 = False

if Test1:
    print("Test1")
elif Test2: 
    print("TEst2")
else:
    print("nothing")