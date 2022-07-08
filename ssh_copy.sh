#!/bin/bash

# Variables
CONNECTION="$(nmcli con show --active | grep -i eduroam)"
DESTINPUT=false
REMOTEINPUT=false
LOCALINPUT=false
DIRINPUT=false

# Connect VPN
if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con up Universität\ Bayreuth
    } &> /dev/null
fi

# Check input destination
while [[ $DESTINPUT = false ]]; do
    read -p "Destination: " DEST DESTADD
    if [[ "$DEST" = "remote" ]] ||  [[ "$DEST" = "r" ]]; then
        DESTINPUT=true
    elif [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
        DESTINPUT=true
    # Context
    else
        echo "Two options available: remote or local"
    fi
done

# When directories are different
if [[ "$DESTADD" = "d" ]];then
    read -p "Remote dir: " REMOTEDIR REMOTEADD
    read -p "Local  dir: " -e LOCALDIR LOCALADD

    # Copy files
    if [[ "$DEST" = "remote" ]] ||  [[ "$DEST" = "r" ]]; then
        rsync -a -P --exclude='data.h5' $LOCALDIR/* bt712347@btrzx1-1.rz.uni-bayreuth.de:$REMOTEDIR
        if [[ "$LOCALADD" = "j" ]] || [[ "$REMOTEADD" = "mj" ]] || [[ "$REMOTEADD" = "jm" ]]; then
            ssh btrzx1-1.rz.uni-bayreuth.de "cp -r gkw/run/* $REMOTEDIR/"
            ssh btrzx1-1.rz.uni-bayreuth.de "cd $REMOTEDIR && sbatch jobscript*"
        fi
    elif [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
        rsync -a -P --exclude={'gkw.*','DM*','slurm*','monitor_job.py'} bt712347@btrzx1-1.rz.uni-bayreuth.de:$REMOTEDIR/* $LOCALDIR
    fi
#When directories are equal
else
    read -p "Dir: " -e DIR DIRADD
    # Copy files
    if [[ "$DEST" = "remote" ]] ||  [[ "$DEST" = "r" ]]; then
        rsync -a -P --exclude='data.h5' $DIR/* bt712347@btrzx1-1.rz.uni-bayreuth.de:$DIR
        if [[ "$DIRADD" = "j" ]]; then
            ssh btrzx1-1.rz.uni-bayreuth.de "cp -r gkw/run/* $DIR/"
            ssh btrzx1-1.rz.uni-bayreuth.de "cd $DIR && nohup python3 -u monitor_job.py &"
        fi
    elif [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
        rsync -a -P --exclude={'gkw.*','DM*','slurm*','monitor_job.py'} bt712347@btrzx1-1.rz.uni-bayreuth.de:$DIR/* $DIR
    fi
fi

# Disconnect VPN
if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con down Universität\ Bayreuth
    } &> /dev/null
fi