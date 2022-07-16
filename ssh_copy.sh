#!/bin/bash

# Variables
CONNECTION="$(nmcli con show --active | grep -i eduroam)"
DESTINPUT=false
REMOTEINPUT=false
LOCALINPUT=false
DIRINPUT=false

# Connect VPN
#if [[ "$CONNECTION" = "" ]]; then
#    {
#    nmcli con up Universität\ Bayreuth
#    } &> /dev/null
#fi

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
    # Check if remote directory exist
    while [[ $REMOTEINPUT = false ]]; do
        read -p "Remote dir: " REMOTEDIR REMOTEADD
        if ssh btrzx1-1.rz.uni-bayreuth.de "[ -d $REMOTEDIR ]"; then
            REMOTEINPUT=true
        # Make folder if necessary
        elif [[ "$REMOTEADD" = "m" ]] || [[ "$REMOTEADD" = "mj" ]] || [[ "$REMOTEADD" = "jm" ]]; then
            if [[ "$DEST" = "remote" ]] || [[ "$DEST" = "r" ]]; then
                ssh btrzx1-1.rz.uni-bayreuth.de "mkdir -p $REMOTEDIR"
                REMOTEINPUT=true
            else
                echo "Directory has to exist on remote machine"
            fi
        # Context
        else
            if [[ "$DEST" = "remote" ]] || [[ "$DEST" = "r" ]]; then
                echo -e "Directory does not exist on remote machine \n-> Create directory with 'dir m'"
            else
                echo "Directory does not exist on remote machine"
            fi
        fi
    done
    # Check if local directory exists
    while [[ $LOCALINPUT = false ]]; do
        read -p "Local  dir: " -e LOCALDIR LOCALADD
        if [ -d "$LOCALDIR" ]; then
            LOCALINPUT=true    
        # Make folder if necessary
        elif [[ "$LOCALADD" = "m" ]]; then
            if [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
                mkdir -p $LOCALDIR
                LOCALINPUT=true
            else
                echo "Directory has to exist on local machine"
            fi
        # Context
        else
            if [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
                echo -e "Directory does not exist on local machine \n-> Create directory with 'dir m'"
            else
                echo "Directory does not exist on local machine"
            fi
        fi
    done
    # Copy files
    if [[ "$DEST" = "remote" ]] ||  [[ "$DEST" = "r" ]]; then
        rsync -a -P --exclude='data.h5' $LOCALDIR/* bt712347@btrzx1-1.rz.uni-bayreuth.de:$REMOTEDIR
        if [[ "$LOCALADD" = "j" ]] || [[ "$REMOTEADD" = "mj" ]] || [[ "$REMOTEADD" = "jm" ]]; then
            ssh btrzx1-1.rz.uni-bayreuth.de "cp -r gkw/run/* $REMOTEDIR/"
            ssh btrzx1-1.rz.uni-bayreuth.de "cd $REMOTEDIR && nohup python3 -u monitor_job.py 1>/dev/null 2>/dev/null &"
        fi
    elif [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
        rsync -a -P --exclude={'gkw.*','DM*','FDS','slurm*','monitor_job.py'} bt712347@btrzx1-1.rz.uni-bayreuth.de:$REMOTEDIR/* $LOCALDIR
    fi
#When directories are equal
else
    # Check if directory exists
    while [[ $DIRINPUT = false ]]; do
        read -p "Dir: " -e DIR DIRADD
        
        if [ -d "$DIR" ] && ssh btrzx1-1.rz.uni-bayreuth.de "[ -d $DIR ]"; then
            DIRINPUT=true
        # Make folder if necessary on local machine
        elif [[ "$DIRADD" = "m" ]] || [[ "$DIRADD" = "mj" ]] || [[ "$DIRADD" = "jm" ]]; then
            # Local
            if [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
                mkdir -p $DIR
                DIRINPUT=true
            # Remote
            elif [[ "$DEST" = "remote" ]] || [[ "$DEST" = "r" ]]; then
                ssh btrzx1-1.rz.uni-bayreuth.de "mkdir -p $DIR"
                DIRINPUT=true
            else
                echo "Directory has to exist on machine"
            fi
        # Context
        else
            # Local
            if [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
                echo -e "Directory does not exist on local machine \n-> Create directory with 'dir m'"
            # Remote
            elif [[ "$DEST" = "remote" ]] || [[ "$DEST" = "r" ]]; then
                echo -e "Directory does not exist on remote machine \n-> Create directory with 'dir m'"
            else
                echo "Directory does not exist"
            fi
        fi
    done
    # Copy files
    if [[ "$DEST" = "remote" ]] ||  [[ "$DEST" = "r" ]]; then
        # Copy files in backup folder
        if [[ "$DIRADD" = "b" ]]; then
            rsync -a -P $DIR/* bt712347@btrzx1-1.rz.uni-bayreuth.de:/scratch/bt712347/backup/$DIR
        else
            rsync -a -P --exclude='data.h5' $DIR/* bt712347@btrzx1-1.rz.uni-bayreuth.de:$DIR
        fi
        # Run job
        if [[ "$DIRADD" = "j" ]] || [[ "$DIRADD" = "mj" ]] || [[ "$DIRADD" = "jm" ]]; then
            ssh btrzx1-1.rz.uni-bayreuth.de "cp -r gkw/run/* $DIR/"
            ssh btrzx1-1.rz.uni-bayreuth.de "cd $DIR && nohup python3 -u monitor_job.py &"
        fi
    elif [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
        # Copy files from backup folder
        if [[ "$DIRADD" = "b" ]]; then
            rsync -a -P --exclude={'gkw.*','DM*','FDS','slurm*','monitor_job.py'} bt712347@btrzx1-1.rz.uni-bayreuth.de:/scratch/bt712347/backup/$DIR/* $DIR
        else
            rsync -a -P --exclude={'gkw.*','DM*','FDS','slurm*','monitor_job.py'} bt712347@btrzx1-1.rz.uni-bayreuth.de:$DIR/* $DIR
        fi
    fi
fi

# Disconnect VPN
#if [[ "$CONNECTION" = "" ]]; then
#    {
#    nmcli con down Universität\ Bayreuth
#    } &> /dev/null
#fi