#!/bin/bash

# Variables
CONNECTION="$(nmcli con show --active | grep -i eduroam)"
DESTINPUT=false
REMOTEINPUT=false
LOCALINPUT=false

# Connect VPN
if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con up Universität\ Bayreuth
    } &> /dev/null
fi

# Check input destination
while [[ $DESTINPUT = false ]]; do
    read -p "Destination: " DEST
    if [[ "$DEST" = "remote" ]] ||  [[ "$DEST" = "r" ]]; then
        DESTINPUT=true
    elif [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
        DESTINPUT=true
    # Context
    else
        echo "Two options available: remote or local"
    fi
done

# Check if remote directory exist
while [[ $REMOTEINPUT = false ]]; do
    read -p "Remote dir: " REMOTEDIR REMOTEADD
    if ssh btrzx1-1.rz.uni-bayreuth.de "[ -d $REMOTEDIR ]"; then
        REMOTEINPUT=true
    # Make folder if necessary
    elif [[ "$REMOTEADD" = "m" ]]; then
        if [[ "$COUNTREMOTE" = 1 ]]; then
            echo "No need to make a directory"
        else
            if [[ "$DEST" = "remote" ]] || [[ "$DEST" = "r" ]]; then
                ssh btrzx1-1.rz.uni-bayreuth.de "mkdir -p $REMOTEDIR"
                REMOTEINPUT=true
            else
                echo "Directory has to exist on remote machine"
            fi
        fi
    elif [[ "$REMOTEDIR" = "s" ]]; then
        REMOTEINPUT=true
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
    read -p "Local  dir: " LOCALDIR LOCALADD
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
    elif [[ "$LOCALDIR" = "s" ]]; then
        LOCALINPUT=true
    # Context
    else
        if [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
            echo -e "Directory does not exist on local machine \n-> Create directory with 'dir m'"
        else
            echo "Directory does not exist on local machine"
        fi
    fi
done

# If local directory should the same as remote directory
if [[ "$LOCALDIR" = "s" ]]; then
    LOCALDIR=$REMOTEDIR
elif [[ "$REMOTEDIR" = "s" ]]; then
    REMOTEDIR=$LOCALDIR
fi

# Copy files
if [[ "$DEST" = "remote" ]] ||  [[ "$DEST" = "r" ]]; then
    scp -r $LOCALDIR/* bt712347@btrzx1-1.rz.uni-bayreuth.de:$REMOTEDIR
    if [[ "$LOCALADD" = "j" ]]; then
        ssh btrzx1-1.rz.uni-bayreuth.de "cp -r gkw/run/* $REMOTEDIR/"
        ssh btrzx1-1.rz.uni-bayreuth.de "cd $REMOTEDIR && sbatch jobscript*"
    fi
elif [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
    scp -r bt712347@btrzx1-1.rz.uni-bayreuth.de:$REMOTEDIR/* $LOCALDIR
    rm -rf $LOCALDIR/gkw.*
fi

# Disconnect VPN
if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con down Universität\ Bayreuth
    } &> /dev/null
fi