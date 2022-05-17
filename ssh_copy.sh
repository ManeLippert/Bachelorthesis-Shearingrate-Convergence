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
    read -p "Remote dir: " REMOTEDIR MAKEDIR
    if ssh btrzx1-1.rz.uni-bayreuth.de "[ -d $REMOTEDIR ]"; then
        REMOTEINPUT=true
    # Make folder if necessary
    elif [[ "$MAKEDIR" = "m" ]]; then
        if [[ "$DEST" = "remote" ]] || [[ "$DEST" = "r" ]]; then
            ssh btrzx1-1.rz.uni-bayreuth.de "mkdir $REMOTEDIR"
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
    read -p "Local  dir: " LOCALDIR MAKEDIR   
    if [ -d "$LOCALDIR" ]; then
        LOCALINPUT=true
    # Make folder if necessary
    elif [[ "$MAKEDIR" = "m" ]]; then
        if [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
            mkdir $LOCALDIR
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
    ssh btrzx1-1.rz.uni-bayreuth.de "mkdir temp"
    scp -r $LOCALDIR bt712347@btrzx1-1.rz.uni-bayreuth.de:temp
    ssh btrzx1-1.rz.uni-bayreuth.de "mv temp/$LOCALDIR/* $REMOTEDIR/"
    ssh btrzx1-1.rz.uni-bayreuth.de "rm -rf temp"
elif [[ "$DEST" = "local" ]] || [[ "$DEST" = "l" ]]; then
    mkdir temp
    scp -r bt712347@btrzx1-1.rz.uni-bayreuth.de:$REMOTEDIR temp
    mv temp/$REMOTEDIR/* $LOCALDIR/
    rm -rf temp
fi

# Disconnect VPN
if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con down Universität\ Bayreuth
    } &> /dev/null
fi