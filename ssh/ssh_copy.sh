#!/bin/bash

# Messages ========================================================================================

help()
{
echo -e "\n
<--------------------- $(basename $0) --------------------->
\n
 -l, --local      : Set destination to local machine
 -r, --remote     : Set destination to remote machine
 -d, --dir        : Set directory 
                    (Same at local and remote machine)
 -rd, --remotedir : Set remote directory
 -ld, --localdir  : Set local directory
 -b               : Copy files from backup folder
 -m               : Create directory (also subfolders)
 -h, --help       : Outpu helper message      
\n"
}

# VARIABLES =======================================================================================

REMOTEBASE="/scratch/bt712347/"
#LOCALBASE="~/Bachelorthesis-Shearingrate-Wavelength"
LOCALBASE=""
VPN="vpn-server.uni-bayreuth.de"
SERVER="btrzx1-1.rz.uni-bayreuth.de"

MAKEDIR=false

# Linux
#CONNECTION="$(nmcli con show --active | grep -i eduroam)"

# MacOS
CONNECTION="$(/System/Library/PrivateFrameworks/Apple80211.framework/Resources/airport -I | awk -F: '/ SSID/{print$2}')"

# PARSER ==========================================================================================

PARAMS=""

if [[ "$#" -eq 0 ]]; then
    echo -e "\nERROR: Arguments needed."
    help
    exit 1;
fi

while (( "$#" )); do

    case "$1" in

        -m|--mkdir)
            MAKEDIR=true
            shift
            ;;

        -b|--backup)
            BACKUPDIR=true
            shift
            ;;

        -l|--local|-r|--remote)
            DEST="$1"
            shift
            ;;

        -d|--dir)
            if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
                LOCALDIR="$2"
                REMOTEDIR="$2"
                DIREQUAL=true
                shift 2
            else
                echo -e "\nERROR: Argument for $1 is missing." >&2
                exit 1
            fi
            ;;

        -ld|--localdir)
            if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
                LOCALDIR="$2"
                DIREQUAL=false
                shift 2
            else
                echo -e "\nERROR: Argument for $1 is missing." >&2
                exit 1
            fi
            ;;

        -rd|--remotedir)
            if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
                REMOTEDIR="$2"
                DIREQUAL=false
                shift 2
            else
                echo -e "\nERROR:   Argument for $1 is missing." >&2
                exit 1
            fi
            ;;

        -h|--help)
            help
            exit 0
            ;;

        -*|--*=)
            echo -e "\nERROR: Invalid command option."
            help
            exit 1
            ;;
    esac
done

eval set -- "$PARAMS"

# Check if local dir is defined
if [ -z ${LOCALDIR+x} ]; then
    echo -e "\nERROR: Directory for local machine has to be defined" 
    help
    exit 1
fi

# Check if remote dir is defined
if [ -z ${REMOTEDIR+x} ]; then 
    echo -e "\nERROR: Directory for remote machine has to be defined" 
    help
    exit 1
fi

# FOLDER ==========================================================================================

#cd $HOME$WORKDIR

# VPN =============================================================================================

# Connect VPN
#if [[ "$CONNECTION" != "eduroam" ]]; then
#    {
#        # Linux (Adjust VPN name if necessary)
#        # nmcli con up $VPN
#
#        # MacOS (anyconnect client (installed) and credentials (in home folder) must be defined)
#        cat ~/.anyconnect_credentials | /opt/cisco/anyconnect/bin/vpn -s connect $VPN
#    } &> /dev/null
#fi

# REOMOTE =========================================================================================

# Check folder exists
if ssh btrzx1-1.rz.uni-bayreuth.de "[ -d $REMOTEBASE$REMOTEDIR ]"; then
    :
# Make folder if necessary
elif [[ "$MAKEDIR" = true ]]; then
    if [[ "$DEST" = "--remote" ]] || [[ "$DEST" = "-r" ]]; then
        ssh btrzx1-1.rz.uni-bayreuth.de "mkdir -p $REMOTEBASE$REMOTEDIR"
    else
        echo -e "\nERROR: Directory has to exist on remote machine"
    fi
# Context
else
    if [[ "$DEST" = "--remote" ]] || [[ "$DEST" = "--r" ]]; then
        echo -e "\nERROR: Directory does not exist on remote machine"
        help
    else
        echo -e "\nERROR: Directory does not exist on remote machine"
    fi
fi

# LOCAL ===========================================================================================

# Check folder exists
if [ -d "$LOCALBASE$LOCALDIR" ]; then
    :   
# Make folder if necessary
elif [[ "$MAKEDIR" = true ]]; then
    if [[ "$DEST" = "--local" ]] || [[ "$DEST" = "-l" ]]; then
        mkdir -p $LOCALBASE$LOCALDIR
    else
        echo -e "\nERROR: Directory has to exist on local machine"
    fi
# Context
else
    if [[ "$DEST" = "--local" ]] || [[ "$DEST" = "-l" ]]; then
        echo -e "\nERROR: Directory does not exist on local machine"
        help
    else
        echo -e "\nERROR: Directory does not exist on local machine"
    fi
fi

# RSYNC ===========================================================================================

# Copy files
if [[ "$DEST" = "--remote" ]] ||  [[ "$DEST" = "-r" ]]; then
    # Copy files in backup folder
    if [[ "$BACKUPDIR" = true ]]; then
        rsync -a -P $LOCALBASE$LOCALDIR/* $SERVER:$REMOTEBASEbackup/$REMOTEDIR
    else
        rsync -a -P --exclude='data.h5' $LOCALBASE$LOCALDIR/* $SERVER:$REMOTEBASE$REMOTEDIR
    fi

elif [[ "$DEST" = "--local" ]] || [[ "$DEST" = "-l" ]]; then
    # Copy files from backup folder
    if [[ "$BACKUPDIR" = true ]]; then
        rsync -a -P --exclude={'gkw.*','DM*','FDS','slurm*','monitor_job.py'} $SERVER:$REMOTEBASEbackup/$REMOTEDIR/* $LOCALBASE$LOCALDIR
    else
        rsync -a -P --exclude={'gkw.*','DM*','FDS','slurm*','monitor_job.py'} $SERVER:$REMOTEBASE$REMOTEDIR/* $LOCALBASE$LOCALDIR
    fi
fi

# VPN =============================================================================================

# Disconnect VPN
#if [[ "$CONNECTION" != "eduroam" ]]; then
#    {
#        # Linux (Adjust VPN name if necessary)
#        #nmcli con down $VPN
#
#        # MacOS
#        /opt/cisco/anyconnect/bin/vpn disconnect
#    } &> /dev/null
#fi