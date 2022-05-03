#!/bin/bash

# This file collects useful functions which are used by several of the
# other scripts shipped with GKW.

# Use it by "sourcing" it:
# Either write
#   . $GKW_HOME/scripts/gkwlib.sh
# or
#   source $GKW_HOME/scripts/gkwlib.sh
# at the beginning of your script.

find_latest_executable(){
    #after calling this function, the variable $latest_executable is set.

    # need gkw home
    if [ -z "${GKW_HOME}" ]; then
	echo "You need to set GKW_HOME!"
	exit 1
    fi

    
    # get the network node hostname and strip everything behind the first dot
    HN=`uname -n | sed -e 's/\..*$//'`
    # get the first three letters
    HN3=$(echo ${HN} | awk '{ string=substr($1, 1, 3); print string; }')
    
    latest_executable=$(ls -t ${GKW_HOME}/run/gkw*${HN3}*.x | head -1)

    if [ -z "${latest_executable}" ]; then
	echo "No suitable executables found for ${HN3}"
	exit 1
    fi
}

