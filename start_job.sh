#!/bin/bash

JOBFINISHED=false

while [[ $JOBFINISHED = false ]]; do
    # Check if restart file is written
    if [[ -f "FDS.dat" ]];then 
        DATA=$(awk '!/^#/&&NF{printf tolower($0" ")}' FDS.dat)

        NTREMAIN="$(cut -d',' -f4 <<<"$DATA")"

        NTREMAINVALUE="$(cut -d'=' -f2 <<<"$NTREMAIN")"
        NTREMAINVALUE="$(echo -e "${NTREMAINVALUE}" | tr -d '[:space:]')"

        if [[ "$NTREMAINVALUE" = "0" ]];then
            JOBFINISHED=true
        else
            sbatch jobscript*
            # Wait 24h (walltime) and 5min (additional) to check again
            sleep 86700 
        fi
    else
        sbatch jobscript*
        # Wait 24h (walltime) and 5min (additional) to check again
        sleep 86700
    fi
done