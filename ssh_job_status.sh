#!/bin/bash

#!/bin/bash

# Variables
CONNECTION="$(nmcli con show --active | grep -i eduroam)"

# Connect VPN
if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con up Universität\ Bayreuth
    } &> /dev/null
fi

ssh btrzx1-1.rz.uni-bayreuth.de "squeue -u bt712347"

# Disconnect VPN
if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con down Universität\ Bayreuth
    } &> /dev/null
fi