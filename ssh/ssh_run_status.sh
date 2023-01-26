#!/bin/bash

# Variables
#CONNECTION="$(nmcli con show --active | grep -i eduroam)"

# Connect VPN
#if [[ "$CONNECTION" = "" ]]; then
#    {
#    nmcli con up Universität\ Bayreuth
#    } &> /dev/null
#fi

#sleep 2
ssh btrzx1-1.rz.uni-bayreuth.de "squeue -u bt712347"
#ssh btrzx1-1.rz.uni-bayreuth.de "cd /scratch/bt712347/data && find . -name status.txt -exec cat {} \;"
#ssh btrzx1-1.rz.uni-bayreuth.de "cd /scratch/bt712347/data && find . -name status.txt -exec tail -8 {} \;"

# Disconnect VPN
#if [[ "$CONNECTION" = "" ]]; then
#    {
#    nmcli con down Universität\ Bayreuth
#    } &> /dev/null
#fi
