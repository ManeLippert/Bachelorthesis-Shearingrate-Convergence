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
#ssh btrzx1-1.rz.uni-bayreuth.de

# ssh in specfic folder
ssh -t btrzx1-1.rz.uni-bayreuth.de "cd /scratch/bt712347/ && bash --login"

# Disconnect VPN
#if [[ "$CONNECTION" = "" ]]; then
#    {
#    nmcli con down Universität\ Bayreuth
#    } &> /dev/null
#fi
