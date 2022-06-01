#!/bin/bash

# Variables
CONNECTION="$(nmcli con show --active | grep -i eduroam)"

# Connect VPN
if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con up Universität\ Bayreuth
    } &> /dev/null
fi

sleep 2
ssh btrzx1-1.rz.uni-bayreuth.de

# Disconnect VPN
if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con down Universität\ Bayreuth
    } &> /dev/null
fi