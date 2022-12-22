#!/bin/bash

# VARIABLES ======================================================================================

# Linux
#CONNECTION="$(nmcli con show --active | grep -i eduroam)"

# MacOS
CONNECTION="$(/System/Library/PrivateFrameworks/Apple80211.framework/Resources/airport -I | awk -F: '/ SSID/{print$2}')"
VPN="vpn-server.uni-bayreuth.de"
DIR="/scratch/bt712347/"

# CONNECT TO SSH =================================================================================

# Connect VPN
if [[ "$CONNECTION" != "eduroam" ]]; then
    {
        # Linux (Adjust VPN name if necessary)
        # nmcli con up $VPN

        # MacOS (anyconnect client (installed) and credentials (in home folder) must be defined)
        cat ~/.anyconnect_credentials | /opt/cisco/anyconnect/bin/vpn -s connect $VPN
    } &> /dev/null
fi

#ssh btrzx2-1.rz.uni-bayreuth.de

# ssh in specfic folder
ssh -t btrzx2-1.rz.uni-bayreuth.de "cd $DIR && bash --login"

# Disconnect VPN
if [[ "$CONNECTION" != "eduroam" ]]; then
    {
        # Linux (Adjust VPN name if necessary)
        #nmcli con down $VPN

        # MacOS
        /opt/cisco/anyconnect/bin/vpn disconnect
    } &> /dev/null
fi
