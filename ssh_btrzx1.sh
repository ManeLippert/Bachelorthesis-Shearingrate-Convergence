#!/bin/bash

CONNECTION="$(nmcli con show --active | grep -i eduroam)"

if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con up Universität\ Bayreuth
    } &> /dev/null
fi

#sleep 1
ssh btrzx1-1.rz.uni-bayreuth.de

if [[ "$CONNECTION" = "" ]]; then
    {
    nmcli con down Universität\ Bayreuth
    } &> /dev/null
fi