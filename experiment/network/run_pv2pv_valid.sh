#!/usr/bin/sh 
simhome=`pwd`
simhome="${simhome}/"
simglia="/home/cluster/stefanos/Documents/Glia"
simglia="${simglia}/"

../../mechanism_simple/x86_64/special -nobanner \
-c "CALIB_SCRIPT=1" \
pv2pv_validation.hoc 
