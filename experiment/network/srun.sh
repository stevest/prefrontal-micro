#!/bin/bash

simhome="/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network"
storagehome="/home/cluster/stefanos/Documents/Glia"
cd $simhome

jobstdout=""

mechanisms="mechanism_simple"

## run in single core:
/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/x86_64/special /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/singlecell.hoc -
