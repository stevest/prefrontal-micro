#!/bin/bash/

parallel="1"
simplified="0"
cluster="0"
exp="1"
state="1"
id="12"
sn="16"
vclamp="0.0"
binary="1"
clustbias="1.0"

if [ "$simplified" == "1"]; then
	mechanisms="mechanism_simple"
else
	mechanisms="mechanism_complex"
fi

if [ "$parallel" == "1" ]; then

mpirun -v -n 8 ../../$mechanisms/x86_64/special -mpi -nobanner -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" final.hoc 

else

../../$mechanisms/x86_64/special -nobanner -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" final.hoc

fi