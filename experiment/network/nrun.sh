#!/bin/bash
mpirun -n 24 -mca btl ^openib ~/Desktop/prefrontal-micro/mechanism_simple/x86_64/special -mpi -nobanner -c "PARALLEL=$1" -c "SIMPLIFIED=$2" -c "CLUSTER_ID=$3" -c "EXPERIMENT=$4" -c "EXPID=$5" -c "AA=$6" -c "VCLAMP=$7" final.hoc
