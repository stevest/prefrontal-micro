#!/bin/sh
##source /opt/gridengine/default/common/settings.csh
##$ -S /bin/sh
##$ -V
###$ -cwd ##Does NOT HAVE the /home/ prefix; causing error
##$ -N JobStefanos
##$ -o /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/flow.out -j y
##$ -l h=!compute-0-1&!compute-0-3&!compute-0-4&!compute-0-6&!compute-0-7&!compute-0-12&!compute-0-13&!compute-0-14&!compute-0-15
##$ -pe orte 90
##$ -R y
##sge_o_workdir="/home/cluster/stefanos/Documents/Github/prefrontal-micro/experiment/network/"

echo "========================================================================================"
echo "NEURON MPIRUN starting at $(date)"
echo "========================================================================================"

qsub -b y -S /bin/sh -V -N JobStefanos -o /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/flow.out -j y -pe orte 90 -R y /opt/openmpi/bin/mpirun /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/hellompi 
echo Run script reached end.
