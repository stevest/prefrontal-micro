#jobstdout="blah \\\nde blah $(date)"
#echo $jobstdout
#t/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/mechanism_simple/x86_64/special -nobanner -c "execute1(\"strdef STDOUT\")" -c "execute1(\"STDOUT = \\\"$jobstdout\\\"\")" /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/testrun.hoc

#get job number inside job:
jobname="testJob"
nodes="1"
qsub -b y -S /bin/bash -V -N $jobname -o /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/${jobname}_${JOB_ID}.out -j y -pe orte $nodes -R y /opt/openmpi/bin/mpirun /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/mechanism_simple/x86_64/special -nobanner -mpi /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/testrun.hoc
