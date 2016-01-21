#!/bin/bash
##Make sure that before each run git repo is clean, so
##one can track each run to its source code:
##Documentation in: https://git-scm.com/docs/git-status
dirtygit="0"
##Files MODIFIED since index:
##Exclude run script(s) from check:
echo "Excluding /run.sh"
dirtygit=$(( $dirtygit + $(git status --porcelain 2>/dev/null| grep "^ M" | grep -v "/run.sh" | wc -l) ))
##Files DELETED since index:
dirtygit=$(( $dirtygit + $(git status --porcelain 2>/dev/null| grep "^ D" | wc -l) ))
##Files added to the index, but uncommitted
dirtygit=$(( $dirtygit + $(git status --porcelain 2>/dev/null| grep "^A" | wc -l) ))
##Check for new/untracked files (its better safe than sorry):
dirtygit=$(( $dirtygit + $(git status --porcelain 2>/dev/null| grep "^\?" | wc -l) ))


if [[ $dirtygit > 0 ]]; then
	echo "Can not continue with run when git repo is dirty. Exiting..."
	exit 1
else
	echo "Git repo is clean. Continue run with SHA1: "
	echo `git rev-parse HEAD`
fi

source /opt/gridengine/default/common/settings.csh
##$ -S /bin/sh
##$ -V
###$ -cwd ##Does NOT HAVE the /home/ prefix; causing error
##$ -N JobStefanos
##$ -o /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/flow.out -j y
##$ -l h=!compute-0-1&!compute-0-3&!compute-0-4&!compute-0-6&!compute-0-7&!compute-0-12&!compute-0-13&!compute-0-14&!compute-0-15
##$ -pe orte 90
##$ -R y
##sge_o_workdir="/home/cluster/stefanos/Documents/Github/prefrontal-micro/experiment/network/"

cd /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network
echo `pwd`



parallel="1"
nodes="24"
##jobname="STR_N100_S6_STC0"
jobstdout=""
simplified="1"
cluster="1"
# 0=Random, 1=Structured
exp="0"
state="7"
id="12"
sn="16"
vclamp="0.0"
binary="1"
clustbias="0.0"
startRun="0"
endRun="99"
if [ "$exp" == "1" ]; then
	jobname="STR_N100_S20_STC${cluster}"
else
	jobname="RND_N100_S20_STC${cluster}"
fi

if [ "$simplified" == "1" ]; then
	mechanisms="mechanism_simple"
else
	mechanisms="mechanism_complex"
fi

echo mechanisms folder is $mechanisms

if [ "$parallel" == "1" ]; then

#nohup mpirun -v -n 8 ../../$mechanisms/x86_64/special -mpi -nobanner -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" final.hoc | tee nohup.out & 

jobstdout="$jobstdout\\\n========================================================================================"
jobstdout="$jobstdout\\\nNEURON MPIRUN starting at $(date)"
jobstdout="$jobstdout\\\n========================================================================================"
jobstdout="$jobstdout\\\nJOB NAME IS: $jobname"
jobstdout="$jobstdout\\\n========================================================================================"

##mpirun -v -n 2 --host compute-0-1,compute-0-3 hostname; source ~/.bash_profile;echo $LD_LIBRARY_PATH; echo $PATH; nrniv testmpi.hoc 
#previous working run#mpirun -v -n 92 --host compute-0-18,compute-0-19,compute-0-20,compute-0-22 /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/x86_64/special -nobanner -mpi -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" final.hoc 


#function sigusr1handler()
#{
#        echo "SIGUSR1 caught by shell script" 1>&2
#}
 
#function sigusr2handler()
#{
#        echo "SIGUSR2 caught by shell script" 1>&2
#}
 
#trap sigusr1handler SIGUSR1
#trap sigusr2handler SIGUSR2

##parentheses1='execute1\(\"strdef STDOUT\"\)'
##qrsh -pe orte $nodes mpirun -v -np $nodes --mca plm_base_verbose 1 hostname; echo $HOSTNAME; /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/x86_64/special -nobanner -mpi -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" final.hoc 

#Working with STDOUT:#qsub -b y -S /bin/bash -V -N $jobname -o /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/$outputFile.out -j y -pe orte $nodes -R y /opt/openmpi/bin/mpirun /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/x86_64/special -nobanner -mpi -c "RUN=$run" -c "execute1\(\\\"'strdef STDOUT, JOBNAME'\\\"\)" -c "execute1\(\\\"'STDOUT = \\\"$jobstdout\\\"'\\\"\)" -c "execute1\(\\\"'JOBNAME = \\\"$jobname\\\"'\\\"\)" -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/final.hoc 
echo $jobstdout
#POSIXLY_CORRECT=0

for run in $(seq $startRun $endRun);
do
	echo $run;
	uniquejobname="${jobname}_${run}"
	outputFile=$uniquejobname.out
	## Submit as Job in Sun Grid Engine:
	##qsub -b y -S /bin/bash -V -N $uniquejobname -o /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/SGEoutput/$outputFile -j y -pe orte $nodes -R y /opt/openmpi/bin/mpirun /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/x86_64/special -nobanner -mpi -c "RUN=$run" -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/final.hoc 
	## Run locally:
	##/opt/openmpi/bin/mpirun -v -n 6 /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/x86_64/special -nobanner -mpi -c "RUN=$run" -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/final.hoc 
	#qsub -b y -S /bin/bash -V -N postjob -pe orte 1 -hold_jid $uniquejobname postjob.sh $outputFile $uniquejobname
done
#function run_neuron()
#{
#	/opt/openmpi/bin/mpirun /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/x86_64/special -nobanner -mpi -c "execute1\(\\\"'strdef STDOUT'\\\"\)" -c "execute1\(\\\"'STDOUT = \\\"$jobstdout\\\"'\\\"\)" -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/final.hoc 
#}
#qsub -b y -S /bin/bash -V -N $jobname -o /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/$jobname.out -j y -pe orte $nodes -R y run_neuron
 
#ssh compute-0-1 'bash -s' < kill_neurons.sh 
#ssh compute-0-3 'bash -s' < kill_neurons.sh 
#ssh compute-0-4 'bash -s' < kill_neurons.sh 
#ssh compute-0-6 'bash -s' < kill_neurons.sh 
#ssh compute-0-7 'bash -s' < kill_neurons.sh 
#ssh compute-0-12 'bash -s' < kill_neurons.sh 
#ssh compute-0-13 'bash -s' < kill_neurons.sh 
#ssh compute-0-14 'bash -s' < kill_neurons.sh 
#ssh compute-0-15 'bash -s' < kill_neurons.sh 
#ssh compute-0-18 'bash -s' < kill_neurons.sh 
#ssh compute-0-19 'bash -s' < kill_neurons.sh 
#ssh compute-0-20 'bash -s' < kill_neurons.sh 
#ssh compute-0-22 'bash -s' < kill_neurons.sh 

else

nohup ../../$mechanisms/x86_64/special -nobanner -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" final.hoc | tee nohup.out &

fi
echo Run script reached end.
