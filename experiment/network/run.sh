#!/usr/bin/sh ##Make sure that before each run git repo is clean, so ##one can track each run to its source code:
## ANSI escape codes:
ERR='\033[0;31m'
WARN='\033[0;33m'
INFO='\033[0;32m'
NOC='\033[0m'

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


	gitsha1=`git rev-parse HEAD`
if [[ $dirtygit > 0 ]]; then
	echo -e "${ERR}Can not continue with run when git repo is dirty. Exiting...${NOC}"
	echo -e "${INFO}Current HEAD is: ${gitsha1}${NOC}"
	#exit 1
else
	gitsha1=`git rev-parse HEAD`
	echo -e ${INFO}"Git repo is clean. Continue run with SHA1: ${gitsha1}${NOC}"
fi

##source /opt/gridengine/default/common/settings.csh
simhome=`pwd`
simhome="${simhome}/"
simglia="/home/cluster/stefanos/Documents/Glia"
simglia="${simglia}/"
## Define neuron repo. This is the git repo inside ~/Libraries.
nrn_repository="nrn"
echo "Currently at directory:"
echo $simhome
echo `pwd`
echo "Using NEURON from repository: ${nrn_repository}."

parallel="1"
## Use scheduler or directly run with mpi:
schedule="1"
#All nodes are:312 
nodes="312" ##jobname="STR_N100_S6_STC0" jobstdout=""
cluster="6"
# 0=Random, 1=Structured
#!!! MAJOR CAUTION CHANGE W.*0.7 values tmp for random experiment!!!
exp="1"
#!!! MAJOR CAUTION CHANGE W.*0.7 values tmp for random experiment!!!
sn="10"
clustbias="1"
excitbias="1"
inhibias="1"
Cl="7"
Fs="1"
startRun="0"
endRun="0"
VARPID="0.5"
custom_jobs=(41 42)
#naming convention in ten characters:

mechanisms="mechanism_simple"


if [ "$parallel" == "1" ]; then

jobstdout="$jobstdout\\\n========================================================================================"
jobstdout="$jobstdout\\\nNEURON MPIRUN starting at $(date)"
jobstdout="$jobstdout\\\n========================================================================================"
jobstdout="$jobstdout\\\nJOB NAME IS: $jobname"
jobstdout="$jobstdout\\\n========================================================================================"


for run in $(seq $startRun $endRun); do
#for run in "${custom_jobs[@]}"
run="0"
for BGe in $(seq 1 1); do
	for BGi in $(seq 1 1); do
#	cluster="${run}"
	if [ "$exp" == "1" ]; then
		jobname="control_excx${excitbias}_inhx${inhibias}_BGE${BGe}_BGI${BGi}_Ss4c${cluster}_SN${sn}_r"
	else
		jobname="control_excx${excitbias}_inhx${inhibias}_BGE${BGe}_BGI${BGi}_Rs4c${cluster}_SN${sn}_r"
	fi
	uniquejobname="${jobname}${run}"
	outputFile=$uniquejobname.out
	outputDir="${simglia}${uniquejobname}"
	echo "Output Job directory is:"
	echo $outputDir
	if [ -d $outputDir ]; then
		echo "Job directory already exists. Stopping before overriding data."
		#exit 1
	else
		mkdir -p $outputDir;
	fi
	## Submit as Job in Sun Grid Engine:
	if [ "$schedule" == "1" ]; then
	echo -e "${INFO}SCHEDULER VERSION IS COMMENCING ${NOC}"
	qsub -b y -S /bin/bash -V -N $uniquejobname -o "${outputDir}/${outputFile}" -j y -pe orte 1-$nodes -p -3 -R y /opt/openmpi/bin/mpirun /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/myspecial ${nrn_repository} -nobanner -mpi \
	-c "RUN=$run" \
	-c "execute1\(\\\"'strdef JOBNAME, JOBDIR, GITSHA1, SN, SIMHOME, SIMGLIA'\\\"\)" \
	-c "execute1\(\\\"'SN = \\\"$sn\\\"'\\\"\)" \
	-c "execute1\(\\\"'SNd = $sn'\\\"\)" \
	-c "execute1\(\\\"'GITSHA1 = \\\"$gitsha1\\\"'\\\"\)" \
	-c "execute1\(\\\"'JOBNAME = \\\"$uniquejobname\\\"'\\\"\)" \
	-c "execute1\(\\\"'JOBDIR = \\\"$outputDir\\\"'\\\"\)" \
	-c "execute1\(\\\"'SIMHOME = \\\"$simhome\\\"'\\\"\)" \
	-c "execute1\(\\\"'SIMGLIA = \\\"$simglia\\\"'\\\"\)" \
	-c "PARALLEL=$parallel" \
	-c "CLUSTER_ID=$cluster" \
	-c "EXPERIMENT=$exp" \
	-c "CLUSTBIAS=$clustbias" \
	-c "EXCITBIAS=$excitbias" \
	-c "INHIBIAS=$inhibias" \
	-c "FS=$Fs" \
	-c "CL=$Cl" \
	-c "BGE=$BGe" \
	-c "BGI=$BGi" \
	-c "VARPID=$VARPID" \
	/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/final.hoc 

	else
## Run local mpi without scheduler:
	echo -e "${INFO}NO SCHEDULED MPI VERSION IS COMMENCING${NOC}"
PATH=/home/stefanos/Libraries/${nrn_repository}/x86_64/bin:$PATH
LD_LIBRARY_PATH=/home/stefanos/Libraries/${nrn_repository}/x86_64/lib:$LD_LIBRARY_PATH
export PATH
export LD_LIBRARY_PATH
echo "NEURON executable in PATH before mpirun is: "
echo `which nrniv`

/opt/openmpi/bin/mpirun -np 1 \
	-x PATH \
	-x LD_LIBRARY_PATH \
	/home/stefanos/Documents/SourceTree/prefrontal-micro/$mechanisms/myspecial ${nrn_repository} -nobanner -mpi \
	-c "RUN=$run" \
	-c 'execute1("strdef JOBNAME, JOBDIR, GITSHA1, SN, SIMHOME, SIMGLIA")' \
	-c 'execute1("SN = \"'$sn'\"")' \
	-c 'execute1("GITSHA1 = \"'$gitsha1'\"")' \
	-c 'execute1("JOBNAME = \"'$uniquejobname'\"")' \
	-c 'execute1("JOBDIR = \"'$outputDir'\"")' \
	-c 'execute1("SIMHOME = \"'$simhome'\"")' \
	-c 'execute1("SIMGLIA = \"'$simglia'\"")' \
	-c "PARALLEL=$parallel" \
	-c "CLUSTER_ID=$cluster" \
	-c "EXPERIMENT=$exp" \
	-c "CLUSTBIAS=$clustbias" \
	-c "EXCITBIAS=$excitbias" \
	-c "INHIBIAS=$inhibias" \
	-c "FS=$Fs" \
	-c "CL=$Cl" \
	-c "BGE=$BGe" \
	-c "BGI=$BGi" \
	-c "VARPID=$VARPID" \
	/home/stefanos/Documents/SourceTree/prefrontal-micro/experiment/network/final.hoc
	fi
	
done
done
done
 
else

	echo -e "${INFO}NO MPI VERSION IS COMMENCING${NOC}"
nohup ../../$mechanisms/myspecial -nobanner -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" final.hoc | tee nohup.out &

fi
echo Run script reached end.
