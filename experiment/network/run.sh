#!/usr/bin/sh ##Make sure that before each run git repo is clean, so ##one can track each run to its source code: ## ANSI escape codes:
ERR='\033[0;31m'
WARN='\033[0;33m'
INFO='\033[0;32m'
NOC='\033[0m'

## I care about auto-commiting before each run. Since between runs I normaly modify some files, check
## ONLY for modified files and commit them. I should remember to manually commit if I make major changes.
##Documentation in: https://git-scm.com/docs/git-status
dirtygit="0"
##Files MODIFIED since index:
##Exclude run script(s) from check:
#echo "Excluding /run.sh"
##Update, include run.sh because why not.
dirtygit=$(( $dirtygit + $(git status --porcelain 2>/dev/null| grep "^ M" | wc -l) ))
#dirtygit=$(( $dirtygit + $(git status --porcelain 2>/dev/null| grep "^ M" | grep -v "/run.sh" | wc -l) ))
##Files DELETED since index:
#dirtygit=$(( $dirtygit + $(git status --porcelain 2>/dev/null| grep "^ D" | wc -l) ))
##Files added to the index, but uncommitted
#dirtygit=$(( $dirtygit + $(git status --porcelain 2>/dev/null| grep "^A" | wc -l) ))
##Check for new/untracked files (its better safe than sorry):
#dirtygit=$(( $dirtygit + $(git status --porcelain 2>/dev/null| grep "^\?" | wc -l) ))

git_branch=`git rev-parse --abbrev-ref HEAD`
#Commit only if in 'runs' branch (to avoid spamming) and if git is dirty:
if [[ ($dirtygit > 0 && $git_branch == "runs") ]]; then
	#echo -e "${WARN}The git repo is dirty. You have been warned...${NOC}"
	#echo -e "${INFO}Current HEAD is: ${gitsha1}${NOC}"
	echo -e "${INFO}Auto-committing dirty repo:${NOC}"
	# Auto commit only modified/deleted files:
	eval "git commit -am 'RUN AUTOCOMMIT'"
	#exit 1
	gitsha1=`git rev-parse HEAD`
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
nodes="78" ##52##jobname="STR_N100_S6_STC0" jobstdout=""
cluster="0"
# 0=Random, 1=Structured
exp="1"
## Serial number of network (RNG) in MATLAB:
sn="11"
## Ean einai clustered oi synapseeis stous dendrites:
clustbias="1"
## Excitation /inhibition bias (multiplier factor) gia PC2PC synapses
## for both NMDA AMPA
excitbias="1"
inhibias="1"
## ONly NMDA bias 
nmdabias="2.0"
ampabias="1.0"
## only GABAb
gababfactor="1"
## Gia to Background bias (alla to exw sbhsei)
BGe="10"
BGi="1"
## Posa stimulus synapses bazw
stimmagnitude="40"
## Stimulus frequency:
stimfreq="50"
## Default Elimination of Reciprocal Factor:
erf="0.0"
## Eliminate reciprocal Specific (defined in networkparameters.hoc)
ers="1"
## Default NMDA decay tau:
## afto grafei sto network.hoc to tau decay tou NMDA!!
nmdatau="90"
## NMDA beta is a factor that shifts the lognormal function to the left 
## if negative (minus sign is inside its mod file) so greater values
## enhance gNMDA.
nmdaflag="0"
## number of dendritic (basal) segments 
dendnseg="5"
## How many clusters I do identify 
Cl="7"
## How much (normalized) the weights are squashed into a narrow uniform distribution ## Zero means original lognormal weights dist; One means quantized 0.5 weights (connectivity of pairs
## 	stays the same.
Fs="1"
startRun="0"
endRun="0"
VARPID="0.5"
custom_jobs=(41 42)
erf_array=(0.0 0.1 0.3 0.5 0.7 0.9)
#naming convention in ten characters:

mechanisms="mechanism_simple"

if [ "$parallel" == "1" ]; then

jobstdout="$jobstdout\\\n========================================================================================"
jobstdout="$jobstdout\\\nNEURON MPIRUN starting at $(date)"
jobstdout="$jobstdout\\\n========================================================================================"
jobstdout="$jobstdout\\\nJOB NAME IS: $jobname"
jobstdout="$jobstdout\\\n========================================================================================"


##for run in $(seq $startRun $endRun); do
#for run in "${custom_jobs[@]}"
run="0"
#Move inhibitory synapses at different dendritic locations to check for more states:
ipid="0.05"
# cluster dendritic input to a single point on dendrite:
# locpid : 1=proximal skewed, 2=median normal, 3=distal skewed
locpid="7"
clpid="0.45"
stimfreq="60"
stimnoise="0.0"
inhibias="1"
excitbias="1"
gababfactor="1"
pv2pc="4"
pc2pc="100"
# Pass simulation stop externally in seconds:
tstop_sec="2.5"
#for pc2pc in $(seq 62 2 120); do
#for pv2pc in $(seq 52 52); do
for cluster in $(seq 0 0); do
##for gababfactor in $(seq 26 34); do
#for excitbias in $(seq 25 25); do
for inhibias  in $(seq 1 1 4); do
##for erf in "${erf_array[@]}"; do
#	cluster="${run}"
	if [ "$exp" == "1" ]; then
		exp_str="S"
		#jobname="distally_EB$(printf '%.3f' $excitbias)_IB$(printf '%.3f' $inhibias)_ST${stimmagnitude}_BGE${BGe}_BGI${BGi}_GBF$(printf '%.3f' $gababfactor)_Fs$(printf '%.3f' $Fs)_Cl${Cl}_NDS${dendnseg}_NMDAFLAG$(printf '%.3f' $nmdaflag)_CLB$(printf '%.3f' $clustbias)_Ss4c${cluster}_SN${sn}_r"
		# ran after army
		#jobname="ctrI50_ERF$(printf '%.1f' $erf)_EB$(printf '%.3f' $excitbias)_IB$(printf '%.3f' $inhibias)_ST${stimmagnitude}_GBF$(printf '%.3f' $gababfactor)_NMDAb$(printf '%.3f' $nmdabias)_${exp_str}s7c${cluster}_SN${sn}_r"
		#jobname="NFAi_ctrI50_EB$(printf '%.3f' $excitbias)_IB$(printf '%.3f' $inhibias)_ST${stimmagnitude}_GBF$(printf '%.3f' $gababfactor)_NMDAb$(printf '%.3f' $nmdabias)_Ab$(printf '%.3f' $ampabias)_${exp_str}s7c${cluster}_SN${sn}_r"
		#jobname="ERS${ers}_FiSF${stimfreq}_ctrI50_EB$(printf '%.3f' $excitbias)_IB$(printf '%.3f' $inhibias)_ST${stimmagnitude}_GBF$(printf '%.3f' $gababfactor)_NMDAb$(printf '%.3f' $nmdabias)_Ab$(printf '%.3f' $ampabias)_${exp_str}s7c${cluster}_SN${sn}_r"
		#jobname="test_SF${stimfreq}_IPID${ipid}ctrI50_EB$(printf '%.3f' $excitbias)_IB$(printf '%.3f' $inhibias)_ST${stimmagnitude}_GBF$(printf '%.3f' $gababfactor)_NMDAb$(printf '%.3f' $nmdabias)_Ab$(printf '%.3f' $ampabias)_${exp_str}s7c${cluster}_SN${sn}_r"
		jobname="noMg_heavy_SF${stimfreq}SN${stimnoise}_pc2pc${pc2pc}ctrI50_EB$(printf '%.3f' $excitbias)_IB$(printf '%.3f' $inhibias)_ST${stimmagnitude}_GBF$(printf '%.3f' $gababfactor)_NMDAb$(printf '%.3f' $nmdabias)_Ab$(printf '%.3f' $ampabias)_${exp_str}s${tstop_sec}c${cluster}_SN${sn}_r"
	else
		exp_str="R"
		jobname="NFiSF${stimfreq}_ctrI50_EB$(printf '%.3f' $excitbias)_IB$(printf '%.3f' $inhibias)_ST${stimmagnitude}_GBF$(printf '%.3f' $gababfactor)_NMDAb$(printf '%.3f' $nmdabias)_Ab$(printf '%.3f' $ampabias)_${exp_str}s7c${cluster}_SN${sn}_r"
		#jobname="ctrI50_ERF$(printf '%.1f' $erf)_EB$(printf '%.3f' $excitbias)_IB$(printf '%.3f' $inhibias)_ST${stimmagnitude}_GBF$(printf '%.3f' $gababfactor)_NMDAb$(printf '%.3f' $nmdabias)_${exp_str}s7c${cluster}_SN${sn}_r"
		#jobname="F_ctrI50_EB$(printf '%.3f' $excitbias)_IB$(printf '%.3f' $inhibias)_ST${stimmagnitude}_GBF$(printf '%.3f' $gababfactor)_NMDAb$(printf '%.3f' $nmdabias)_Rs7c${cluster}_SN${sn}_r"
	fi
	uniquejobname="${jobname}${run}"
	outputFile=$uniquejobname.out
	outputDir="${simglia}${uniquejobname}"
	echo "Output Job directory is:"
	echo $outputDir
	if [ -d $outputDir ]; then
		echo "${WARN}Job directory already exists. Stopping before overriding data.${NOC}"
		#exit 1
	else
		mkdir -p $outputDir;
	fi
	## Submit as Job in Sun Grid Engine:
	if [ "$schedule" == "1" ]; then
	echo -e "${INFO}SCHEDULER VERSION IS COMMENCING ${NOC}"
	qsub -b y -S /bin/bash -V -N $uniquejobname -o "${outputDir}/${outputFile}" -j y -pe orte 78-$nodes -p 0 -R y /opt/openmpi/bin/mpirun /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/myspecial ${nrn_repository} -nobanner -mpi \
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
	-c "TSTOP_SEC=$tstop_sec" \
	-c "CLUSTER_ID=$cluster" \
	-c "EXPERIMENT=$exp" \
	-c "CLUSTBIAS=$clustbias" \
	-c "EXCITBIAS=$excitbias" \
	-c "INHIBIAS=$inhibias" \
	-c "NMDABIAS=$nmdabias" \
	-c "AMPABIAS=$ampabias" \
	-c "ST=$stimmagnitude" \
	-c "SF=$stimfreq" \
	-c "STIMNOISE=$stimnoise" \
	-c "FS=$Fs" \
	-c "CL=$Cl" \
	-c "BGE=$BGe" \
	-c "BGI=$BGi" \
	-c "ERF=$erf" \
	-c "ERS=$ers" \
	-c "IPID=$ipid" \
	-c "CLPID=$clpid" \
	-c "LOCPID=$locpid" \
	-c "NMDATAU=$nmdatau" \
	-c "NMDA_FLAG=$nmdaflag" \
	-c "DEND_NSEG=$dendnseg" \
	-c "PV2PCsyns=$pv2pc" \
	-c "PC2PCsyns=$pc2pc" \
	-c "GABABFACTOR=$gababfactor" \
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
	/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/myspecial ${nrn_repository} -nobanner -mpi \
	-c "RUN=$run" \
	-c 'execute1("strdef JOBNAME, JOBDIR, GITSHA1, SN, SIMHOME, SIMGLIA")' \
	-c 'execute1("SN = \"'$sn'\"")' \
	-c 'execute1("SNd = '$sn'")' \
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
	-c "NMDABIAS=$nmdabias" \
	-c "ST=$stimmagnitude" \
	-c "FS=$Fs" \
	-c "CL=$Cl" \
	-c "BGE=$BGe" \
	-c "BGI=$BGi" \
	-c "NMDA_FLAG=$nmdaflag" \
	-c "DEND_NSEG=$dendnseg" \
	-c "GABABFACTOR=$gababfactor" \
	-c "VARPID=$VARPID" \
	/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/final.hoc
	fi
	
done
done
##done
##done
 
else

	echo -e "${INFO}NO MPI VERSION IS COMMENCING${NOC}"
nohup ../../$mechanisms/myspecial -nobanner -c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" final.hoc | tee nohup.out &

fi
echo Run script reached end.
