#!/bin/bash
##Make sure that before each run git repo is clean, so 
##one can track each run to its source code:
## ANSI escape codes:
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
echo "Currently at directory:"
echo $simhome
echo `pwd`

parallel="1"
## Use scheduler or directly run with mpi:
schedule="1"
#All nodes are:312 
nodes="12" ##52##jobname="STR_N100_S6_STC0" 
cluster="0"
# 0=Random, 1=Structured
exp="1"
## Serial number of network (RNG) in MATLAB:
sn="1"
## Ean einai clustered oi synapseeis stous dendrites:
clustbias="1"
## Excitation /inhibition bias (multiplier factor) gia PC2PC synapses
## for both NMDA AMPA
excitbias="1"
inhibias="1"
## ONly NMDA bias (default is 10)
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


#==============================================================================
#==============================================================================
# RUN PARAMETERS:

run="0"
#Move inhibitory synapses at different dendritic locations to check for more states:
ipid="0.05"
# cluster dendritic input to a single point on dendrite:
# locpid : 1=proximal skewed, 2=median normal, 3=distal skewed
locpid="7"
clpid="0.45"
stimfreq="60"
stimnoise="0.5"
inhibias="3"
excitbias="1"
gababfactor="15"
pv2pc="4"
pc2pc="10"
no_mg="0"
# Pass simulation stop externally in seconds:
tstop_sec="0.01"
loccl="1"
cluster="0"
if [ "$exp" == "1" ]; then
	exp_str="S"
else
	exp_str="R"
fi

jobname="NEW_LC${loccl}_SN${stimnoise}_pc2pc${pc2pc}_pv2pc${pv2pc}_EB$(printf '%.3f' $excitbias)_IB$(printf '%.3f' $inhibias)_GBF$(printf '%.3f' $gababfactor)_NMDAb$(printf '%.3f' $nmdabias)_AMPAb$(printf '%.3f' $ampabias)_${exp_str}s${tstop_sec}c${cluster}_SN${sn}_r"

uniquejobname="${jobname}${run}"
outputFile=$uniquejobname.out
outputDir="${simglia}${uniquejobname}"
echo "Output Job directory is:"
echo $outputDir
if [ -d $outputDir ]; then
	echo "${WARN}Job directory already exists. Possible data overwrite.${NOC}"
else
	mkdir -p $outputDir;
fi



jobstdout=""
jobstdout="${jobstdout}\\\n========================================================================================"
jobstdout="${jobstdout}\\\nNEURON SIMULATION starting at $(date)"
jobstdout="${jobstdout}\\\n========================================================================================"
jobstdout="${jobstdout}\\\nJOB NAME IS: ${jobname}"
jobstdout="${jobstdout}\\\n========================================================================================"

# Group all run parameters in one variable to use across running scenarios:
run_variables=( 
"execute1(\"strdef JOBNAME, JOBDIR, GITSHA1, SN, SIMHOME, SIMGLIA\")" 
"execute1(\"SN = \\\"$sn\\\"\")" 
"execute1(\"GITSHA1 = \\\"$gitsha1\\\"\")" 
"execute1(\"JOBNAME = \\\"$uniquejobname\\\"\")" 
"execute1(\"JOBDIR = \\\"$outputDir\\\"\")" 
"execute1(\"SIMHOME = \\\"$simhome\\\"\")" 
"execute1(\"SIMGLIA = \\\"$simglia\\\"\")" 
"RUN=$run" 
"PARALLEL=$parallel" 
"TSTOP_SEC=$tstop_sec" 
"CLUSTER_ID=$cluster" 
"EXPERIMENT=$exp" 
"CLUSTBIAS=$clustbias" 
"EXCITBIAS=$excitbias" 
"INHIBIAS=$inhibias" 
"NMDABIAS=$nmdabias" 
"AMPABIAS=$ampabias" 
"ST=$stimmagnitude" 
"SF=$stimfreq" 
"STIMNOISE=$stimnoise" 
"FS=$Fs" 
"CL=$Cl" 
"LOCCL=$loccl" 
"BGE=$BGe" 
"BGI=$BGi" 
"ERF=$erf" 
"ERS=$ers" 
"IPID=$ipid" 
"CLPID=$clpid" 
"LOCPID=$locpid" 
"NMDATAU=$nmdatau" 
"NMDA_FLAG=$nmdaflag" 
"DEND_NSEG=$dendnseg" 
"PV2PCsyns=$pv2pc" 
"PC2PCsyns=$pc2pc"
"GABABFACTOR=$gababfactor" 
"VARPID=$VARPID" 
"NO_MG=$no_mg" 
)

# Interleave variables with NEURON -c flag:
wrapped_vars=()
for E in "${run_variables[@]}"; do
    wrapped_vars+=("-c")
    wrapped_vars+=("${E}")
done

# Parenthesis need escaping in SGE version:
wrapped_vars_sge=()
for E in "${run_variables[@]}"; do
    wrapped_vars_sge+=("-c")
    wrapped_vars_sge+=("'${E}'")
done


##for gababfactor in $(seq 26 34); do
for loccl in $(seq 1 1); do
#for stimnoise in $(seq 0.1 0.1 0.5); do
for cluster in $(seq 0 0); do

#==============================================================================
#==============================================================================
## EXECUTE SIMULATION:
if [ "$parallel" == "1" ]; then
		## Submit as Job in Sun Grid Engine:
		if [ "$schedule" == "1" ]; then
			echo -e "${INFO}SCHEDULED MPI VERSION IS COMMENCING ${NOC}"
			#qsub -b y -S /bin/bash -V -N $uniquejobname -o "${outputDir}/${outputFile}" -j y -pe orte 24-$nodes -p 0 -R y /opt/openmpi/bin/mpirun /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/$mechanisms/myspecial ${nrn_repository} -nobanner -mpi \
			#$run_variables \
			#/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/final.hoc 
			qsub -b y -S /bin/bash -V -N ${uniquejobname} -o ${outputDir}/${outputFile} -j y -pe orte 24-${nodes} -p 0 -R y /opt/openmpi/bin/mpirun /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/mechanism_simple/myspecial -nobanner -mpi "${wrapped_vars_sge[@]}" /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/final.hoc

		else
			## Run local mpi without scheduler:
			echo -e "${INFO}NO SCHEDULED MPI VERSION IS COMMENCING${NOC}"
			# Are those necessary?
			#PATH=/home/stefanos/Libraries/${nrn_repository}/x86_64/bin:$PATH
			#LD_LIBRARY_PATH=/home/stefanos/Libraries/${nrn_repository}/x86_64/lib:$LD_LIBRARY_PATH
			#export PATH
			#export LD_LIBRARY_PATH
			#echo "NEURON executable in PATH before mpirun is: "
			#echo `which nrniv`

			/opt/openmpi/bin/mpirun -np ${nodes} -x PATH -x LD_LIBRARY_PATH /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/mechanism_simple/myspecial -nobanner -mpi "${wrapped_vars[@]}" /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/final.hoc
		fi
		
else
	# If not $parallel
	echo -e "${INFO}NO MPI VERSION IS COMMENCING${NOC}"
	#-c "PARALLEL=$parallel" -c "SIMPLIFIED=$simplified" -c "CLUSTER_ID=$cluster" -c "EXPERIMENT=$exp" -c "ST=$state" -c "ID=$id" -c "SN=$sn" -c "VCLAMP=$vclamp" -c "ISBINARY=$binary" -c "CLUSTBIAS=$clustbias" 
	#cmd="nohup ../../${mechanisms}/myspecial -nobanner ${run_variables} final.hoc | tee nohup.out &"
	#echo "${run_variables[@]/#/-c }" -

	#/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/mechanism_simple/myspecial ${nrn_repository} -nobanner -mpi ${wrapped_vars[@]} /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/final.hoc
	/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/mechanism_simple/myspecial "${wrapped_vars[@]}" test.hoc -
fi


done
done

#==============================================================================
#==============================================================================
# SCRIPT END
echo Run script reached end.
