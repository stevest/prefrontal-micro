#!/usr/bin/sh 
simhome=`pwd`
simhome="${simhome}/"
simglia="/home/cluster/stefanos/Documents/Glia"
simglia="${simglia}/"

schedule="0"

for excitb in $(seq 1 1); do
for nsyn in $(seq 5 5 70); do
	
	if [ "$schedule" == "1" ]; then
	jobname="NMDA_fast_nsyn${nsyn}_excitb$(printf '%.3f' $excitb)_"
	uniquejobname="${jobname}${run}"
	outputFile=$uniquejobname.out
	outputDir="${simglia}${uniquejobname}"
	echo "Output Job directory is:"
	echo $outputDir
	if [ -d $outputDir ]; then
		echo "Job directory already exists!"
	else
		mkdir -p $outputDir;
	fi
	qsub -b y -S /bin/bash -V -N $uniquejobname -o "${outputDir}/${outputFile}" -j y -pe orte 1 -p 0 -R y /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/mechanism_simple/x86_64/special nrn -nobanner \
	-c "NSYN=$nsyn" \
	-c "EXCITB=$excitb" \
	-c "FAST=1" \
	-c "HAVEAMPA=1"  \
	/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/NMDA_fast.hoc
	else
	../../mechanism_simple/x86_64/special -nobanner \
	-c "NSYN=$nsyn"  \
	-c "EXCITB=$excitb"  \
	-c "FAST=0" \
	-c "HAVEAMPA=1" \
	NMDA_fast.hoc
	fi	
done
done


