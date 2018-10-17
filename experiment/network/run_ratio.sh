#!/usr/bin/sh 
simhome=`pwd`
simhome="${simhome}/"
simglia="/home/cluster/stefanos/Documents/Glia"
simglia="${simglia}/"

nmda_fast="1"
schedule="0"
for synno in $(seq 5 5); do
	for expampa in $(seq 0 1); do
		if [ "$schedule" == "1" ]; then
			jobname="NMDA_AMPA_ratio_nsyn${synno}_expampa${expampa}"
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
			-c "SYN_NO=$synno" \
			-c "EXPAMPA=$expampa" \
			-c "NMDA_FAST=$nmda_fast" \
			-c "CALIB_SCRIPT=1" \
			/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/NMDA_AMPA_ratio.hoc
		else
			../../mechanism_simple/x86_64/special -nobanner \
			-c "SYN_NO=$synno" \
			-c "EXPAMPA=$expampa" \
			-c "NMDA_FAST=$nmda_fast" \
			-c "CALIB_SCRIPT=1" \
			NMDA_AMPA_ratio.hoc
		fi	
	done
done


