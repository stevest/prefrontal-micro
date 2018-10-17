#!/usr/bin/sh 
for excitb in $(seq 0.3 0.2 0.9); do
for nsyn in $(seq 5 5 200); do
	../../mechanism_simple/x86_64/special -nobanner -c "NSYN=$nsyn" -c "EXCITB=$excitb" -c "NOJUMP=1" -c "HAVEAMPA=1" NMDA_no_Mg.hoc
done
done


