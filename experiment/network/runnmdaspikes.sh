#!/usr/bin/sh 
for nmdatau in $(seq 10 10); do
for isi in $(seq 1 1); do
for number in $(seq 1 1); do
for synapses in $(seq 1 1); do
../../mechanism_simple/x86_64/special -nobanner -c "SYNAPSES=$synapses" -c "ISI=$isi" -c "NUMBER=$number" -c "NMDATAU=$nmdatau" nmda_spikes.hoc 
done
done
done
done

