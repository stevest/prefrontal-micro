# !/bin/bash
ssh compute-0-1 'bash -s' < kill_neurons.sh
ssh compute-0-3 'bash -s' < kill_neurons.sh
ssh compute-0-4 'bash -s' < kill_neurons.sh
ssh compute-0-6 'bash -s' < kill_neurons.sh
ssh compute-0-7 'bash -s' < kill_neurons.sh
#ssh compute-0-12 'bash -s' < kill_neurons.sh
ssh compute-0-13 'bash -s' < kill_neurons.sh
ssh compute-0-14 'bash -s' < kill_neurons.sh
ssh compute-0-15 'bash -s' < kill_neurons.sh
ssh compute-0-18 'bash -s' < kill_neurons.sh
ssh compute-0-19 'bash -s' < kill_neurons.sh
ssh compute-0-20 'bash -s' < kill_neurons.sh
ssh compute-0-22 'bash -s' < kill_neurons.sh

