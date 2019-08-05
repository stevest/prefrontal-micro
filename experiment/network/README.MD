To run the model:
1. Open the run.sh
2. Edit your specific configuration details. Instructions exist in each section
3. Run the simulation with:$ sh run.sh
4. Collect the results inside the simulation output folder ($jobname variable)
5. Convert the simulation output to HDF5, utilizing the load_raw_batch.py and
   passing as the only argument the full simulation output path.
6. For results analysis please refer to the github repo:
   https://github.com/stamatiad/prefrontal_analysis.git

Regarding any questions please email me at stamatiad.st@gmail.com
