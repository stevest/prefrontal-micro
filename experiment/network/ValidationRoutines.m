%% Multiple automated routines for validation


%% Validate NMDA/AMPA ratio and save it in .hoc file

% Paralelize the task in multiple jobs:
matlabpool open 12

% run jobs in background (paralellized by OS)
tic;
parfor cj = 0:55 % total morphologies
    %              validateRatio( MODE, PYID, RUN_NO, MINVAR, MAXVAR, MAXERR, TARGET, MAXITER)
    [weightNMDA(cj+1),scaleNMDA(cj+1)] = validateRatio( 1, cj, 1, 0, 5, 0.002, 0.060, 500) % run_no ALWAYS 1! (BUG)
    [weightAMPA(cj+1),scaleAMPA(cj+1)] = validateRatio( 2, cj, 1, 0, 5, 0.002, 0.045, 500)
end
sprintf('Time taken in seconds %f\n\n',toc)
save('ValidatedNMDA_AMPA_RATIO.mat', 'scaleNMDA', 'scaleAMPA', 'weightNMDA','weightAMPA');
matlabpool close

% Write NMDA results to a .hoc file:
fid = fopen('NMDA_Array.hoc','w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref nmdaArr\n');
fprintf(fid,'nmdaArr = new Vector(56)\n\n');
fprintf(fid,'// Validated Weights:\n');
for id=1:56
    fprintf(fid,'nmdaArr.x[%d] = NMDAfactor * %f\t//current is %f\n',id-1, weightNMDA(id), scaleNMDA(id));
end
fprintf(fid,'// END OF FILE\n');
fclose(fid);

% Write AMPA results to a .hoc file:
fid = fopen('AMPA_Array.hoc','w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref ampaArr\n');
fprintf(fid,'ampaArr = new Vector(56)\n\n');
fprintf(fid,'// Validated Weights:\n');
for id=1:56
    fprintf(fid,'ampaArr.x[%d] = AMPAfactor * %f\t//current is %f\n',id-1, weightAMPA(id), scaleAMPA(id));
end
fprintf(fid,'// END OF FILE\n');
fclose(fid);

