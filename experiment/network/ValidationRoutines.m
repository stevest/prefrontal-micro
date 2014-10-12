%% Multiple automated routines for validation


%% Validate NMDA/AMPA ratio and save it in .hoc file

% Paralelize the task in multiple jobs:
matlabpool open 4

% run jobs in background (paralellized by OS)
% Ratios: 1,5 = 0.0675
% 1,33 = 0.06
% 1.11 = 0.05
tic;
parfor cj = 0:55 % total morphologies
    %              validateRatio( MODE, PYID, RUN_NO, MINVAR, MAXVAR, MAXERR, TARGET, MAXITER)
    [weightNMDA(cj+1),currentNMDA(cj+1)] = validateRatio( 1, cj, 50, 0, 5, 0.002, 0.0675, 20) % run_no ALWAYS 1! (BUG!)
    [weightAMPA(cj+1),currentAMPA(cj+1)] = validateRatio( 2, cj, 50, 0, 5, 0.002, 0.045, 20) % giati to ampa current to fit itan : 0.0164???
end
sprintf('Time taken in seconds %f\n\n',toc)
% save('ValidatedNMDA_AMPA_RATIO.mat', 'scaleNMDA', 'scaleAMPA', 'weightNMDA','weightAMPA');
matlabpool close

% Write NMDA results to a .hoc file:
fid = fopen('NMDA_Array_1_1.hoc','w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref nmdaArr\n');
fprintf(fid,'nmdaArr = new Vector(56)\n\n');
fprintf(fid,'// Validated Weights:\n');
for id=1:56
    fprintf(fid,'nmdaArr.x[%d] = NMDAfactor * %f\t//current is %f\n',id-1, weightNMDA(id), currentNMDA(id));
end
fprintf(fid,'// END OF FILE\n');
fclose(fid);

% % Write AMPA results to a .hoc file:
% fid = fopen('AMPA_Array_1_5.hoc','w');
% fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
% fprintf(fid,'// Object decleration:\n');
% fprintf(fid,'objref ampaArr\n');
% fprintf(fid,'ampaArr = new Vector(56)\n\n');
% fprintf(fid,'// Validated Weights:\n');
% for id=1:56
%     fprintf(fid,'ampaArr.x[%d] = AMPAfactor * %f\t//current is %f\n',id-1, weightAMPA(id), currentAMPA(id));
% end
% fprintf(fid,'// END OF FILE\n');
% fclose(fid);

