%% Parse stimulation matrix! WTF!


% Parse the old spike train data:
fid = fopen( sprintf('stimdata.txt',ru), 'r' );
cac = textscan( fid, '%s', 'Delimiter', '\n' );
% size(cac{1},1) must equals the # of rows in your data file.
totalRows = size(cac{1},1);
fclose( fid );

stimdata_str = zeros(100,7);
stimdata_rnd = zeros(100,7);


% for total rows in data:
fprintf('Parsing file.\n');
for k=1:totalRows
    fprintf('%f\n',k/totalRows);
    currentRow = cac{1}{k,1};
    
    if strfind(currentRow,'PcellStimListSTR.x[')
        tmp = strsplit(currentRow, ' = ');
        tmploc = strsplit(tmp{1}, '[');
        cellid = str2double(tmploc{2}(1:end-1))+1;
        clusterid = str2double(tmploc{3}(1:end-1))+1;
        stimdata_str(cellid,clusterid) = str2double(tmp{2})+1;
    elseif strfind(currentRow,'PcellStimListRND.x[')
        tmp = strsplit(currentRow, ' = ');
        tmploc = strsplit(tmp{1}, '[');
        cellid = str2double(tmploc{2}(1:end-1))+1;
        clusterid = str2double(tmploc{3}(1:end-1))+1;
        stimdata_rnd(cellid,clusterid) = str2double(tmp{2})+1;
    end

end
fprintf('DONE parsing file.\n');
