%% Script to convert from previous BG spike train input format to the new:
pathprefix = '\\139.91.162.90\cluster\stefanos\Documents\Glia\OLDBGST';
% Before is the original code up to 21/1/2016:
% Now create a single file/array for each cell to aid round-robin algorithm
% in cluster and memmory efficiency:
num_inha=2;
num_exc=2;
nruns = 100;
nPC = 700;
nPV = 233;
%disconnect from run struct to save memory
for ru = 63:4:nruns
    % Parse the old spike train data:
    fid = fopen( sprintf('\\\\139.91.162.90\\cluster\\stefanos\\Documents\\GitHub\\prefrontal-micro\\experiment\\network\\experiment_12\\importBackgroundStimParams_run_%04d.hoc',ru), 'r' );
    cac = textscan( fid, '%s', 'Delimiter', '\n' );
    % size(cac{1},1) must equals the # of rows in your data file.
    totalRows = size(cac{1},1);
    fclose( fid );

    % for total rows in data:
    BGb = cell(6,nPC);
    BGa = cell(6,nPC);
    BGap = cell(6,nPC);
    BGs = cell(2,nPV);
    tmpvectmat = [];
    insideVector = false;
    fprintf('Parsing file.\n');
    for k=1:totalRows
        fprintf('%f\n',k/totalRows);
    %     fprintf('Parsing data on row %d of %d...\n',k,totalRows);
        currentRow = cac{1}{k,1};
    %     fprintf('Row contains:\n%s\n',currentRow);
        %if line creates a new vector:
        if (~isempty(strfind(currentRow,'new Vector('))) && (insideVector)
            if ~isempty(strfind(eachRowElement{1},'BG_Stim_basal'))
                BGb{synid+1,cellid+1} = tmpvectmat;
            elseif ~isempty(strfind(eachRowElement{1},'BG_Stim_Apicpr'))
                BGap{synid+1,cellid+1} = tmpvectmat;
            elseif ~isempty(strfind(eachRowElement{1},'BG_Stim_Apic'))
                BGa{synid+1,cellid+1} = tmpvectmat;
            elseif ~isempty(strfind(eachRowElement{1},'BG_Stim_SomaPV'))
                BGs{synid+1,cellid+1+nPC} = tmpvectmat;
            end
            % parse location:
            tmp = strsplit(currentRow, ' = ');
            tmploc = strsplit(tmp{1}, '[');
            cellid = tmploc{3};
            cellid = str2double(cellid(1:end-1));
            synid = tmploc{4};
            synid = str2double(synid(1:end-1));

            tmpvectmat = [];
            insideVector = true;
        elseif ~isempty(strfind(currentRow,'new Vector('))
            insideVector = true;
            % parse location:
            tmp = strsplit(currentRow, ' = ');
            tmploc = strsplit(tmp{1}, '[');
            cellid = tmploc{3};
            cellid = str2double(cellid(1:end-1));
            synid = tmploc{4};
            synid = str2double(synid(1:end-1));
            continue;
        end
        if insideVector
            % pare thn ep omenh grammh:
            eachRowElement = strsplit(currentRow, ' = ');
            if length(eachRowElement)>1
                tmpvectmat = [tmpvectmat str2num(eachRowElement{2})];        
            end        
        end
    end
    fprintf('DONE parsing file.\n');
    for c=1:nPC+nPV
        % different number for pyramidals
        if c > nPC
            numberOfSyns = num_inha;
        else
            numberOfSyns = num_exc*3;
        end
        % Create, name and save a single file for that cell:
        fid = fopen(sprintf('%s/BG_cellgid%03d_run_%03d.hoc',pathprefix,c-0,ru),'W');
        fprintf(fid,'// Object decleration:\n');
        if c > nPC
            fprintf(fid,'objref BGs%03d[%d]\n',....
                c-1,num_inha);
            for k=1:numberOfSyns
                spikes = BGs{k,c};
                fprintf(fid,'BGs%03d[%d]=new Vector(%d)\n',c-1,k-1,length(spikes) );
                for j=1:length(spikes)
                    fprintf(fid,'BGs%03d[%d].x[%d]=%d\n',c-1,k-1,j-1, abs(spikes(j)));
                end
            end
        else
            fprintf(fid,'objref BGb%03d[%d],BGap%03d[%d],BGa%03d[%d]\n',....
                c-1,num_exc,...
                c-1,num_exc,...
                c-1,num_exc);
            for k=1:round(numberOfSyns/3)
                spikes = BGb{k,c};
                fprintf(fid,'BGb%03d[%d]=new Vector(%d)\n',c-1,k-1,length(spikes) );
                for j=1:length(spikes)
                    fprintf(fid,'BGb%03d[%d].x[%d]=%d\n',c-1,k-1,j-1, abs(spikes(j)));
                end
            end
            for k=1:round(numberOfSyns/3)
                spikes = BGap{k,c};
                fprintf(fid,'BGap%03d[%d]=new Vector(%d)\n',c-1,k-1,length(spikes) );
                for j=1:length(spikes)
                    fprintf(fid,'BGap%03d[%d].x[%d]=%d\n',c-1,k-1,j-1, abs(spikes(j)));
                end
            end
            for k=1:round(numberOfSyns/3)
                spikes = BGa{k,c};
                fprintf(fid,'BGa%03d[%d]=new Vector(%d)\n',c-1,k-1,length(spikes) );
                for j=1:length(spikes)
                    fprintf(fid,'BGa%03d[%d].x[%d]=%d\n',c-1,k-1,j-1, abs(spikes(j)));
                end
            end
        end
        fclose(fid);
    end
    save(sprintf('%s/RUN%04d.mat',pathprefix,ru),'BGb','BGa','BGap','BGs');
end
