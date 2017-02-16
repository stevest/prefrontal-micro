function [vsoma]=load_raw_batch_agnostic(pathto,varargin)
% Load experiment batch..
% ...GIVEN that one run is inside the pathto folder.
% One run is equivalent to one MPI batch job (in its own folder!).

% show all warning messages:
warning('on','all');

r = dir(pathto);
if isempty(r)
    return;
end
rf = {r(~[r.isdir]).name};

% Handle pyramidal somatic voltage:
files_v_soma = sort(rf(cellfun(@(x) ~isempty(strfind(x,'v_soma')),rf)));
vsoma = cell(length(files_v_soma),1);
tic;
for fn = 1:length(files_v_soma)
    filename = fullfile(pathto,files_v_soma{fn});
    vtrace = load(filename);
    % dt equals 0.1:
    vtrace = vtrace(1:10:end);
    vsoma{fn} = vtrace;
end
fprintf('Time: %f secs\n',toc);
outputfile = fullfile(pathto,'vsoma.mat');
save(outputfile, 'vsoma');
% remove files to save space and sanity:
tic;
for fn = 1:length(files_v_soma)
    filename = fullfile(pathto,files_v_soma{fn});
    delete(filename);	
end
fprintf('Files v_soma deleted (Time: %f secs)!\n',toc);

% Handle pyramidal dendritic voltage:
files_v_dend = sort(rf(cellfun(@(x) ~isempty(strfind(x,'v_dend')),rf)));
if ~isempty(files_v_dend)
    vdend = cell(length(files_v_soma),1);
    tic;
    for fn = 1:length(files_v_dend)
        %filter the dend segment!
        tmp = strsplit(files_v_dend{fn}, {'_','.'} );
        seg = str2double(tmp{3})+1;
        dend = str2double(tmp{4})+1;
        filename = fullfile(pathto,files_v_dend{fn});
        vtrace = load(filename);
        % dt equals 0.1:
        vtrace = vtrace(1:10:end);
        vdend{dend,seg} = vtrace;
    end
    fprintf('Time: %f secs\n',toc);
    outputfile = fullfile(pathto,'vdend.mat');
    save(outputfile, 'vdend');
    tic;
    % remove files to save space and sanity:
    for fn = 1:length(files_v_dend)
        filename = fullfile(pathto,files_v_dend{fn});
        delete(filename);
    end
    fprintf('Files v_dend deleted (Time: %f secs)!\n',toc);
end

% Handle pyramidal to pyramidal synaptic locations:
files_pid_pyramidal = sort(rf(cellfun(@(x) ~isempty(strfind(x,'pid_pyramidal')),rf)));
if ~isempty(files_pid_pyramidal)
    pidpyramidal = cell(length(files_v_soma),length(files_v_soma));
    tic;
    for fn = 1:length(files_pid_pyramidal)
        %get files per node id:
        filename = fullfile(pathto,files_pid_pyramidal{fn});
        % Parse the old spike train data:
        fid = fopen( filename, 'r' );
        fcont = textscan( fid, '%s', 'Delimiter', '\n' );
        totalRows = size(fcont{1},1);
        fclose( fid );
        for crn = 1:totalRows
            cr = fcont{1}{crn,1};
            if strfind(cr,'src=')
                crdump = strsplit(cr, {'src=','trg=',' '} );
                src = str2double(crdump{2})+1;
                trg = str2double(crdump{3})+1;
                continue;
            end
            % each normal line MUST be just a float!
            if( any(isstrprop(cr, 'digit')) )
                pid = str2double(cr);
                % check for empty lines:
                %         if pid >= 0 && ~isnan(pid)
                if ~exist('src','var') || ~exist('trg','var')
                    warning('src/trg nonexisting in line %d, file %s !',crn,filename);
                else
                    pidpyramidal{src,trg} = [pidpyramidal{src,trg} pid];
                end
            end
        end
    end
    fprintf('Time: %f secs\n',toc);
    outputfile = fullfile(pathto,'pidpyramidal.mat');
    save(outputfile, 'pidpyramidal');
    % remove files to save space and sanity:
    tic;
    for fn = 1:length(files_pid_pyramidal)
        filename = fullfile(pathto,files_pid_pyramidal{fn});
        delete(filename);
    end
    fprintf('Files files_pid_pyramidal deleted (Time: %f secs)!\n',toc);
end

% Handle pyramidal to opyramidal delays:
files_delay_pyramidal = sort(rf(cellfun(@(x) ~isempty(strfind(x,'delay_pyramidal')),rf)));
if ~isempty(files_delay_pyramidal)
    delaypyramidal = cell(length(files_v_soma),length(files_v_soma));
    tic;
    for fn = 1:length(files_delay_pyramidal)
        %get files per node id:
        filename = fullfile(pathto,files_delay_pyramidal{fn});
        % Parse the old spike train data:
        fid = fopen( filename, 'r' );
        fcont = textscan( fid, '%s', 'Delimiter', '\n' );
        totalRows = size(fcont{1},1);
        fclose( fid );
        for crn = 1:totalRows
            cr = fcont{1}{crn,1};
            if strfind(cr,'src=')
                crdump = strsplit(cr, {'src=','trg=',' '} );
                src = str2double(crdump{2})+1;
                trg = str2double(crdump{3})+1;
                continue;
            end
            % each normal line MUST be just a float!
            if( any(isstrprop(cr, 'digit')) )
                delay = str2double(cr);
                % check for empty lines:
                %         if delay >= 0 && ~isnan(delay)
                delaypyramidal{src,trg} = [delaypyramidal{src,trg} delay];
            else
                %             if isfloat(delay) && ~isnan(delay)
                %                 warning('Negative delay read from file for cellid %d!!',trg);
                %             end
            end
        end
    end
    fprintf('Time: %f secs\n',toc);
    outputfile = fullfile(pathto,'delaypyramidal.mat');
    save(outputfile, 'delaypyramidal');
    % remove files to save space and sanity:
    tic;
    for fn = 1:length(files_delay_pyramidal)
        filename = fullfile(pathto,files_delay_pyramidal{fn});
        delete(filename);
    end
    fprintf('Files files_delay_pyramidal deleted (Time: %f secs)!\n',toc);
end

% Handle pyramidal dendritic iNMDA:
files_iNMDA_dend = sort(rf(cellfun(@(x) ~isempty(strfind(x,'inmda_dend')),rf)));
if ~isempty(files_iNMDA_dend)
    inmda = cell(length(files_iNMDA_dend),1);
    tic;
    for fn = 1:length(files_iNMDA_dend)
        filename = fullfile(pathto,files_iNMDA_dend{fn});
        iNMDA = load(filename);
        tmp = strsplit(files_iNMDA_dend{fn}, {'inmda_dend_','.'} );
        trg = str2double(tmp{2})+1;
        % dt equals 0.1:
        iNMDA = iNMDA(1:10:end);
        inmda{trg,1} = iNMDA;
    end
    fprintf('Time: %f secs\n',toc);
    outputfile = fullfile(pathto,'inmda.mat');
    save(outputfile, 'inmda');
    % remove files to save space and sanity:
    tic;
    for fn = 1:length(files_iNMDA_dend)
        filename = fullfile(pathto,files_iNMDA_dend{fn});
        delete(filename);
    end
    fprintf('Files iNMDA_dend deleted (Time: %f secs)!\n',toc);
end

end
