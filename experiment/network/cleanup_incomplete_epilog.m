function cleanup_incomplete_epilog()
% cd /home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/
% matlab -nojvm -nodisplay -nosplash -r "cleanup_incomplete_epilog; exit();";

% pathbase = '\\139.91.162.90\cluster\stefanos\Documents\Glia\';
pathbase = '/home/cluster/stefanos/Documents/Glia/';
diary('blah1.txt')
diary on;

d = dir(fullfile(pathbase));
rfn = {d([d.isdir]).name};
fn_idx = false(size(rfn));
for ib=2:4
    for myerf = [0,0.1:0.2:0.9]
fn_idx = fn_idx | cellfun(@(x) ~isempty(strfind(x,sprintf('ctrI50_ERF%0.1f_EB25.000_IB%0.3f_ST40_GBF15.000_NMDAb2.000_',myerf,ib))),rfn);
    end
end

check_flds = sort(rfn(fn_idx));

% folderto = sprintf('distallyRNDcl_EB25.000_IB2.000_ST40_GBF15.000_NMDAb2.000_Rs10c17_SN11_r1');

for foldr = 1:length(check_flds)
    folderto = check_flds{foldr};
    outfile =  sprintf('%s.out',folderto);
    
    % Parse output file to make sure that the epilog matlab script was
    % terminated improperly:
    fprintf('---------------------------------------------------------\n');
    fprintf('Output file:%s\n',outfile);
    fprintf('Loading output file...');
    fid = fopen(fullfile(pathbase,folderto,outfile),'r');
    cac = textscan( fid, '%s', 'Delimiter', '\n' );
    totalRows = size(cac{1},1);
    fclose( fid );
    fprintf('DONE\n');
    
    epilog_finished_properly = false;
    for k=totalRows:-1:1
        % search only the part of the (Huge!) output file where the epilog.m
        % ought to have run.
        currentRow = cac{1}{k,1};
        
        % Look for epilog terminating string:
        if strfind(currentRow,'Cell batch loaded')
            epilog_finished_properly = true;
            fprintf('epilog_finished_properly = true\n');
            break;
        end
        
        if strfind(currentRow,'Running epilog script...')
            fprintf('Reached end of search.\n');
            fprintf('epilog_finished_properly = false\n');
            break;
        end
        
    end
    
    
    
    
    if ~epilog_finished_properly
        fprintf('Reloading raw files if possible.\n');
        N=250;
        d = dir(fullfile(pathbase,folderto));
        rf = {d(~[d.isdir]).name};
        
        %     reload if possible
%         if vsoma.mat does not exist:
        if ~exist(fullfile(pathbase,folderto,'vsoma.mat'),'file')
            % Handle pyramidal somatic voltage:
            files_v_soma = sort(rf(cellfun(@(x) ~isempty(strfind(x,'v_soma')),rf)));
            vsoma = cell(length(files_v_soma),1);
            tic;
            for fn = 1:length(files_v_soma)
                filename = fullfile(pathbase,folderto,files_v_soma{fn});
                vtrace = load(filename);
                % dt equals 0.1:
            %     vtrace = vtrace(1:10:end);
                vsoma{fn} = vtrace;
            end
            fprintf('Time: %f secs\n',toc);
            outputfile = fullfile(pathbase,folderto,'vsoma.mat');
            save(outputfile, 'vsoma');
        end

        if ~exist(fullfile(pathbase,folderto,'nmdaA_pair.mat'),'file')
            % Handle pyramidal to pyramidal NMDA_A for each pair:
            files_NMDA_A_pair = sort(rf(cellfun(@(x) ~isempty(strfind(x,'nmdaA_src_')),rf)));
            if ~isempty(files_NMDA_A_pair)
                nmdaA_pair = cell(N,N);
                tic;
                for fn = 1:length(files_NMDA_A_pair)
                    %get files per node id:
                    filename = fullfile(pathbase,folderto,files_NMDA_A_pair{fn});
                    NMDA_A_pair = load(filename);
                    tmp = strsplit(files_NMDA_A_pair{fn}, {'nmdaA_src_','_trg_','.txt'} );
                    src = str2double(tmp{2})+1;
                    trg = str2double(tmp{3})+1;
                    NMDA_A_pair = NMDA_A_pair(1:10:end);
                    nmdaA_pair{src,trg} = NMDA_A_pair;
                end
                fprintf('Time: %f secs\n',toc);
                outputfile = fullfile(pathbase,folderto,'nmdaA_pair.mat');
                save(outputfile, 'nmdaA_pair');
            else
                fprintf('No NMDA A raw files detected.\n');
            end
        end
        
        if ~exist(fullfile(pathbase,folderto,'inmda.mat'),'file')
            % Handle pyramidal dendritic iNMDA:
            files_iNMDA_dend = sort(rf(cellfun(@(x) ~isempty(strfind(x,'inmda_dend')),rf)));
            if ~isempty(files_iNMDA_dend)
                inmda = cell(N,1);
                tic;
                for fn = 1:length(files_iNMDA_dend)
                    filename = fullfile(pathbase,folderto,files_iNMDA_dend{fn});
                    iNMDA = load(filename);
                    tmp = strsplit(files_iNMDA_dend{fn}, {'inmda_dend_','.'} );
                    trg = str2double(tmp{2})+1;
                    % dt equals 0.1:
                    iNMDA = iNMDA(1:10:end);
                    inmda{trg,1} = iNMDA;
                end
                fprintf('Time: %f secs\n',toc);
                outputfile = fullfile(pathbase,folderto,'inmda.mat');
                save(outputfile, 'inmda');
            end
        end
        
        if ~exist(fullfile(pathbase,folderto,'vdend.mat'),'file')
            % Handle pyramidal dendritic voltage:
            files_v_dend = sort(rf(cellfun(@(x) ~isempty(strfind(x,'v_dend')),rf)));
            vdend = cell(N,5);
            tic;
            for fn = 1:length(files_v_dend)
                %filter the dend segment!
                tmp = strsplit(files_v_dend{fn}, {'_','.'} );
                seg = str2double(tmp{3})+1;
                dend = str2double(tmp{4})+1;
                filename = fullfile(pathbase,folderto,files_v_dend{fn});
                vtrace = load(filename);
                % dt equals 0.1:
                %         vtrace = vtrace(1:10:end);
                vdend{dend,seg} = vtrace;
            end
            fprintf('Time: %f secs\n',toc);
            outputfile = fullfile(pathbase,folderto,'vdend.mat');
            save(outputfile, 'vdend');
        end
    
        
    end
    
    %     clear up the text files and you are OK to go.
    
    % txtfn = sort(rf(cellfun(@(x) ~isempty(strfind(x,'.txt')),rf)));
    %
    % for fn = 1:length(txtfn)
    %     filename = fullfile(pathbase,folderto,txtfn{fn});
    %     delete(filename);
    % end
    % fprintf('Files TXT deleted (Time: %f secs)!\n',toc);
    
    lscmd = ['find ',fullfile(pathbase,folderto),'/ -type f -name ''*.txt'' '];
    delcmd = ['find ',fullfile(pathbase,folderto),'/ -type f -name ''*.txt'' -delete'];
    
    tic;
    [status,cmdout] = system(lscmd);
    fprintf('Status returned: %g\n',status);
    fprintf('Deleting files:\n%s\n',cmdout);
    [status,cmdout] = system(delcmd);
    % cmdout
    fprintf('Files TXT deleted (Time: %f secs)!\n',toc);
    fprintf('\n');
    
end
diary off;
