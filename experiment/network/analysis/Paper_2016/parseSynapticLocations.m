%% Parse stimulation matrix! WTF!
close all;clear all; clc;
SN=7
stc = 5;
for fi = 1:100
    foldername = sprintf('\\\\139.91.162.90\\cluster\\stefanos\\Documents\\Glia\\OUTONLYupdatedStimGABAb01NEWBGST_Rs10c%d_SN%d_r%d',stc,SN,fi-1)
    filename = sprintf('OUTONLYupdatedStimGABAb01NEWBGST_Rs10c%d_SN%d_r%d.out',stc,SN,fi-1)
    if exist(fullfile(foldername,filename),'file')
        % Parse the old spike train data:
        fid = fopen( fullfile(foldername,filename), 'r' );
        cac = textscan( fid, '%s', 'Delimiter', '\n' );
        totalRows = size(cac{1},1);
        fclose( fid );
        
        synapticLocations = cell(700,700);
        synapticDelays = cell(700,700);
        
        
        % for total rows in data:
        fprintf('Parsing file.\n');
        ctr = 1;
        while ctr <= totalRows
            fprintf('%f\n',ctr/totalRows);
            currentRow = cac{1}{ctr,1};
            
            if strfind(currentRow,'Synaptic locations src: ')
                tmp = strsplit(currentRow, {'Synaptic locations src: ','  ,target: ','  : '});
                locmat = zeros(1,5);
                for k=1:5
                    locmat(k) = str2double(cac{1}{ctr+k,1});
                end
                synapticLocations{str2double(tmp{2})+1,str2double(tmp{3})+1} = locmat;
                ctr = ctr + 5;
            elseif strfind(currentRow,'Delays src ')
                tmp = strsplit(currentRow, {'Delays src ',', trg ',': ',' '});
                if str2double(tmp{2})+1 == str2double(tmp{3})+1
                    synapticDelays{str2double(tmp{2})+1,str2double(tmp{3})+1} = ...
                        [str2double(tmp{4})];
                else
                    synapticDelays{str2double(tmp{2})+1,str2double(tmp{3})+1} = ...
                        [str2double(tmp{4}),str2double(tmp{5}),str2double(tmp{6}),str2double(tmp{7}),str2double(tmp{8})];
                end
            end
            ctr = ctr+ 1;
        end
        fprintf('DONE parsing file.\n');
        
        % Save data:
        save(fullfile(sprintf('synapticLocDel_Rs10c%d_SN%d_r%d.mat',stc,SN,fi-1)),'synapticLocations','synapticDelays','-v7.3');
    end
end

