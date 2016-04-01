%% Parse network parameters! WTF!

% Parse the old spike train data:
fid = fopen( sprintf('X:\\Documents\\GitHub\\prefrontal-micro\\experiment\\network\\experiment_12\\_importNetworkParametersSTR.hoc'), 'r' );
cac = textscan( fid, '%s', 'Delimiter', '\n' );
% size(cac{1},1) must equals the # of rows in your data file.
totalRows = size(cac{1},1);
fclose( fid );

C = zeros(933,933);
W = zeros(933,933);

% Regular way (can read anything) for total rows in data:
fprintf('Parsing file.\n');
for k=1:totalRows
    fprintf('%f\n',k/totalRows);
    currentRow = cac{1}{k,1};
    
    if strfind(currentRow,'connMatrix.x[')
        tmp = strsplit(currentRow, ' = ');
        tmploc = strsplit(tmp{1}, '[');
        m = str2double(tmploc{2}(1:end-1))+1;
        n = str2double(tmploc{3}(1:end-1))+1;
        C(m,n) = logical(str2double(tmp{2}));
    elseif strfind(currentRow,'weightsMatrix.x[')
        tmp = strsplit(currentRow, ' = ');
        tmploc = strsplit(tmp{1}, '[');
        m = str2double(tmploc{2}(1:end-1))+1;
        n = str2double(tmploc{3}(1:end-1))+1;
        W(m,n) = str2double(tmp{2});
    end

end
fprintf('DONE parsing file.\n');

%%
orgC = C;
orgW = W;
%% Parfor: requires the data to be written in particular order!
N = 700;
begining = 0;
for k=1:totalRows
    currentRow = cac{1}{k,1};
    if strfind(currentRow,'connMatrix.x[')
        begining = k;
        break;
    end
end
ending = begining+N^2-1;

 % slice input array:
fprintf('Parsing connections.\n');
parfor k=begining:ending
    currentRow = cac{1}{k,1};
    tmp = strsplit(currentRow, '=');
    m = mod(k-begining,N);
    n = floor((k-begining)/N)
    C(m,n) = logical(str2double(tmp{2}));
    
end
fprintf('DONE parsing file.\n');
