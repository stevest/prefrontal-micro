%% Export vertex coordinates in mathematica GraphPlot3D format:
mystr = '{';
for k=1:933
    if k<=700
        mystr = [mystr,sprintf('%d->{%f,%f,%f},',k,PCsomata(k,1),PCsomata(k,2),PCsomata(k,3))];
    else
        mystr = [mystr,sprintf('%d->{%f,%f,%f},',k,PVsomata(k-700,1),PVsomata(k-700,2),PVsomata(k-700,3))];
    end
    
end
mystr = [mystr(1:end-1),'}'];

fid = fopen('vcoords.txt','w');
fprintf(fid,mystr);
fclose(fid);

%% Export edges (connectivity) in mathematica GraphPlot3D format:
% mystr = blanks((933+935)*933+2);
mystr = cell(1,935);
mystr{1} = '{';
for k=1:933
    k
    mystr{k+1} = [mystr{k+1},'{'];
    for j=1:933
        mystr{k+1} = [mystr{k+1},sprintf('%d,',run.state_str(k,j))];
    end
    mystr{k+1} = [mystr{k+1}(1:end-1),'},'];
end
 mystr{k+1} = mystr{k+1}(1:end-1);
mystr{k+2} = [mystr{k+2},'}'];

A = strjoin(mystr,'');


fid = fopen('gadjmat.txt','w');
fprintf(fid,A);
fclose(fid);