function exportNetworkPositions(PCpos,PVpos,outPath)
% Export cells position in 3D space for later reference.

fid = fopen([outPath,'/networkPositionsPyramidals.txt'],'w');
for i=1:size(PCpos,1)
    fprintf(fid,'%f %f %f\n',PCpos(i,1),PCpos(i,2),PCpos(i,3));
end
fclose(fid);

fid = fopen([outPath,'/networkPositionsParvalbumin.txt'],'w');
for i=1:size(PVpos,1)
    fprintf(fid,'%f %f %f\n',PVpos(i,1),PVpos(i,2),PVpos(i,3));
end
fclose(fid);

end