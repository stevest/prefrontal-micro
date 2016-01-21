function exportNetworkCluster(PCclusters,outPath)
% Export clustering : values are cell ID exemplars corresponding to PC
% positions.

fid = fopen([outPath,'/networkPyramidalClustersStructured.txt'],'w');
for i=1:size(PCclusters,1)
    fprintf(fid,'%d\n',PCclusters(i,1));
end
fclose(fid);

fid = fopen([outPath,'/networkPyramidalClustersRandom.txt'],'w');
for i=1:size(PCclusters,1)
    fprintf(fid,'%d\n',PCclusters(i,2));
end
fclose(fid);

end