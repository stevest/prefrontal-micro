% nPC = 75 ;
allCombs3 = combntns(1:nPC,3);
motifs = cell(1,16);
for i=1:length(allCombs3)
    
    mt = AllConnMat(allCombs3(i,:),allCombs3(i,:));
    
    if sum(sum(mt .* ~eye(3))) == 0
        motifs(end,1) = {allCombs3(i,:)};
    elseif sum(sum(mt .* ~eye(3))) == 1
        motifs(end,2) = {allCombs3(i,:)};
    elseif sum(sum(mt .* ~eye(3))) == 5
        motifs(end+1,15) = {allCombs3(i,:)};
    elseif sum(sum(mt .* ~eye(3))) == 6
        motifs(end+1,16) = {allCombs3(i,:)};
    elseif sum(sum(mt .* ~eye(3))) == 2
          if sum(sum( mt.*mt' .* ~eye(3) )) == 2
              motifs(end+1,3) = {allCombs3(i,:)}; %3
          elseif any(sum(mt .* ~eye(3),2)==2)
              motifs(end+1,3) = {allCombs3(i,:)}; %3
            
    elseif sum(sum(mt .* ~eye(3))) == 3
        % sub categories:
    elseif sum(sum(mt .* ~eye(3))) == 4
        % subcategories:
    end
end
