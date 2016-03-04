function [NCN] = m_commonNeighbors(adj)
    
% NCN = zeros(size(adj));

% Common neighbors (any) algorithm is just matrix multiplication!
% Theoretical computer science FTW!
adj = double(adj | adj');
adj = adj*adj';
adj(logical(eye(length(adj)))) = 0;
NCN = adj;


% start=1;
% for i=1:length(adj)
%     for j=start:length(adj)   
%         if (i~=j)
%             NCN(i,j) = sum(all(adj(:,[i,j]),2));
%         else
%            NCN(i,j)=0;
%         end
%     end
% %     start = start +1;
% end


end