function [NCN] = m_commonNeighbors(adj)
% Common neighbors (any) algorithm is just matrix multiplication!
% Theoretical computer science FTW!
adj = double(adj | adj');
adj = adj*adj';
adj(logical(eye(length(adj)))) = 0;
NCN = adj;

end