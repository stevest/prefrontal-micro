% Generate an instance of a random graph

function [E]=create_graph_WS(N, prob,flag)

E = zeros(N,N);

if strcmp(flag,'ind')
    % As in Song et al, 2005
    for i=1:N
        for j=1:N
            if i~=j
                randomNum = rand(1);
                if randomNum<prob
                    E(i,j)=1;
                else
                    E(i,j)=0;
                end
            else
                E(i,j)=0;
            end
        end
    end
end

end
