deg_1 = degrees(PC2PC(:,:,1));
deg_2 = degrees(PC2PC(:,:,2));
deg_3 = degrees(PC2PC(:,:,3));
deg_4 = degrees(PC2PC(:,:,4));
deg_5 = degrees(PC2PC(:,:,5));

% tmpMat = PC2PC(:,:,1);
% tmpMat = tmpMat | tmpMat';
% % [components,cliques,CC] = k_clique(3,tmpMat);
% maximalcliques_subgraphs = maximalCliques( tmpMat, 3 );

plot(histc(sort(deg_1,'descend'),0:10:200),'r');hold on;
plot(histc(sort(deg_2,'descend'),0:10:200),'g');hold on;
plot(histc(sort(deg_3,'descend'),0:10:200),'b');hold on;
plot(histc(sort(deg_4,'descend'),0:10:200),'c');hold on;
plot(histc(sort(deg_5,'descend'),0:10:200),'m');hold on;
