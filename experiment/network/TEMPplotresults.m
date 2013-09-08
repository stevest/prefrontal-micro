for i=1:100
   load(sprintf('soma_%d_run_%d.txt',1,i)); 
   eval(['plot(',sprintf('soma_%d_run_%d',1,i),')']);
%    pause(0.2);
pause();
   get(gcf);cla;
end