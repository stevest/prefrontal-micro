cd /home/gkastel/Desktop/Stefanos/PFC_C3_5_alltogether/experiment/network/VALIDATION_NMDA_ATTEN/2/ 

for i=0:15
   name = sprintf('NMDA_ATTEN_%d.txt',i);
   load(name) ;
   eval(['plot(',name(1:end-4),'(4900:end))']);hold on;
   eval(['plot(',name(1:end-4),'(4900:end))']);hold on;
%    pause
end

%%
cd /home/gkastel/Desktop/Stefanos/PFC_C3_5_alltogether/experiment/network/VALIDATION_NMDA/Spikes/

for i=0:11
   name = sprintf('NMDA_%d.txt',i);
   load(name) ;
   eval(['plot(',name(1:end-4),'(4900:end))']);hold on;
%    eval(['plot(',name(1:end-4),'(4900:end))']);hold on;
%    pause
end