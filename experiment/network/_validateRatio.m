function [SG,current] = validateRatio( MODE, PYID, RUN_NO, MINVAR, MAXVAR, MAXERR, TARGET, MAXITER)
stopIdiff=MAXERR; %10;
targetI=TARGET; %45;
maxiter=MAXITER; %100;

minSG=MINVAR; %0.0000138;
maxSG=MAXVAR; %0.1;
SG=adaptvalue(minSG,maxSG);

Idiff=100;
iter=0;
while (SG>minSG) && (SG<maxSG) && (iter<maxiter)
%     setG(SG);   
    sim(SG, MODE, PYID, RUN_NO)    
    current=mean(getminmax(MODE,PYID));
    Idiff=current-targetI;
    if (abs(Idiff)<stopIdiff) % finished
        sprintf('current is: %04.15f\n',current)
        sprintf('Error is: %04.15f\n',Idiff)
        break;
    end
    
    if Idiff<0 % we want a bigger value
        minSG=SG;
    else    % we want a smaller value
        maxSG=SG;
    end
    SG=adaptvalue(minSG,maxSG);
    
    sprintf('Error: %04.15f\n',Idiff)
    iter=iter+1;
end
disp('Optimization finished')

end


function newval=adaptvalue(minval,maxval)
    newval=minval+abs(minval-maxval)/2;
end


function sim(OPTVAR, MODE, PYID, RUN_NO)
    [status,result] = system(['rm data/1/AMPA*']);
    [status,result] = system(['../../mechanism/$(arch)/special -nobanner ',sprintf('-c "PYID=%d" -c "MODE=%d" -c "RUN_NO=%d" -c "OPTVAR=%f"',PYID,MODE,RUN_NO,OPTVAR),' NMDA_AMPA.hoc']);
end

% function setG(weight)
%     [status,result] = system(['echo ' sprintf('ampaweight=%04.15f',weight) ' > ampaweight.dat']);
%     sprintf('new ampaweight: %04.15f\n',weight)
% end

function acc=getminmax(MODE, PYID)
    indrange=200:600; % HARDCODED RANDE!
    vals={'NMDA','AMPA'};
    s=MODE;
    f=dir(sprintf('data/NMDA_AMPA_ratio/%d/%s_%d_*.txt',MODE, vals{s}, PYID));
    fileList = {f(:).name}';
    for i=1:numel(f)
        dataSTR = sprintf('data/NMDA_AMPA_ratio/%d/%s',MODE,fileList{i});
        NMDA=load(dataSTR);
        if(isempty(NMDA))
            continue;
        end
%         NMDA(1:1)=[];NMDA=NMDA*1000;
        acc(i)=abs(nanmin(NMDA(indrange))-nanmax(NMDA(indrange)));
    end
end
