function AMPA_params()
stopIdiff=10;
targetI=45;
maxiter=100;

minampa=0.0000138;
maxampa=0.00138;
ampaweight=adaptvalue(minampa,maxampa);

Idiff=100;
iter=0;
while (ampaweight>minampa) && (ampaweight<maxampa) && (iter<maxiter)
    setampa(ampaweight);   
    sim()    
    scale=mean(getminmax());
    Idiff=scale-targetI;
    if (abs(Idiff)<stopIdiff) % finished
        break;
    end
    
    if Idiff<0 % we want a bigger value
        minampa=ampaweight;
    else    % we want a smaller value
        maxampa=ampaweight;
    end
    ampaweight=adaptvalue(minampa,maxampa);
    
    sprintf('Error: %04.15f\n',Idiff)
    iter=iter+1;
end
disp('Optimization finished')

end


function newval=adaptvalue(minval,maxval)
    newval=minval+abs(minval-maxval)/2;
end


function sim()
    [status,result] = system(['rm data/1/AMPA*']);
    [status,result] = system('../../mechanism/$(arch)/special NMDA_AMPA.hoc');
end

function setampa(ampaweight)
    [status,result] = system(['echo ' sprintf('ampaweight=%04.15f',ampaweight) ' > ampaweight.dat']);
    sprintf('new ampaweight: %04.15f\n',ampaweight)
end

function acc=getminmax()
    vals={'NMDA','AMPA'};
    s=2;
    f=dir(['data/1/' vals{s} '_*.txt']);
    fileList = {f(:).name}';
    for i=1:numel(f)
        dataSTR = sprintf('data/1/%s',fileList{i});
        NMDA=load(dataSTR);
        NMDA(1:1)=[];NMDA=NMDA*1000;
        indrange=1900:2100;
        acc(i)=abs(nanmin(NMDA(indrange))-nanmax(NMDA(indrange)));
    end
end
