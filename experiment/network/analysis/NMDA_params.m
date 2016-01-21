function NMDA_params()
stopIdiff=10;
targetI=50;
maxiter=100;

minnmda=0.00001;
maxnmda=0.001; 
nmdaweight=adaptvalue(minnmda,maxnmda);

Idiff=100;
iter=0;
while (nmdaweight>minnmda) && (nmdaweight<maxnmda) && (iter<maxiter)
    setnmda(nmdaweight);   
    sim()    
    scale=mean(getminmax());
    sprintf('scale: %04.15f\n',scale)
    
    Idiff=scale-targetI;
    if (abs(Idiff)<stopIdiff) % finished
        break;
    end    
    if Idiff<0 % we want a bigger value
        minnmda=nmdaweight;
    else    % we want a smaller value
        maxnmda=nmdaweight;
    end
    nmdaweight=adaptvalue(minnmda,maxnmda);
    
    sprintf('Error: %04.15f\n',Idiff)
    iter=iter+1;
end
disp('Optimization finished')

end


function newval=adaptvalue(minval,maxval)
    newval=minval+abs(minval-maxval)/2;
end


function sim()
    [status,result] = system(['rm data/1/NMDA*']);
    [status,result] = system('../../mechanism/$(arch)/special NMDA_AMPA.hoc');
end

function setnmda(nmdaweight)
    [status,result] = system(['echo ' sprintf('nmdaweight=%04.15f',nmdaweight) ' > nmdaweight.dat']);
    sprintf('new nmdaweight: %04.15f\n',nmdaweight)
end

function acc=getminmax()
    clear acc
    vals={'NMDA','AMPA'};s=1;
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
