function GABA_params()
% be sure to check the result variable in the sim() function to make sure
% everything is fine. 


% change this: 
paramname='gabaaweight';
hocfile='GABAA.hoc';
ofbn='gabaa'; % outputfile basename
targetI=10; % gabaa ampl bout half of gabab, which should be 20 pico amps

% this is more or less generic: 
maxiter=100;
stopIdiff=0.1*targetI; % maximal permittable error
maxgaba=0.08;
mingaba=0.000008;
if mingaba>=maxgaba
    error('maxvalue should be larger than minvalue')
end
gabaweight=adaptvalue(mingaba,maxgaba);


Idiff=9999999999;
iter=0;
while (gabaweight>mingaba) && (gabaweight<maxgaba) && (iter<maxiter)
    setparam(gabaweight,paramname);   
    sim(paramname,hocfile)    
    scale=mean(getminmax(ofbn,1));  
    sprintf('scale: %f\n',scale)
    
    Idiff=scale-targetI;
    if (abs(Idiff)<stopIdiff) % finished
        disp('Fit criterion reached')
        break;
    end    
    if Idiff<0 % we want a bigger value
        mingaba=gabaweight;
    else    % we want a smaller value
        maxgaba=gabaweight;
    end
    gabaweight=adaptvalue(mingaba,maxgaba);
    
    sprintf('Error: %f\n',Idiff)
    iter=iter+1;
end
disp('Optimization finished')
sprintf('%d iterations',iter)
end


function newval=adaptvalue(minval,maxval)
    newval=mean([maxval minval]);
end


function sim(paramname,hocfile)
    [status,~] = system(['rm data/1/' paramname '*']);
    [status,result] = system(['../../mechanism/$(arch)/special ' hocfile]);
    if status~=0
        error('system call returns error')        
    end    
end

function setparam(gabaweight,paramname)
    [status,result] = system(['echo ' sprintf([paramname '=%04.15f'],gabaweight) ' > ' paramname '.dat']);
    sprintf(['new ' paramname ' : %f\n'],gabaweight)
end

function acc=getminmax(ofbn,miden)    
    clear acc    
    f=dir(['data/1/' ofbn '_*.txt']);
    fileList = {f(:).name}';
    for i=1:numel(f)
        dataSTR = sprintf('data/1/%s',fileList{i});
        gaba=load(dataSTR);
        %gaba(1:1)=[];
        gaba=gaba*1000;  % converting to picoA
        indrange=2000:2100;
        if nargin>1 && miden==1
            acc(i)=abs(nanmax(gaba(indrange)));
        else
            acc(i)=abs(nanmin(gaba(indrange))-nanmax(gaba(indrange)));
        end
    end
end
