function validate_PV(mypath3)

    %% NMDA/AMPA ratio:
    
    % NMDA:
    EXPERn=[];
    i=1
    while(exist([mypath3,sprintf('PVNMDA_curr_%d.txt',i-1)]) == 2)
        EXPERn(i,:) = load([mypath3,sprintf('PVNMDA_curr_%d.txt',i-1)]);
%         plot(EXPERn(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        NMDApv(j) = EXPERn(j,3900) - max(EXPERn(j,3900:end));
    end
    NMDApv
    
    
    
    % AMPA:
    EXPERa=[];
    i=1
    while(exist([mypath3,sprintf('PVAMPA_curr_%d.txt',i-1)]) == 2)
        EXPERa(i,:) = load([mypath3,sprintf('PVAMPA_curr_%d.txt',i-1)]);
%         plot(EXPERa(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        AMPApv(j) = EXPERa(j,3900) - min(EXPERa(j,3900:end));
    end
    AMPApv
    
    
    mean(NMDApv./AMPApv) % must be ~0.2
    
    

end