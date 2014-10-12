function validate_CR(mypath3)

    %% NMDA/AMPA ratio:
    
    % NMDA:
    EXPERn=[];
    i=1
    while(exist([mypath3,sprintf('CRNMDA_curr_%d.txt',i-1)]) == 2)
        EXPERn(i,:) = load([mypath3,sprintf('CRNMDA_curr_%d.txt',i-1)]);
%         plot(EXPERn(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        NMDAcr(j) = EXPERn(j,3900) - max(EXPERn(j,3900:end));
    end
    NMDAcr
    
    % AMPA:
    EXPERa=[];
    i=1
    while(exist([mypath3,sprintf('CRAMPA_curr_%d.txt',i-1)]) == 2)
        EXPERa(i,:) = load([mypath3,sprintf('CRAMPA_curr_%d.txt',i-1)]);
%         plot(EXPERa(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        AMPAcr(j) = EXPERa(j,3900) - min(EXPERa(j,3900:end));
    end
    AMPAcr
    
    mean(NMDAcr./AMPAcr) % must be ~0.2
    
    %% Input resistance:
    EXPERir=[];
    i=1
    while(exist([mypath3,sprintf('CR_inputResistance_%d.txt',i-1)]) == 2)
        EXPERir(i,:) = load([mypath3,sprintf('CR_inputResistance_%d.txt',i-1)]);
%         plot(EXPERir(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        IR(j) = (EXPERir(j,3900) - EXPERir(j,4800)) / 0.1;
    end
    IR
    
    mean(IR) 

end