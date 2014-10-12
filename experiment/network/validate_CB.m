function validate_CB(mypath3)

    %% NMDA/AMPA ratio:
    
    % NMDA:
    EXPERn=[];
    i=1
    while(exist([mypath3,sprintf('CBNMDA_curr_%d.txt',i-1)]) == 2)
        EXPERn(i,:) = load([mypath3,sprintf('CBNMDA_curr_%d.txt',i-1)]);
%         plot(EXPERn(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        NMDAcb(j) = EXPERn(j,3900) - max(EXPERn(j,3900:end));
    end
    NMDAcb
    
    
    % AMPA:
    EXPERa=[];
    i=1
    while(exist([mypath3,sprintf('CBAMPA_curr_%d.txt',i-1)]) == 2)
        EXPERa(i,:) = load([mypath3,sprintf('CBAMPA_curr_%d.txt',i-1)]);
%         plot(EXPERa(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        AMPAcb(j) = EXPERa(j,3900) - min(EXPERa(j,3900:end));
    end
    AMPAcb
    
    mean(NMDAcb./AMPAcb) 
    
    
    %% CR2CB GABAa validation:
     % AMPA:
    EXPERg=[];
    i=1
    while(exist([mypath3,sprintf('CRCB_GABAa_%d.txt',i-1)]) == 2)
        EXPERg(i,:) = load([mypath3,sprintf('CRCB_GABAa_%d.txt',i-1)]);
%         plot(EXPERg(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        GABAacb(j) = EXPERg(j,3900) - min(EXPERg(j,3900:end));
    end
    GABAacb
    
    mean(GABAacb) 
    
    %% Input resistance:
    EXPERir=[];
    i=1
    while(exist([mypath3,sprintf('CB_inputResistance_%d.txt',i-1)]) == 2)
        EXPERir(i,:) = load([mypath3,sprintf('CB_inputResistance_%d.txt',i-1)]);
%         plot(EXPERir(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        IR(j) = (EXPERir(j,3900) - EXPERir(j,4800)) / 0.1;
    end
    IR
    
    mean(IR) 
end