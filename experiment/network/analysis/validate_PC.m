function validate_PC(mypath3)

    %% mIPSC (Woo, 2007)
    EXPER=[];
    i=1
    while(exist([mypath3,sprintf('mIPSC_%d.txt',i-1)]) == 2)
        EXPER(i,:) = load([mypath3,sprintf('mIPSC_%d.txt',i-1)]);
        plot(EXPER(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        mIPSC(j) = EXPER(j,3900) - max(EXPER(j,3900:end));
    end
    mean(mIPSC) % must be ~10pA
    
    %% GABAb/GABAa ratio:
    
    % GABAa:
    EXPERa=[];
    i=1
    while(exist([mypath3,sprintf('GABAa_mV_%d.txt',i-1)]) == 2)
        EXPERa(i,:) = load([mypath3,sprintf('GABAa_mV_%d.txt',i-1)]);
        plot(EXPERa(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        GABAa(j) = EXPERa(j,3900) - min(EXPERa(j,7800));
    end
    
    % GABAb:
    EXPERb=[];
    i=1
    while(exist([mypath3,sprintf('GABAb_mV_%d.txt',i-1)]) == 2)
        EXPERb(i,:) = load([mypath3,sprintf('GABAb_mV_%d.txt',i-1)]);
        plot(EXPERb(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        GABAb(j) = EXPERb(j,3900) - min(EXPERb(j,7800));
    end
    mean(GABAb./GABAa) % must be ~0.2
    
    
    %% RS GABAa current (compared to FS GABAa current):
    EXPERfs=[];
    i=1
    while(exist([mypath3,sprintf('FSGABAa_curr_%d.txt',i-1)]) == 2)
        EXPERfs(i,:) = load([mypath3,sprintf('FSGABAa_curr_%d.txt',i-1)]);
        plot(EXPERfs(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        FSc(j) = EXPERfs(j,3900) - min(EXPERfs(j,3900:end));
    end
    
    EXPERrs=[];
    i=1
    while(exist([mypath3,sprintf('RSGABAa_curr_%d.txt',i-1)]) == 2)
        EXPERrs(i,:) = load([mypath3,sprintf('RSGABAa_curr_%d.txt',i-1)]);
        plot(EXPERrs(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        RSc(j) = EXPERrs(j,3900) - min(EXPERrs(j,3900:end));
    end
    mean(RSc./FSc) % must be ~0.1
    
    %% IS GABAa current (compared to Xenia's IS GABAa current):
    EXPERis=[];
    i=1
    while(exist([mypath3,sprintf('ISGABAa_curr_%d.txt',i-1)]) == 2)
        EXPERis(i,:) = load([mypath3,sprintf('ISGABAa_curr_%d.txt',i-1)]);
        plot(EXPERis(i,3900:end));hold on;
        i = i + 1;
    end
    
    for j=1:i-1
        ISc(j) = EXPERis(j,3900) - min(EXPERis(j,3900:end));
    end
    
    
    mean(ISc) % must be ~0.012.75 (Xenia)
end