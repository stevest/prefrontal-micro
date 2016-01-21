PVvalid = [];
IC_AMP = -0.1:0.05:6 ;
for i = 1:10
    PVvalid(:,i) = load(sprintf('%svalid_PV/PV_inputResistance_%d.txt',mypath,i-1));
    VI(i) = (PVvalid(4000,i) - min(PVvalid(4001:9000,i))) / IC_AMP(i);
end

(PVvalid(4000) - min(PVvalid(4001:9000))) / 0.1

plot(PVvalid)
plot(PVvalid(4000:5000)) % nA
plot(PVvalid(10000:13000)) % nA

PVvalid = [];
for i = 1:10
    PVvalid(:,i) = load(sprintf('%svalid_PV/PV_inputResistance_%d.txt',mypath,i-1));
    AMPA(i)=(min(PVvalid(10000:11000,i)) - PVvalid(9999,i))*1000 % nA
end


spike_timing = find( ([0;diff(sign(diff(PVvalid)))<0;0] & [sign(PVvalid)==1]) );
numel(spike_timing)

(3.2e-4/6) / (7.5e-4/6)

%%
CBvalid = [];
IC_AMP = -0.1:0.05:6 ;
for i = 1:10
    CBvalid(:,i) = load(sprintf('%svalid_CB/CB_inputResistance_%d.txt',mypath,i-1));
    VI(i) = (CBvalid(4000,i) - min(CBvalid(4001:9000,i))) / IC_AMP(i);
end

(CBvalid(4000) - min(CBvalid(4001:9000))) / 0.1

plot(CBvalid)
plot(CBvalid(4000:5000)) % nA
plot(CBvalid(10000:13000)) % nA

CBvalid = [];
for i = 1:1
    CBvalid(:,i) = load(sprintf('%svalid_CB/CB_inputResistance_%d.txt',mypath,i-1));
    AMPA(i)=(min(CBvalid(10000:11000,i)) - CBvalid(9999,i))*1000 % pA
end

spike_timing = find( ([0;diff(sign(diff(CBvalid(:,i))))<0;0] & [sign(CBvalid(:,i))==1]) );
numel(spike_timing)

(3.2e-4/6) / (7.5e-4/6)

%%

CRvalid = [];
IC_AMP = -0.1:0.05:6 ;
for i = 1:10
    CRvalid(:,i) = load(sprintf('%svalid_CR/CR_inputResistance_%d.txt',mypath,i-1));
    VI(i) = (CRvalid(4000,i) - min(CRvalid(4001:9000,i))) / IC_AMP(i);
end

(CRvalid(4000) - min(CRvalid(4001:9000))) / 0.1

plot(CRvalid(:)) % nA
(max(CRvalid(10000:11000)) - CRvalid(9999))*1000 % pA

spike_timing = find( ([0;diff(sign(diff(CRvalid(:,i))))<0;0] & [sign(CRvalid(:,i))==1]) );
numel(spike_timing)

(3.2e-4/6) / (7.5e-4/6)


%%
PCvalid = [];
IC_AMP = -0.1:0.05:6 ;
for i = 1:10
    PCvalid(:,i) = load(sprintf('%svalid_PC/PC_inputResistance_%d.txt',mypath,i-1));
    VI(i) = (PCvalid(4000,i) - min(PCvalid(4001:9000,i))) / IC_AMP(i);
end

(PCvalid(4000) - min(PCvalid(4001:9000))) / 0.1

(PCvalid(4000) - PCvalid(10000)) / 0.1

plot(PCvalid(:)) % nA
(max(PCvalid(10000:11000)) - PCvalid(9999))*1000 % pA

spike_timing = find( ([0;diff(sign(diff(PCvalid(:,i))))<0;0] & [sign(PCvalid(:,i))==1]) );
numel(spike_timing)

(3.2e-4/6) / (7.5e-4/6)

