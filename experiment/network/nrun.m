classdef nrun < handle
    %NRUN class to run NEURON remotely...
    properties
        id
        sn
        path
        pathToHere
        nruns
        state
        gparams
        stimstart
        stimend
        stimdur
        tstop
        dt
        nPC
        nPV
        nCB
        nCR
        nAll
        PCsomata
        PVsomata
        CBsomata
        CRsomata
        connBinsPC2PC
        connProbsPC2PC
        connBinsPV2PC
        connProbsPV2PC
        connBinsCB2PC
        connProbsCB2PC
        recipBinsPC2PC
        recipProbsPC2PC
        NconnBins
        NincomingProbs
        NoutgoingProbs
        distPC2PC
        distPV2PC
        distCB2PC
        ConnMatPC2PV
        ConnMatPV2PC
        ConnMatCB2PC
        connmatrix
        
        
        state_str
        state_rnd
        weights_str
        weights_rnd
        
        mergedCN_str
        mergedCN_rnd
        NC_str
        NC_rnd
        cellsPerCluster_str
        cellsPerCluster_rnd
        stimCellsPerCluster
        
        labels_str
        labels_rnd
        
        StimMat_str
        StimMat_rnd
        
        SpikesStimDend
        SpikesStimApic
        
        SpikesBgBasal
        SpikesBgProximal
        SpikesBgApical
        SpikesBgPV
        SpikesBgCB
        SpikesBgCR
        
       
        SLASH
        ISBINARY
        CLUSTBIAS
    end
    
    methods %(Access = public)
        function obj = nrun(id,npyrs,serialno,state,exprun,stop)
            % Constructor:
            if(strfind(computer,'PC'))
                obj.SLASH = '\';
            else
                obj.SLASH = '/';
            end
            % Use binary format when exporting files for speed of transfer
            % (achieves x500 less filesize!!)
            obj.ISBINARY = 1;
            
%             rng(id,'twister');
            % seed the random generator with the id of the experiment:
            % Affects stimulation, background etc..
            RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',id));
            obj.id = id;
            obj.sn = serialno;
            obj.dt = 10;
            obj.state = state;
            obj.path = sprintf('experiment_%d',obj.id);
            obj.tstop = stop;
            obj.nruns=exprun;
            obj.stimstart = 500;
            obj.stimend = 1500;
            obj.stimdur = obj.stimend - obj.stimstart;
            obj.nPC = npyrs;

        end
        function init(obj,varargin)
            system(sprintf('mkdir experiment_%d',obj.id));
            
            if nargin > 1
                states = varargin{1};
                %load precomputed params from states file:
                obj.gparams = states.gparams(obj.state);
                % treat cell no like local parameter (overrided from the
                % pre-generated state):
                obj.nPC = obj.gparams.nPC;
                obj.nPV = obj.gparams.nPV;
                obj.nCB = obj.gparams.nCB;
                obj.nCR = obj.gparams.nCR;
                obj.nAll = obj.gparams.nAll;
                % load states:
                obj.state_str = states.state_str(:,:,obj.state);
                obj.state_rnd = states.state_rnd(:,:,obj.state);
                
                obj.labels_str = obj.gparams.labels_str;
                obj.labels_rnd = obj.gparams.labels_rnd;
                
                obj.cellsPerCluster_str = obj.gparams.cellsPerCluster_str ;
                obj.cellsPerCluster_rnd = obj.gparams.cellsPerCluster_rnd ;
                obj.stimCellsPerCluster = obj.gparams.stimCellsPerCluster ;
                
                obj.NC_str = obj.gparams.NC_str;
                obj.NC_rnd = obj.gparams.NC_rnd;
                
                
            else
%             obj.gparams.connBinsPC2PC = 0:30:500;
%             obj.gparams.connProbsPC2PC = 0.25 .* exp(-0.006 * obj.gparams.connBinsPC2PC);
%             obj.gparams.connBinsPV2PC = 0:20:500;
%             obj.gparams.connProbsPV2PC = linspace(1,0,length(obj.gparams.connBinsPV2PC)) * 0.3; % Takeshi Otsuka, 2009:
%             obj.gparams.connBinsCB2PC = 0:20:500;
%             obj.gparams.connProbsCB2PC = linspace(1,0,length(obj.gparams.connBinsCB2PC)) ;
%             obj.gparams.recipBinsPC2PC = 0:30:500;
%             obj.gparams.recipProbsPC2PC = 0.12 .* exp(-0.006 * obj.gparams.recipBinsPC2PC);
%             % Update probabilities from Perin et al. (Figure S4)
%             obj.gparams.NconnBins = [0,1,2,3,4];
%             obj.gparams.NincomingProbs = [0.1, 0.2, 0.25, 0.4, 0.52] ;
%             obj.gparams.NoutgoingProbs = [0.12, 0.25, 0.24, 0.2, 0.23] ;
%             
%             %  ---- Initialize Pyramidal cells ----
%             obj.gparams.PCsomata = nrun.CreateRandomNetwork(obj.nPC, 90, 3); % was 50, to much...
%             obj.gparams.distPC2PC = generateDistanceMat(obj.gparams.PCsomata', 0);
%             %  ---- Rest of the cells ----
%             obj.nPV = round(obj.nPC*13/75);
%             obj.nCB = round(obj.nPC*6/75);
%             obj.nCR = round(obj.nPC*6/75);
%             if(obj.nPV==0)% NEURON ISSUE:
%                 obj.nPV = 1;
%             end
%             if(obj.nCB == 0)
%                 obj.nCB = 1;
%             end
%             if(obj.nCR == 0)
%                 obj.nCR=1;
%             end
%             obj.nAll = sum([obj.nPC,obj.nPV,obj.nCB,obj.nCR]);
%             % PV somata: about 40 per mm2 for 50um depth of slice as in:
%             %      Ohira, K., Takeuchi, R., Iwanaga, T., & Miyakawa, T. (2013).
%             %      Chronic fluoxetine treatment reduces parvalbumin expression and
%             %      perineuronal nets in gamma-aminobutyric acidergic interneurons of
%             %      the frontal cortex in adult mice. Molecular Brain, 6(1), 43.
%             %      doi:10.1186/1756-6606-6-43
%             
%             % Ara exw peripou 12.5 PV se enan kybo 250x250x250 um
%             obj.gparams.PVsomata = CreateCubeNetworkPV(90, obj.nPV); % 226.6 per mm squared (?!?) paper???
%             obj.gparams.CBsomata = CreateCubeNetworkPV(90, obj.nCB);
%             obj.gparams.CRsomata = CreateCubeNetworkPV(0, obj.nCR);
%             
%             % Distance of each interneuron from Pyramidal and connect based on
%             % probability from Yuste 2011:
%             
%             % gia to connectivity twn PV 2 PC  o Packer et al., 2011:
%             % elegxei gia ena region gyrw apo to ka8e PC. Epomenws ta histograms einai
%             % swsta gia ta PV2PC (ka8ws briskontai sto idio layer). Oso gia ta
%             % CB(SOM)2PC pou briskontai se diaforetika layers mporoume na eikasoume
%             % oti:
%             % * Apo to sxhma twn SOM (Packer et al., 2013), oso metakineisai sto
%             % transverse layer oi tuft dendrites enws PC 8a exoun tis idies pi8anotites
%             % gia overlap opos exoun kai ta PC pou briskontai sto idio layer. Epomenws
%             % metraw san distance CB2PC mono to distance sto transverse plane kai
%             % xrisimopoiw tis pi8anotites pou dinei o Packer et al., 2013, Fig4D
%             obj.gparams.distPV2PC = distancePV2PC(obj.gparams.PVsomata,obj.gparams.PCsomata);
%             obj.gparams.distCB2PC = distancePV2PC(obj.gparams.CBsomata,obj.gparams.PCsomata);
%             
%             obj.gparams.ConnMatPV2PC = connectPV2PC(obj.gparams.distPV2PC,obj.gparams.connBinsPV2PC,obj.gparams.connProbsPV2PC);
%             obj.gparams.ConnMatCB2PC = connectCB2PC(obj.gparams.distCB2PC,obj.gparams.connBinsCB2PC,obj.gparams.connProbsCB2PC);
%             obj.gparams.CRsomata = CreateCubeNetworkPV(0, obj.nCR);
%             obj.gparams.ConnMatPC2PV = obj.connectPC2PV(obj.gparams.distPV2PC')  ;
%             obj.constraintPVreciprocity();
%             
%             % Set indices for ease of mind:
%             pc = size(obj.gparams.PCsomata,1);
%             pv = size(obj.gparams.PCsomata,1) + size(obj.gparams.PVsomata,1);
%             cb = size(obj.gparams.PCsomata,1) + size(obj.gparams.PVsomata,1) + size(obj.gparams.CBsomata,1);
%             cr = size(obj.gparams.PCsomata,1) + size(obj.gparams.PVsomata,1) + size(obj.gparams.CBsomata,1) + size(obj.gparams.CRsomata,1);
%             
%             % Populate the final all-to-all connectivity matrix:
%             obj.connmatrix = zeros(obj.nAll);
%             
%             %         % Pyramidals to all types of interneurons:
%             %         AllConnMat(1:pc,pc+1:end) = 1; % Connected to all
%             % PCs to PVs
%             obj.connmatrix(1:pc,pc+1:pv) = obj.ConnMatPC2PV;
%             %     AllConnMat(1:pc,pv+1:end) = 1; % Connected to all interneurons: too much. maybe less connectivity?
%             obj.connmatrix(1:pc,pv+1:end) = rand(obj.nPC, obj.nAll-pv) > 0.9; % Connected to all interneurons: too much. maybe less connectivity?
%             % PVs connect to all other PVs + autapses (need to gap junctions?)
%             obj.connmatrix(pc+1:pv,pc+1:pv) = 1;
%             % PVs to PCs based on above connectivity (Yuste 2011):
%             obj.connmatrix(pc+1:pv,1:pc) = obj.ConnMatPV2PC;
%             % CB only connect to PC (Xenia) na to psa3w...
%             obj.connmatrix(pv+1:cb,1:pc) = obj.ConnMatCB2PC;
%             % CR only connect to PC and CB (Xenia) na to psa3w...
%             obj.connmatrix(cb+1:cr,1:pc) = 1;
%             obj.connmatrix(cb+1:cr,pv+1:cb) = 1;
%             
% %             % Load precomputed network states (both random and structured):
% %             load('states.mat');
%             
%             obj.state_str = obj.connmatrix;
%             obj.state_str(1:obj.nPC,1:obj.nPC) = states.PC2PC_str(:,:,obj.state);
%             obj.state_rnd = obj.connmatrix;
%             obj.state_rnd(1:obj.nPC,1:obj.nPC) = states.PC2PC_rnd(:,:,obj.state);
%             
%             % Find nearest neighbors of Pyramidals only
%             [CNi_str,CNo_str] = m_commonNeighbors(obj.state_str(1:obj.nPC,1:obj.nPC));
%             [CNi_rnd,CNo_rnd] = m_commonNeighbors(obj.state_rnd(1:obj.nPC,1:obj.nPC));
%             obj.mergedCN_str = [CNi_str + CNo_str] .* obj.state_str(1:obj.nPC,1:obj.nPC);
%             obj.mergedCN_rnd = [CNi_rnd + CNo_rnd] .* obj.state_rnd(1:obj.nPC,1:obj.nPC);
%             mergedCN_strT = obj.mergedCN_str';
%             obj.mergedCN_str(logical(tril(ones(size(obj.mergedCN_str)),-1))) = mergedCN_strT(logical(tril(ones(size(obj.mergedCN_str)),-1)));
%             mergedCN_rndT = obj.mergedCN_rnd';
%             obj.mergedCN_rnd(logical(tril(ones(size(obj.mergedCN_rnd)),-1))) = mergedCN_strT(logical(tril(ones(size(obj.mergedCN_rnd)),-1)));
% 
%             % Affinity Propagation:
%             obj.cellsPerCluster_str = [];
%             obj.cellsPerCluster_rnd = [];
%             % Performing affinity propagation inside class depricated
% %             performAffinityPropagation(obj);
%             Sid = 1;
%             obj.NC_rnd(Sid) = length(unique(states.Labels_rnd(:,obj.state)));
%             [~,~,obj.labels_rnd(:,Sid)] = unique(states.Labels_rnd(:,obj.state));
%             obj.NC_str(Sid) = length(unique(states.Labels_str(:,obj.state)));
%             [~,~,obj.labels_str(:,Sid)] = unique(states.Labels_str(:,obj.state));
%             %how many cells in each structured cluster?
%             obj.cellsPerCluster_str = histc(obj.labels_str(:,Sid),1:obj.NC_str(Sid))';
%             %how many cells in each random cluster?
%             obj.cellsPerCluster_rnd = histc(obj.labels_rnd(:,Sid),1:obj.NC_rnd(Sid))';
%             obj.stimCellsPerCluster = min(min(obj.cellsPerCluster_str),min(obj.cellsPerCluster_rnd));
            end
            % Generate weights matrices:
            % Kayaguchi : more potent synapses in reciprocal pairs
            % Perin: more potent synapses with more clustering coeff:
            % We try to combine them
            obj.weights_str = ones(obj.nAll);
            obj.weights_rnd = ones(obj.nAll);
            
%             [~,~,cc] = clust_coeff(obj.state_str(1:obj.nPC,1:obj.nPC));
%             obj.weights_str(repmat(nrun.expamp(cc)',[],obj.nPC))
            for cl = 1:obj.NC_str
                idx = find(obj.labels_str(:,1) == cl);
                obj.weights_str(idx,idx) = 1.3; 
            end
            for cl = 1:obj.NC_rnd
                idx = find(obj.labels_rnd(:,1) == cl);
                obj.weights_rnd(idx,idx) = 1.3;
            end
            
            
            % Create Stimulation Matrices:
            obj.StimMat_str = zeros(obj.stimCellsPerCluster,obj.NC_str);
            for t = 1:obj.NC_str
                [r,~,~] = find(obj.labels_str == t);
                tmpr = randperm(length(r));
                obj.StimMat_str(:,t) = r(tmpr(1:obj.stimCellsPerCluster));
            end
            
            obj.StimMat_rnd = zeros(obj.stimCellsPerCluster,obj.NC_rnd);
            for t = 1:obj.NC_rnd
                [r,~,~] = find(obj.labels_rnd == t);
                tmpr = randperm(length(r));
                obj.StimMat_rnd(:,t) = r(tmpr(1:obj.stimCellsPerCluster));
            end

            % Generate input spiketrains:
            obj.generateStimulation();
            obj.generateBackgroundActivity();
            
            % Export Parameters to NEURON:
%             obj.exportParams2Neuron();
                        
        end
        function push(obj)
            % Loacl copy (duh):
            copyfile(sprintf('%s%s%s/',obj.pathToHere,obj.SLASH,obj.path),sprintf('C:\\cygwin64\\home\\steve\\%s',obj.path));
            % Upload experiment folder to biosrv cluster:
            system(sprintf('C:\\cygwin64\\bin\\bash.exe -l -c "/home/steve/jarvis2cluster.exp %s"',obj.path));
        end
        function pull(obj)
             % Download experiment folder:
            system(sprintf('C:\\cygwin64\\bin\\bash.exe -l -c "/home/steve/cluster2jarvis.exp /home/stefanos/Desktop/prefrontal-micro/experiment/network/%s %s"',obj.path,'/home/steve/'));
            
            % Loacl copy - again (duh):
            copyfile(sprintf('C:\\cygwin64\\home\\steve\\%s',obj.path),sprintf('%s%s%s/',obj.pathToHere,obj.SLASH,obj.path));
        end
        function run(obj,cid,EXPERIMENT,SIMPLIFIED,PARALLEL)
            
%             % Push data to cluster:
%             obj.push();
            
            % Run the experiment remotely:
            VCLAMP = 0; % Depricated
            if (SIMPLIFIED)
                neuroncmd = sprintf(['mpirun -v -n 24 '...
                    '../../mechanism_simple/x86_64/special -mpi -nobanner '...
                    '-c "PARALLEL=%d" '...
                    '-c "SIMPLIFIED=%d" '...
                    '-c "CLUSTER_ID=%d" '...
                    '-c "EXPERIMENT=%d" '...
                    '-c "ST=%d" '...
                    '-c "ID=%d" '...
                    '-c "SN=%d" '...
                    '-c "VCLAMP=%f" '...
                    '-c "ISBINARY=%d" '...
                    '-c "CLUSTBIAS=%f" '...
                    'final.hoc'], ...
                    PARALLEL,SIMPLIFIED,cid-1,EXPERIMENT,obj.state,obj.id,obj.sn,VCLAMP,obj.ISBINARY,obj.CLUSTBIAS);
            else
                neuroncmd = sprintf(['mpirun -n 24 -mca btl ^openib '...
                    '../../mechanism_complex/x86_64/special -mpi -nobanner '...
                    '-c "PARALLEL=%d" '...
                    '-c "SIMPLIFIED=%d" '...
                    '-c "CLUSTER_ID=%d" '...
                    '-c "EXPERIMENT=%d" '...
                    '-c "ST=%d" '...
                    '-c "ID=%d" '...
                    '-c "SN=%d" '...
                    '-c "VCLAMP=%f" '...
                    '-c "ISBINARY=%d" '...
                    '-c "CLUSTBIAS=%f" '...
                    'final.hoc'], ...
                    PARALLEL,SIMPLIFIED,cid-1,EXPERIMENT,obj.state,obj.id,obj.sn,VCLAMP,obj.ISBINARY,obj.CLUSTBIAS);
            end
            
            % run neuron:
            tic;
            system(sprintf('C:\\cygwin64\\bin\\bash.exe -l -c "/home/steve/nexecute.exp %s"',neuroncmd))
%             system(sprintf('C:\\cygwin64\\bin\\bash.exe -l -c "/home/steve/nexecute.exp %s %s %s %s %s %s %s"',))
            runtime = toc
            
%            % pull results from cluster:
%            obj.pull();
           
        end
        
        
        function constraintPVreciprocity(obj)
            % Takeshi Otsuka, 2009:
            % Test if the PC2PV reciprocals are 20% of PC2PV pairs:
            ctr = 10000;
            while ctr && (1.5 < abs(((sum(sum((obj.ConnMatPC2PV & obj.ConnMatPV2PC')))*100) / numel(obj.ConnMatPC2PV)) - 20))
                % if reciprocals are more/less, move connections to reach 20%:
                foundIt = 0;
                while ~foundIt
                    ip = ceil(rand(1)*size(obj.ConnMatPC2PV,1));
                    jp = ceil(rand(1)*size(obj.ConnMatPC2PV,2));
                    if obj.ConnMatPC2PV(ip,jp) && ~obj.ConnMatPV2PC(jp,ip)
                        foundIt = 1;
                    end
                end
                foundIt = 0;
                while ~foundIt
                    ii = ceil(rand(1)*size(obj.ConnMatPC2PV,1));
                    ji = ceil(rand(1)*size(obj.ConnMatPC2PV,2));
                    if ~obj.ConnMatPC2PV(ii,ji) && obj.ConnMatPV2PC(ji,ii)
                        foundIt = 1;
                    end
                end
                % move the projection of the PC to form reciprocal with the PV:
                obj.ConnMatPC2PV(ip,jp) = 0;
                obj.ConnMatPC2PV(ii,ji) = 1;
                ctr = ctr - 1 ;
            end
        end
        function performAffinityPropagation(obj)
            %perform affinity propagation algorithm:
            Sid=1; % depricated variable from another clustering algo..
            
            % try to force different No of cluster to get as many clusters as you
            % can from the populations with high min cells per cluster:
            srtClNo = 10;
            setFlg = 1;
            while setFlg
                for i=1:10
                    [idx_str,~,~,~,~]=apclusterK(obj.mergedCN_str,srtClNo);
                    obj.NC_str(Sid) = length(unique(idx_str));
                    [~,~,obj.labels_str(:,Sid)] = unique(idx_str);
                    if(min(histc(obj.labels_str(:,Sid),1:obj.NC_str(Sid))') > 6)
                        setFlg = 0;
                        break;
                    end
                    
                end
                srtClNo = srtClNo - 1;
            end
            
            
            % try to force different nNo of cluster to get as many clusters as you
            % can from the populations with high min cells per cluster:
            srtClNo = 10;
            setFlg = 1;
            while setFlg
                for i=1:10
                    [idx_rnd,~,~,~,~]=apclusterK(obj.mergedCN_rnd,srtClNo);
                    obj.NC_rnd(Sid) = length(unique(idx_rnd));
                    [~,~,obj.labels_rnd(:,Sid)] = unique(idx_rnd);
                    if(min(histc(obj.labels_rnd(:,Sid),1:obj.NC_rnd(Sid))') > 6)
                        setFlg = 0;
                        break;
                    end
                    
                end
                srtClNo = srtClNo - 1;
            end
            
        end
        function generateStimulation(obj)
            % Generate stimulation data:
            nIncoming=[200,200];
            randShift = 0.1;
            nStim = 30; % how many Hz of stimulation?
            obj.SpikesStimDend = {};
            obj.SpikesStimApic = {};
            for c=1:obj.nPC
                i=1;
                while i <= nIncoming(1)
                    obj.SpikesStimDend(c,i)={linspace(0,obj.stimdur,nStim)};
                    obj.SpikesStimDend{c,i} = round(obj.SpikesStimDend{c,i} + ((rand(1,nStim)-0.5)*randShift));
                    obj.SpikesStimDend{c,i} = obj.SpikesStimDend{c,i}(obj.SpikesStimDend{c,i}>0);
                    obj.SpikesStimDend{c,i} = obj.SpikesStimDend{c,i}(obj.SpikesStimDend{c,i}<obj.stimdur);
                    obj.SpikesStimDend{c,i} = sort((obj.SpikesStimDend{c,i}),'ascend');
                    if(length(unique(obj.SpikesStimDend{c,i})) ~= length(obj.SpikesStimDend{c,i}) )
                        continue;
                    end
                    i = i +1;
                end
                i=1;
                while i <= nIncoming(2)
                    obj.SpikesStimApic(c,i)={linspace(0,obj.stimdur,nStim)};
                    obj.SpikesStimApic{c,i} = round(obj.SpikesStimApic{c,i} + ((rand(1,nStim)-0.5)*randShift));
                    obj.SpikesStimApic{c,i} = obj.SpikesStimApic{c,i}(obj.SpikesStimApic{c,i}>0);
                    obj.SpikesStimApic{c,i} = obj.SpikesStimApic{c,i}(obj.SpikesStimApic{c,i}<obj.stimdur);
                    obj.SpikesStimApic{c,i} = sort((obj.SpikesStimApic{c,i}),'ascend');
                    if(length(unique(obj.SpikesStimApic{c,i})) ~= length(obj.SpikesStimApic{c,i}) )
                        continue;
                    end
                    i = i +1;
                end
            end
        end
        function generateBackgroundActivity(obj)
            % Background activity must result in ~1.2Hz in the post-synaptic cell as in
            % Carl C.H.Petersen Neuron, 2013
            
            num_exc = 160;%100;
            num_inha = 75;%75
            num_inhb = 30;%30
            num_inhc = 30;%30
            
            obj.SpikesBgBasal = {};
            obj.SpikesBgProximal = {};
            obj.SpikesBgApical = {};
            obj.SpikesBgPV = {};
            obj.SpikesBgCB = {};
            obj.SpikesBgCR = {};
            
            % for many_times=1:runs+1
            %     rng(many_times)
            stimPoisson = find(poissrnd(0.001,obj.tstop,1))  ;%0.006 mporei na einai poly gia to structured
            %     stimPoisson = sort([stimPoisson,rand(tstop]);
            % stimPoisson = [200:400:3000]' ;
            TID=100;
            
            % POUSTIA terastiwn diastasewn:
            % Apo to Poisson trace bgazw ta events pou einai pio konta apo time TID
            stimPoisson(diff(stimPoisson)<(TID*2)) = [];
            
            % stimPoisson = [TID:TID:tstop]' ;
            
            % Mipws oi synapseis pou dinw gia background na einai anti-clustered?
            for c=1:obj.nPC
                poil = length(stimPoisson) ;
                for i = 1:num_exc
                    tmp = sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(round(obj.tstop/10000),1)*obj.tstop],'ascend') ;
                    obj.SpikesBgBasal(obj.nruns,c,i) = {unique(tmp(tmp<obj.tstop))};
                    obj.SpikesBgBasal{obj.nruns,c,i}(obj.SpikesBgBasal{obj.nruns,c,i}<0) = [];
                end
                for i = 1:num_exc
                    tmp =  sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(round(obj.tstop/10000),1)*obj.tstop],'ascend') ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
                    obj.SpikesBgProximal(obj.nruns,c,i) = {unique(tmp(tmp<obj.tstop))};
                    obj.SpikesBgProximal{obj.nruns,c,i}(obj.SpikesBgProximal{obj.nruns,c,i}<0) = [];
                end
                for i = 1:num_exc
                    tmp =  sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(round(obj.tstop/10000),1)*obj.tstop],'ascend') ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
                    obj.SpikesBgApical(obj.nruns,c,i) = {unique(tmp(tmp<obj.tstop))};
                    obj.SpikesBgApical{obj.nruns,c,i}(obj.SpikesBgApical{obj.nruns,c,i}<0) = [];
                end
            end
            
            for c = 1:obj.nPV
                for i = 1:num_inha
                    tmp =  sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(round(obj.tstop/10000),1)*obj.tstop],'ascend') ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
                    obj.SpikesBgPV(obj.nruns,c,i) =  {unique(tmp(tmp<obj.tstop))};
                    obj.SpikesBgPV{obj.nruns,c,i}(obj.SpikesBgPV{obj.nruns,c,i}<0) = [];
                end
            end
            
            for c = 1:obj.nCB
                for i = 1:num_inhb
                    tmp =  sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(round(obj.tstop/10000),1)*obj.tstop],'ascend') ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
                    obj.SpikesBgCB(obj.nruns,c,i) =  {unique(tmp(tmp<obj.tstop))};
                    obj.SpikesBgCB{obj.nruns,c,i}(obj.SpikesBgCB{obj.nruns,c,i}<0) = [];
                end
            end
            
            for c = 1:obj.nCR
                for i = 1:num_inhc
                    tmp =  sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(round(obj.tstop/10000),1)*obj.tstop],'ascend') ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
                    obj.SpikesBgCR(obj.nruns,c,i) =  {unique(tmp(tmp<obj.tstop))};
                    obj.SpikesBgCR{obj.nruns,c,i}(obj.SpikesBgCR{obj.nruns,c,i}<0) = [];
                end
            end
            % end
        end
        function exportNetworkPositions(obj)
            % Export cells position in 3D space for later reference.
            
            fid = fopen([obj.path,obj.SLASH,'networkPositionsPyramidals.txt'],'w');
            for i=1:size(PCpos,1)
                fprintf(fid,'%f %f %f\n',obj.PCsomata(i,1),obj.PCsomata(i,2),obj.PCsomata(i,3));
            end
            fclose(fid);
            
            fid = fopen([obj.path,obj.SLASH,'networkPositionsParvalbumin.txt'],'w');
            for i=1:size(PVpos,1)
                fprintf(fid,'%f %f %f\n',obj.PVsomata(i,1),obj.PVsomata(i,2),obj.PVsomata(i,3));
            end
            fclose(fid);
            
        end
        
        function exportParams2Neuron(obj)
            % Export stimulation parameters in .hoc file:
            fid = fopen([obj.path,obj.SLASH,'importStimulationParameters.hoc'],'W');
            fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
            fprintf(fid,'// Override variables\n');
            
            fprintf(fid,'// Object decleration:\n');
            fprintf(fid,'objref PcellStimListSTR\n');
            fprintf(fid,'objref PcellStimListRND\n');
            fprintf(fid,sprintf('PcellStimListSTR = new Matrix(%d,%d)\n',obj.stimCellsPerCluster,obj.NC_str));
            fprintf(fid,sprintf('PcellStimListRND = new Matrix(%d,%d)\n',obj.stimCellsPerCluster,obj.NC_rnd));
%             fprintf(fid,'PcellStimList = new Vector(nPCcells)\n');
            
            fprintf(fid,'\n\n// Import parameters:\n\n');
            % Structured network stimulation:
            for i=1:obj.stimCellsPerCluster
                for j=1:obj.NC_str
                fprintf(fid,'PcellStimListSTR.x[%d][%d] = %d\n',i-1,j-1, obj.StimMat_str(i,j)-1); % NEURON idx
                end
            end
            % Random network stimulation:
            for i=1:obj.stimCellsPerCluster
                for j=1:obj.NC_rnd
                fprintf(fid,'PcellStimListRND.x[%d][%d] = %d\n',i-1,j-1, obj.StimMat_rnd(i,j)-1); % NEURON idx
                end
            end
            fclose(fid);
            
            % Export Network Connectivity:
            % Export STRUCTURED parameter matrices in .hoc file:
            fid = fopen([obj.path,obj.SLASH,'importNetworkParametersSTR.hoc'],'W');
            fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
            fprintf(fid,'// Override variables\n');
            
            fprintf(fid,sprintf('nPCcells=%d\n',obj.nPC));
            fprintf(fid,sprintf('nPVcells=%d\n',obj.nPV));
            fprintf(fid,sprintf('nCBcells=%d\n',obj.nCB));
            fprintf(fid,sprintf('nCRcells=%d\n',obj.nCR));
            fprintf(fid,sprintf('nAllCells=%d\n\n',obj.nAll));
            fprintf(fid,sprintf('tstop=%d\n\n',obj.tstop));
            fprintf(fid,sprintf('steps_per_ms=%d\n\n',obj.dt));
            fprintf(fid,sprintf('TOTALRUNS=%d\n\n',obj.nruns));
            
            
            fprintf(fid,'// Object decleration:\n');
            fprintf(fid,'objref connMatrix, weightsMatrix\n');
            fprintf(fid,'connMatrix = new Matrix(nAllCells, nAllCells)\n');
            fprintf(fid,'weightsMatrix = new Matrix(nAllCells, nAllCells)\n');
            
            fprintf(fid,'\n\n// Import parameters: (long-long text following!)\n\n');
            % network connectivity:
            for i=1:length(obj.state_str)
                for j=1:length(obj.state_str)
                    fprintf(fid,'connMatrix.x[%d][%d] = %d\n',i-1,j-1, obj.state_str(i,j));
                end
            end
            % Network synaptic weights
            for i=1:length(obj.weights_str)
                for j=1:length(obj.weights_str)
                    fprintf(fid,'weightsMatrix.x[%d][%d] = %f\n',i-1,j-1, obj.weights_str(i,j));
                end
            end
            fclose(fid);
            
            % Export RANDOM parameter matrices in .hoc file:
            fid = fopen([obj.path,obj.SLASH,'importNetworkParametersRND.hoc'],'W');
            fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
            fprintf(fid,'// Override variables\n');
            
            fprintf(fid,sprintf('nPCcells=%d\n',obj.nPC));
            fprintf(fid,sprintf('nPVcells=%d\n',obj.nPV));
            fprintf(fid,sprintf('nCBcells=%d\n',obj.nCB));
            fprintf(fid,sprintf('nCRcells=%d\n',obj.nCR));
            fprintf(fid,sprintf('nAllCells=%d\n\n',obj.nAll));
            fprintf(fid,sprintf('tstop=%d\n\n',obj.tstop));
            fprintf(fid,sprintf('steps_per_ms=%d\n\n',obj.dt));
            fprintf(fid,sprintf('TOTALRUNS=%d\n\n',obj.nruns));
            
            fprintf(fid,'// Object decleration:\n');
            fprintf(fid,'objref connMatrix, weightsMatrix\n');
            fprintf(fid,'connMatrix = new Matrix(nAllCells, nAllCells)\n');
            fprintf(fid,'weightsMatrix = new Matrix(nAllCells, nAllCells)\n');
            
            fprintf(fid,'\n\n// Import parameters: (long-long text following!)\n\n');
            % network connectivity:
            for i=1:length(obj.state_rnd)
                for j=1:length(obj.state_rnd)
                    fprintf(fid,'connMatrix.x[%d][%d] = %d\n',i-1,j-1, obj.state_rnd(i,j));
                end
            end
            % Network synaptic weights
            for i=1:length(obj.weights_rnd)
                for j=1:length(obj.weights_rnd)
                    fprintf(fid,'weightsMatrix.x[%d][%d] = %f\n',i-1,j-1, obj.weights_rnd(i,j));
                end
            end
            fclose(fid);

            % Export stimulation spikes in .hoc file:
            fid = fopen([obj.path,obj.SLASH,'importNetworkStimulation.hoc'],'W');
            fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
            
            fprintf(fid,'// Object decleration:\n');
            fprintf(fid,'objref Stim_Dend[%d][%d], Stim_Apic[%d][%d]\n',...
                length(obj.SpikesStimDend),length(obj.SpikesStimApic),size(obj.SpikesStimDend,2),size(obj.SpikesStimApic,2));
            
            fprintf(fid,'\n\n// Import parameters: \n\n');
            % Only for Pyramidals:
            for c=1:size(obj.SpikesStimDend,1)
                for i=1:size(obj.SpikesStimDend,2)
                    fprintf(fid,'Stim_Dend[%d][%d] = new Vector(%d)\n',c-1,i-1,length(obj.SpikesStimDend{c,i}));
                    for j=1:length(obj.SpikesStimDend{c,i})
                        fprintf(fid,'Stim_Dend[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(obj.SpikesStimDend{c,i}(j)));
                    end
                end
            end
            
            
            for c=1:size(obj.SpikesStimApic,1)
                for i=1:size(obj.SpikesStimApic,2)
                    fprintf(fid,'Stim_Apic[%d][%d] = new Vector(%d)\n',c-1,i-1,length(obj.SpikesStimApic{c,i}));
                    for j=1:length(obj.SpikesStimApic{c,i})
                        fprintf(fid,'Stim_Apic[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(obj.SpikesStimApic{c,i}(j)));
                    end
                end
            end
            
            fclose(fid);
            
            % Export Background spikes in .hoc file:
            fid = fopen([obj.path,obj.SLASH,'importBackgroundStimParams.hoc'],'W');
            
            fprintf(fid,'// Object decleration:\n');
            fprintf(fid,'objref BG_Stim_basal[%d][%d][%d], BG_Stim_Apicpr[%d][%d][%d], BG_Stim_Apic[%d][%d][%d], BG_Stim_SomaPV[%d][%d][%d], BG_Stim_SomaCB[%d][%d][%d], BG_Stim_SomaCR[%d][%d][%d]\n',...
                size(obj.SpikesBgBasal,1),size(obj.SpikesBgBasal,2),size(obj.SpikesBgBasal,3),...
                size(obj.SpikesBgProximal,1),size(obj.SpikesBgProximal,2),size(obj.SpikesBgProximal,3),...
                size(obj.SpikesBgApical,1), size(obj.SpikesBgApical,2),size(obj.SpikesBgApical,3),...
                size(obj.SpikesBgPV,1),size(obj.SpikesBgPV,2),size(obj.SpikesBgPV,3),...
                size(obj.SpikesBgCB,1),size(obj.SpikesBgCB,2),size(obj.SpikesBgCB,3),...
                size(obj.SpikesBgCR,1),size(obj.SpikesBgCR,2),size(obj.SpikesBgCR,3));
            
            fprintf(fid,'BG_dendSyn = %d\n',size(obj.SpikesBgBasal,3));
            fprintf(fid,'BG_apicSyn = %d\n',size(obj.SpikesBgApical,3));
            fprintf(fid,'BG_apicprSyn = %d\n',size(obj.SpikesBgProximal,3));
            fprintf(fid,'BG_PVSyn = %d\n',size(obj.SpikesBgPV,3));
            fprintf(fid,'BG_CBSyn = %d\n',size(obj.SpikesBgCB,3));
            fprintf(fid,'BG_CRSyn = %d\n',size(obj.SpikesBgCR,3));
            
            
            
            fprintf(fid,'\n\n// Import parameters: \n\n');
            for many_times=1:size(obj.SpikesBgBasal,1)
                for c=1:size(obj.SpikesBgBasal,2)
                    for i=1:size(obj.SpikesBgBasal,3)
                        fprintf(fid,'BG_Stim_basal[%d][%d][%d] = new Vector(%d)\n',many_times-1, c-1,i-1,length(obj.SpikesBgBasal{many_times, c,i}));
                        for j=1:length(obj.SpikesBgBasal{many_times, c,i})
                            fprintf(fid,'BG_Stim_basal[%d][%d][%d].x[%d] = %d\n',many_times-1, c-1,i-1,j-1, abs(obj.SpikesBgBasal{many_times,c,i}(j)));
                        end
                    end
                end
            end
            
            for many_times=1:size(obj.SpikesBgProximal,1)
                for c=1:size(obj.SpikesBgProximal,2)
                    for i=1:size(obj.SpikesBgProximal,3)
                        fprintf(fid,'BG_Stim_Apicpr[%d][%d][%d] = new Vector(%d)\n',many_times-1,c-1,i-1,length(obj.SpikesBgProximal{many_times, c,i}));
                        for j=1:length(obj.SpikesBgProximal{many_times, c,i})
                            fprintf(fid,'BG_Stim_Apicpr[%d][%d][%d].x[%d] = %d\n',many_times-1,c-1,i-1,j-1, abs(obj.SpikesBgProximal{many_times, c,i}(j)));
                        end
                    end
                end
            end
            
            for many_times=1:size(obj.SpikesBgApical,1)
                for c=1:size(obj.SpikesBgApical,2)
                    for i=1:size(obj.SpikesBgApical,3)
                        fprintf(fid,'BG_Stim_Apic[%d][%d][%d] = new Vector(%d)\n',many_times-1,c-1,i-1,length(obj.SpikesBgApical{many_times, c,i}));
                        for j=1:length(obj.SpikesBgApical{many_times, c,i})
                            fprintf(fid,'BG_Stim_Apic[%d][%d][%d].x[%d] = %d\n',many_times-1,c-1,i-1,j-1, abs(obj.SpikesBgApical{many_times, c,i}(j)));
                        end
                    end
                end
            end
            
            for many_times=1:size(obj.SpikesBgPV,1)
                for c=1:size(obj.SpikesBgPV,2)
                    for i=1:size(obj.SpikesBgPV,3)
                        fprintf(fid,'BG_Stim_SomaPV[%d][%d][%d] = new Vector(%d)\n',many_times-1,c-1,i-1,length(obj.SpikesBgPV{many_times, c,i}));
                        for j=1:length(obj.SpikesBgPV{many_times, c,i})
                            fprintf(fid,'BG_Stim_SomaPV[%d][%d][%d].x[%d] = %d\n',many_times-1,c-1,i-1,j-1, abs(obj.SpikesBgPV{many_times, c,i}(j)));
                        end
                    end
                end
            end
            
            for many_times=1:size(obj.SpikesBgCB,1)
                for c=1:size(obj.SpikesBgCB,2)
                    for i=1:size(obj.SpikesBgCB,3)
                        fprintf(fid,'BG_Stim_SomaCB[%d][%d][%d] = new Vector(%d)\n',many_times-1,c-1,i-1,length(obj.SpikesBgCB{many_times, c,i}));
                        for j=1:length(obj.SpikesBgCB{many_times, c,i})
                            fprintf(fid,'BG_Stim_SomaCB[%d][%d][%d].x[%d] = %d\n',many_times-1,c-1,i-1,j-1, abs(obj.SpikesBgCB{many_times, c,i}(j)));
                        end
                    end
                end
            end
            
            for many_times=1:size(obj.SpikesBgCR,1)
                for c=1:size(obj.SpikesBgCR,2)
                    for i=1:size(obj.SpikesBgCR,3)
                        fprintf(fid,'BG_Stim_SomaCR[%d][%d][%d] = new Vector(%d)\n',many_times-1,c-1,i-1,length(obj.SpikesBgCR{many_times, c,i}));
                        for j=1:length(obj.SpikesBgCR{many_times, c,i})
                            fprintf(fid,'BG_Stim_SomaCR[%d][%d][%d].x[%d] = %d\n',many_times-1,c-1,i-1,j-1, abs(obj.SpikesBgCR{many_times, c,i}(j)));
                        end
                    end
                end
            end
            
            fclose(fid);
        end
        
        function S = saveobj(obj)
%             % Put everything in a struct and save to mat:
%             S.id = obj.id;
%             S.path = obj.path;
%             S.pathToHere = obj.pathToHere;
%             S.nruns = obj.nruns;
%             S.state = obj.state;
%             S.stimstart = obj.stimstart;
%             S.stimend = obj.stimend;
%             S.stimdur = obj.stimdur;
%             S.tstop = obj.tstop;
%             S.dt = obj.dt;
%             S.nPC = obj.nPC;
%             S.nPV = obj.nPV;
%             S.nCB = obj.nCB;
%             S.nCR = obj.nCR;
%             S.nAll = obj.nAll;
%             S.PCsomata = obj.PCsomata;
%             S.PVsomata = obj.PVsomata;
%             S.CBsomata = obj.CBsomata;
%             S.CRsomata = obj.CRsomata;
%             S.connBinsPC2PC = obj.connBinsPC2PC;
%             S.connProbsPC2PC = obj.connProbsPC2PC;
%             S.connBinsPV2PC = obj.connBinsPV2PC;
%             S.connProbsPV2PC = obj.connProbsPV2PC;
%             S.connBinsCB2PC = obj.connBinsCB2PC;
%             S.connProbsCB2PC = obj.connProbsCB2PC;
%             S.recipBinsPC2PC = obj.recipBinsPC2PC;
%             S.recipProbsPC2PC = obj.recipProbsPC2PC;
%             S.NconnBins = obj.NconnBins;
%             S.NincomingProbs = obj.NincomingProbs;
%             S.NoutgoingProbs = obj.NoutgoingProbs;
%             S.distPC2PC = obj.distPC2PC;
%             S.distPV2PC = obj.distPV2PC;
%             S.distCB2PC = obj.distCB2PC;
%             S.ConnMatPC2PV = obj.ConnMatPC2PV;
%             S.ConnMatPV2PC = obj.ConnMatPV2PC;
%             S.ConnMatCB2PC = obj.ConnMatCB2PC;
%             S.connmatrix = obj.connmatrix;
%                         
%             S.state_str = obj.state_str;
%             S.state_rnd = obj.state_rnd;
%             S.weights_str = obj.weights_str;
%             S.weights_rnd = obj.weights_rnd;
%             
%             S.mergedCN_str = obj.mergedCN_str;
%             S.mergedCN_rnd = obj.mergedCN_rnd;
%             S.NC_str = obj.NC_str;
%             S.NC_rnd = obj.NC_rnd;
%             S.cellsPerCluster_str = obj.cellsPerCluster_str;
%             S.cellsPerCluster_rnd = obj.cellsPerCluster_rnd;
%             S.stimCellsPerCluster = obj.stimCellsPerCluster;
%             
%             S.labels_str = obj.labels_str;
%             S.labels_rnd = obj.labels_rnd;
%             
%             S.StimMat_str = obj.StimMat_str;
%             S.StimMat_rnd = obj.StimMat_rnd;
%             
%             S.SpikesStimDend = obj.SpikesStimDend;
%             S.SpikesStimApic = obj.SpikesStimApic;
%             
%             S.SpikesBgBasal = obj.SpikesBgBasal;
%             S.SpikesBgProximal = obj.SpikesBgProximal;
%             S.SpikesBgApical = obj.SpikesBgApical;
%             S.SpikesBgPV = obj.SpikesBgPV;
%             S.SpikesBgCB = obj.SpikesBgCB;
%             S.SpikesBgCR = obj.SpikesBgCR;
%             
%             S.SLASH = obj.SLASH;
S = obj;
        end
        %

    end
    
    methods (Static)
%         function newObj = loadobj(S)
%             % Put everything in a struct and save to mat:
%             S.id = obj.id;
%             S.path = obj.path;
%             S.pathToHere = obj.pathToHere;
%             S.nruns = obj.nruns;
%             S.state = obj.state;
%             S.stimstart = obj.stimstart;
%             S.stimend = obj.stimend;
%             S.stimdur = obj.stimdur;
%             S.tstop = obj.tstop;
%             S.dt = obj.dt;
%             S.nPC = obj.nPC;
%             S.nPV = obj.nPV;
%             S.nCB = obj.nCB;
%             S.nCR = obj.nCR;
%             S.nAll = obj.nAll;
%             S.PCsomata = obj.PCsomata;
%             S.PVsomata = obj.PVsomata;
%             S.CBsomata = obj.CBsomata;
%             S.CRsomata = obj.CRsomata;
%             S.connBinsPC2PC = obj.connBinsPC2PC;
%             S.connProbsPC2PC = obj.connProbsPC2PC;
%             S.connBinsPV2PC = obj.connBinsPV2PC;
%             S.connProbsPV2PC = obj.connProbsPV2PC;
%             S.connBinsCB2PC = obj.connBinsCB2PC;
%             S.connProbsCB2PC = obj.connProbsCB2PC;
%             S.recipBinsPC2PC = obj.recipBinsPC2PC;
%             S.recipProbsPC2PC = obj.recipProbsPC2PC;
%             S.NconnBins = obj.NconnBins;
%             S.NincomingProbs = obj.NincomingProbs;
%             S.NoutgoingProbs = obj.NoutgoingProbs;
%             S.distPC2PC = obj.distPC2PC;
%             S.distPV2PC = obj.distPV2PC;
%             S.distCB2PC = obj.distCB2PC;
%             S.ConnMatPC2PV = obj.ConnMatPC2PV;
%             S.ConnMatPV2PC = obj.ConnMatPV2PC;
%             S.ConnMatCB2PC = obj.ConnMatCB2PC;
%             S.connmatrix = obj.connmatrix;
%                         
%             S.state_str = obj.state_str;
%             S.state_rnd = obj.stateRND;
%             S.weights_str = obj.weights_str;
%             S.weights_rnd = obj.weights_rnd;
%             
%             S.mergedCN_str = obj.mergedCN_str;
%             S.mergedCN_rnd = obj.mergedCN_rnd;
%             S.NC_str = obj.NC_str;
%             S.NC_rnd = obj.NC_rnd;
%             S.cellsPerCluster_str = obj.cellsPerCluster_str;
%             S.cellsPerCluster_rnd = obj.cellsPerCluster_rnd;
%             S.stimCellsPerCluster = obj.stimCellsPerCluster;
%             
%             S.labels_str = obj.labels_str;
%             S.labels_rnd = obj.labels_rnd;
%             
%             S.StimMat_str = obj.StimMat_str;
%             S.StimMat_rnd = obj.StimMat_rnd;
%             
%             S.SpikesStimDend = obj.SpikesStimDend;
%             S.SpikesStimApic = obj.SpikesStimApic;
%             
%             S.SpikesBgBasal = obj.SpikesBgBasal;
%             S.SpikesBgProximal = obj.SpikesBgProximal;
%             S.SpikesBgApical = obj.SpikesBgApical;
%             S.SpikesBgPV = obj.SpikesBgPV;
%             S.SpikesBgCB = obj.SpikesBgCB;
%             S.SpikesBgCR = obj.SpikesBgCR;
%             
%             S.SLASH = obj.SLASH;
%         end
        function [y]=expamp(xq)
            x0 = linspace(0,1,20);
            y0 = [0,0.75,0.95,1.1,1.25,1.4,1.45,1.49,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5];
            y = interp1(x0,y0,xq);
        end
        function Points = CreateRandomNetwork(cellNo, maxSeperationDistance, RandomDimension)
            initPointsNo = cellNo;
            DistTolerance = 2;
            Points = rand(initPointsNo,RandomDimension);
            mDist = zeros((initPointsNo+1)^2,1);
            for i=1:initPointsNo
                for j=1:initPointsNo
                    %         My mean is affected from zero values (diagonal)
                    mDist((i-1)*(initPointsNo)+j) = mDist((i-1)*(initPointsNo)+j) + sqrt((Points(i,1)-Points(j,1))^2 + (Points(i,2)-Points(j,2))^2 + (Points(i,3)-Points(j,3))^2);
                end
            end
            scrDist = max(mDist) ;
            redoDist = (scrDist < maxSeperationDistance+DistTolerance) || (scrDist > maxSeperationDistance-DistTolerance);
            
            if redoDist
                DistFactor = scrDist / maxSeperationDistance;
                Points = Points / DistFactor;
            end
        end
        function somata = CreateCubeNetworkPV(cubeSize, cellNo)
            % Initializes neuronal PV somata in a cube of given dimentions
            % arg_1 the value of cube side in ?m.
            % arg_2 is the number of cells
            somata = rand(cellNo,3);
            somata = somata .* cubeSize;
        end
        function distMat = distancePV2PC(PV,PC)
            distMat = zeros(size(PV,1),size(PC,1));
            for i=1:size(PV,1)
                distMat(i,1:size(PC,1)) = reshape(sqrt( sum((PC(:,1:3)-repmat(PV(i,1:3),[size(PC,1),1])).^2,2) ),1,[]);
            end
        end
        function ConnMatCB2PC = connectCB2PC(dist,connBins,connProbs)
            % rng('shuffle');
            ConnMatCB2PC = zeros(size(dist));
            for i=1:size(dist,1)
                for j=1:size(dist,2)
                    [~,bin] = histc(dist(i,j),connBins);
                    if(rand(1) <= connProbs(bin))
                        ConnMatCB2PC(i,j) = 1;
                    end
                end
            end
        end
        function connmat = connectPC2PV(dist)
            % randomly (until newer data) connect PC 2 PV with a chance of 23%.
            % as in:
            % Kawaguchi, 2009
            connmat = zeros(size(dist));
            for i=1:size(dist,1) % PC
                for j=1:size(dist,2) % PV
                    if(rand(1) <= 0.23)
                        connmat(i,j) = 1;
                    end
                end
            end
        end

        
    end
    
end

