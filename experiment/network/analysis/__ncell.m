classdef ncell
    %NCELL, class to hold neuron info
    %   ie firing frequency (mean or windowed), position etc..
    % for bugs and comments: stamatiad.st at gmail.com
    
    properties
        mv; % in milivolts
        tstop; % in miliseconds
        dt; % steps per milisecond
        freq; % in Hz
        nspikes; % Number of spikes
        spikes; % timings of spikes in miliseconds
        position; % 3d position (3d vector)
        volume;
        meanLength;
        meanDiam;
        clusterID;
        stBins; % spike train bins
        persistentActivity; % if persistent activity
        tree % hierarchical tree (binary) structure
        soma % array (cell) of structures for each section of the soma
        axon % array (cell) of structures for each section of the axon
        dends % array (cell) of structures for each dendritic section
        % (from branch point to branch point)
        apics % array (cell) of structures for each apical section
        
    end
    
    methods %(Access = public)
        function obj = ncell(varargin)
            % Constructon of the class
            % membrane,dt,swc
            
            % initialize:
            mv=[]; % in milivolts
            tstop=[]; % in miliseconds
            dt=[]; % steps per milisecond
            freq=[]; % in Hz
            nspikes=[]; % Number of spikes
            spikes=[]; % timings of spikes in miliseconds
            position=[]; % 3d position (3d vector)
            volume=[];
            meanLength=[];
            meanDiam=[];
            clusterID=[];
            stBins=[]; % spike train bins
            persistentActivity=[]; % if persistent activity
            tree = struct; % hierarchical tree (binary) structure
            soma = {};
            axon = {};
            dends = {}; % array (cell) of structures for each dendritic section
            apics  = {};% array (cell) of structures for each apical section
            
            if nargin > 0
                membrane = varargin{1};
                if(mod(length(varargin{1}),2))
                    membrane = varargin{1}(1:end-1);
                end
                obj.mv=membrane; % in milivolts
            end
            
            if nargin > 1
                obj.tstop = length(membrane) / varargin{2};
                obj.dt = varargin{2};
                [number_of_spikes, spike_timing] = obj.spike_count(1,obj.tstop);
                obj.freq = number_of_spikes / (obj.tstop / 1000);
                obj.nspikes = number_of_spikes;
                obj.spikes = spike_timing / varargin{2};
            end
            if nargin  > 2
                % Load tree from swc file:
                obj.tree = load_tree(varargin{3});
                % disect soma:
                obj.soma = obj.detectBranches(1);
                % disect axon:
                obj.axon = obj.detectBranches(2);
                % disect basal dendrites:
                obj.dends = obj.detectBranches(3);
                % disect apical dendrites:
                obj.apics = obj.detectBranches(4);
                
            end
            
        end
        function [number_of_spikes, spike_timing] = spike_count(obj,on,off)
            % Returns spikes by detecting zero-crossings that occur when the membrane
            % potential is positive.
            if on ~= 1
                on = on * obj.dt;
            end
            if off ~= 1
                off = off * obj.dt;
            end
            responce = obj.mv(on:off);
            if size(responce,2) > size(responce,1)
                responce = responce';
            end
            spike_timing = find( ([0;diff(sign(diff(responce)))<0;0] & [sign(responce)==1]) );
            number_of_spikes = length(spike_timing);
        end
        % Function to return branch IDXs and points
        function branches = detectBranches(obj,dendType)
            %             mycell = {};
            % Verify trees structure ( must be binary tree!):
            ver_tree (obj.tree);
            % Extract BCT:
            RAW_BTC = sum(full(obj.tree.dA),1) ;
            % trace backwards terminal/branch points to extract a tree
            % branch
            branch = struct();
            branches={};
            ctr=1;
            for i= find( [ ((RAW_BTC==0) | (RAW_BTC>1)) &...
                    (cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == dendType) ] )
                %Add the terminal node to the list:
                branch(1).id= i;
                branch(1).x = obj.tree.X(i);
                branch(1).y = obj.tree.Y(i);
                branch(1).z = obj.tree.Z(i);
                branch(1).d = obj.tree.D(i);
                %         Recursively add the rest of the terminal branch nodes:
                branch = obj.recursionIn(branch,i);
                branches(ctr) = {branch(:)};
                ctr = ctr+1;
            end
        end
        function branch = recursionIn(obj,branch,idx)
            % Let recursion kick in to save lines of code: continue until hit a
            % bifurcation node.
            
            % Add the index of the parent node;
            tmpIdx = find(obj.tree.dA(idx,:) );
            tmpLoc = size(branch,2)+1;
            branch(tmpLoc).id=tmpIdx;
            branch(tmpLoc).x = obj.tree.X(tmpIdx);
            branch(tmpLoc).y = obj.tree.Y(tmpIdx);
            branch(tmpLoc).z = obj.tree.Z(tmpIdx);
            branch(tmpLoc).d = obj.tree.D(tmpIdx);
            % If the parent is continuoum point, continue the recursion:
            if sum(obj.tree.dA(:,branch(end).id))==1
                branch = obj.recursionIn(branch,branch(end).id);
            else
                % If the parent is bifurcation, STOP
                return;
            end
        end
        function subtree = getSubtree(obj,strID)
            % decode names (SWC naming convention here)
            switch strID
                case 'soma'
                    sID = 1;
                case 'axon'
                    sID = 2;
                case 'dend'
                    sID = 3;
                case 'apic'
                    sID = 4;
                otherwise
                    error('Not known subregion name! Exiting getSubtree()');
            end
            
            % Keep branches indicated by subtreeID (ex basal == 3):
            keep = find(cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == sID);
            
            % initialize the temp tree:
            subtree.dA = obj.tree.dA(keep,keep);
            subtree.X = obj.tree.X(keep);
            subtree.Y = obj.tree.Y(keep);
            subtree.Z = obj.tree.Z(keep);
            subtree.D = obj.tree.D(keep);
            subtree.R = obj.tree.R(keep);
            subtree.rnames=obj.tree.rnames;
            subtree.name = [obj.tree.name,'_',strID];
            
        end
        function obj = swapSubtree(obj,subtree,strID)
            %given an empty subtree, exit:
            if( isempty(subtree.dA) )
                warning('EMPTY subtree: no swapping done.');
                return;
            end
            
            % Temprorary: deleting the root (ie soma) is not implemented
            if(strID == 'soma')
                warning('Deleting the root (ie soma) is not implemented. No swapping done.');
                return;
            end
            
            % decode names (SWC naming convention here)
            switch strID
                case 'soma'
                    sID = 1;
                case 'axon'
                    sID = 2;
                case 'dend'
                    sID = 3;
                case 'apic'
                    sID = 4;
                otherwise
                    error('Not known subregion name! Exiting swapSubtree()');
            end
            
            % Verify subtree structure ( must be binary tree!):
            ver_tree (subtree);
            
            % remove old subtree:
            % TREESTOOLBOX HAS BUG IN delete_tree()
            obj.tree = delete_tree(obj.tree,...
                find(cellfun(@str2num,{obj.tree.rnames{obj.tree.R}}) == sID) );
            
            % Concatenate trees:
            somaXYZ = [obj.tree.X(1), obj.tree.Y(1), obj.tree.Z(1)];
            subtreeXYZ = [subtree.X(1), subtree.Y(1), subtree.Z(1)];
            
            obj.tree = cat_tree(obj.tree,tran_tree(subtree,somaXYZ-subtreeXYZ),1,1);
            
            % Update properties:
            % disect soma:
            obj.soma = obj.detectBranches(1);
            % disect axon:
            obj.axon = obj.detectBranches(2);
            % disect basal dendrites:
            obj.dends = obj.detectBranches(3);
            % disect apical dendrites:
            obj.apics = obj.detectBranches(4);
            
        end
        function obj = hasPersistent(obj,stim,freq,dur)
            %returns true if responce similars persistent activity
            %duration in ms
            %IGNORES the stimulus duration!
            for i=stim:obj.tstop - dur
                if( (spike_count(obj,i,i+dur)/ ((dur)/1000)) >= freq )
                    obj.persistentActivity = 1;
                    return;
                end
            end
            obj.persistentActivity = 0;
        end
        function obj = binSpikes(obj,binSize)
            % Function to return bins containing the spikes of the cell
            % The binSize to be given is in miliseconds
            
            % check if obj has freq (must have due to constructor):
            if(isempty(obj.freq))
                error('spike_count have not run forncell obj!');
            end
            obj.stBins = histc(obj.spikes,1:binSize:(obj.tstop*obj.dt));
        end
    end
    
end

