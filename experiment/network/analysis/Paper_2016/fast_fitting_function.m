function [f, gof] = fast_fitting_function(histo, fcell, fparams,edges,nr, varargin)
%
%   author stamatiad.st@gmail.com

p = false;
for a=1:length(varargin)
    if strcmp(varargin{a},'plot');p=true;end
end

rng default;  % For reproducibility
gof = cell(1,length(fcell));
f = cell(1,length(fcell));
tic;
for k = 1:length(fcell)
    fitfun = fcell{k};
    myfittype = fittype(fitfun);
    
    
    % Try to guess distribution options:
%     nr = 10;
    StartPoints = zeros(nr,length(fparams{k,1}));
    
    for j=1:length(fparams{k,1})
        a = fparams{k,1}(j);b = fparams{k,2}(j);
        StartPoints(:,j) = (b-a).*rand(nr,1) + a ;
    end
    
%     bestfit = Inf(nr,1);
    f_tmp = cell(1,nr);
    gof_tmp = cell(1,nr);
    gof_rmse = ones(1,nr);
    edges_just = edges(1:end-1)';
    AllFitOoptions = cell(nr,1);
    for j=1:nr
        AllFitOoptions{j} = fitoptions(myfittype);
        %     myfitoptions.Robust = 'LAR';
        %     myfitoptiond.Algorithm = 'Levenberg-Marquardt';
        %     myfitoptions.Display = 'iter';
        %     myfitoptions.TolFun = 1.0e-09;
        AllFitOoptions{j}.Lower = fparams{k,1};
        AllFitOoptions{j}.Upper = fparams{k,2};
        AllFitOoptions{j}.StartPoint = StartPoints(j,:);
%         AllFitOoptions{j} = myfitoptions;
    end
    parfor j=1:nr
        %         myfitoptions.StartPoint = StartPoints(j,:);
        [f_tmp{j}, gof_tmp{j}] = fit(edges_just,histo,myfittype,AllFitOoptions{j});
        gof_rmse(j) = gof_tmp{j}.rmse;
        %         if gof_tmp.rmse < bestfit(j)
        %             bestfit = gof_tmp.rmse
        %             f{k} = f_tmp; gof{k} = gof_tmp;
        %         end
    end
    % get best fit:
    [~,idx] = min(gof_rmse);
    f{k} = f_tmp{idx}; gof{k} = gof_tmp{idx};
    if(p)
        figure;plot(f{k},edges_just,histo,'k');
        xlabel('Values');ylabel('Relative Frequency');
        title(func2str(fitfun));
    end
end

fprintf('Parameter estimation took: %f seconds\n',toc);



return