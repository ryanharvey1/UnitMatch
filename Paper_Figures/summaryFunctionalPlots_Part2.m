function summaryFunctionalPlots_Part2(UMFiles, whichMetric, groupVector, UseKSLabels)
    %% Will plot summary plots: distribution, ROC and AUC. 
    % UMFiles: list cells contains path to UnitMatch.m files
    % whichMetric: will compute distributions/ROC/AUC on either 'Corr', 'Rank', or 'Sig'. 
    % groupVector: same size as UMFiles, to group e.g. recordings from the
    % same mouse together.
    % UseKSLabels: use KS labels for matching

    %% Define parameters
    
    % Initialize
    if ~exist('whichMetric','var') || isempty(whichMetric)
        whichMetric = 'Rank';
    end

    if ~exist('groupVector','var')
        groupVector = 1:length(UMFiles);
    end
    groups = unique(groupVector);
    groupColor = gray(length(groups)+1);

    if ~exist('UseKSLabels','var')
        UseKSLabels = 0;
    end

    switch whichMetric
        case 'Corr'
            fprintf("Taking the correlation values!\n")
            FPNames = {'FRDiff','ACGCorr','refPopCorr','natImRespCorr'};
            stepsz = [0.1 0.1 0.1 0.1];
            minVal = [0 -1 -1 -1];
            maxVal = [15 1 1 1];
            flipROC = [0 1 1 1];
        case 'Rank'
            fprintf("Taking the rank!\n")
            FPNames = {'FRRank','ACGRank','refPopRank','natImRespRank'};
            stepsz = [1 1 1 1];
            minVal = [0.5 0.5 0.5 0.5];
            maxVal = [20.5 20.5 20.5 20.5];
            flipROC = [0 0 0 0];
        case 'Sig'
            fprintf("Taking the sig!\n")
            FPNames = {'FRSig','ACGSig','refPopSig','natImRespSig'};
            stepsz = [0.1 0.1 0.1 0.1];
            minVal = [0 0 0 0];
            maxVal = [1 1 1 1];
            flipROC = [0 0 0 0];
    end

    histBins = cell(1,numel(FPNames));
    histBinsCenter = cell(1,numel(FPNames));
    for fpIdx = 1:numel(FPNames)
        histBins{fpIdx} = minVal(fpIdx):stepsz(fpIdx):maxVal(fpIdx);
        histBinsCenter{fpIdx} = histBins{fpIdx}(1:end-1) + diff(histBins{fpIdx})/2;
        if strcmp(whichMetric,'Rank')
            histBins{fpIdx}(end+1) = inf;
            histBinsCenter{fpIdx}(end+1) = histBinsCenter{fpIdx}(end)+1;
        end
    end
    ROCBins = 0:0.01:1;
    minMatches = 20;
    durLim = 10*60;

    
    %% Loop over mice to get all Distributions / ROCs / AUCs
    
    FPSum = struct();
    deltaDays = cell(1, length(UMFiles));
    numMatchedUnits = cell(1, length(UMFiles));
    InitialDrift = cell(1,length(UMFiles));
    FixedDrift =  cell(1,length(UMFiles));
    maxAvailableUnits = cell(1, length(UMFiles));
    for midx = 1:length(UMFiles)
        %% Load data
    
        fprintf('Reference %s...\n', UMFiles{midx})
    
        tmpfile = dir(UMFiles{midx});
        if isempty(tmpfile)
            continue
        end
    
        fprintf('Loading the data...\n')
        tic
        load(fullfile(tmpfile.folder, tmpfile.name), 'MatchTable', 'UMparam');
        toc
    
        sessIDs = unique(MatchTable.RecSes1);
    
        if ~isfield(UMparam,'RawDataPaths')
            UMparam.RawDataPaths = UMparam.AllRawPaths;
        end
        %%% HACK -- Can remove later
        if ~iscell(UMparam.RawDataPaths)
            for ii = 1:numel(UMparam.RawDataPaths)
                tmp{ii} = UMparam.RawDataPaths(ii);
            end
            UMparam.RawDataPaths = tmp;
        end
    
        %% Loop through pairs of sessions
    
        fprintf('Looping through days...\n')
        tic
        days{midx} = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
        days{midx} = cell2mat(days{midx}) - days{midx}{1};
        deltaDays{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
        numMatchedUnits{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
        InitialDrift{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
        FixedDrift{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
        for sess1Idx = 1:numel(sessIDs)-1
    
            sess1 = sessIDs(sess1Idx);
            day1 = days{midx}(sess1Idx);
            meta = ReadMeta2(UMparam.RawDataPaths{sess1}.folder);
            durSess1 = str2double(meta.fileTimeSecs);
            if durSess1 < durLim 
                continue
            end

            for sess2Idx = sess1Idx+1:numel(sessIDs)
    
                sess2 = sessIDs(sess2Idx);
                day2 = days{midx}(sess2Idx);
                deltaDays{midx}(sess1Idx,sess2Idx) = abs(day2 - day1);
                meta = ReadMeta2(UMparam.RawDataPaths{sess2}.folder);
                durSess2 = str2double(meta.fileTimeSecs);
                if durSess2 < durLim
                    continue
                end
    
                %% Cut table to specific days
    
                MatchTable_2sess = MatchTable(ismember(MatchTable.RecSes1, [sess1 sess2]) & ismember(MatchTable.RecSes2, [sess1 sess2]), :);
    
                %% Number of matches
    
                if ~UseKSLabels
                    %%% CHECK THAT THIS MAKES SENSE
                    %%% CHOOSE BASED ON UID
                    matchedUnitsIdx = (MatchTable_2sess.UID1 == MatchTable_2sess.UID2) & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2); % using Unique ID
                    %%% OR RECOMPUTE
%                     [~,~,idx,~] = getPairsAcross2Sess(MatchTable_2sess, UMparam.ProbabilityThreshold);
%                     matchedUnitsIdx = zeros(size(MatchTable_2sess,1),1);
%                     matchedUnitsIdx(idx) = 1;
                else
                    matchedUnitsIdx = (MatchTable_2sess.ID1 == MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2);
                end
                numMatchedUnits{midx}(sess1Idx,sess2Idx) = sum(matchedUnitsIdx)/2; % Divided by two because looking both ways -- can be non-integer
                
                maxAvailableUnits{midx}(sess1Idx,sess2Idx) = min([length(unique(MatchTable_2sess.ID1(MatchTable_2sess.RecSes1 == sess1))) length(unique(MatchTable_2sess.ID1(MatchTable_2sess.RecSes1 == sess2)))]);%
                
                %% Extract drift if present
    
                if isfield(UMparam,'drift')
                    InitialDrift{midx}(sess1Idx,sess2Idx) =  vecnorm(UMparam.drift(sess2Idx-1,:,1),2); % Drift in recording 1 is 1 vs 2, etc.
                    FixedDrift{midx}(sess1Idx,sess2Idx) =  vecnorm(UMparam.drift(sess2Idx-1,:,2),2); % Drift in recording 1 is 1 vs 2, etc.
                end
            end
        end

        %% qParams
        if ~exist(fullfile(tmpfile.folder, 'qMetricAUCs.mat'))
            try
                QualityMetricsROCs(UMparam.SaveDir)
                close all
            catch ME
                keyboard
            end
        end
        load(fullfile(tmpfile.folder, 'qMetricAUCs.mat'))

        if ~exist('qParamNames')
            qParamNames = AUCqParams.qParamNames;
            AUCPerqParam = nan(length(qParamNames),0);
        end
        [takethese,puthere] = ismember(AUCqParams.qParamNames,qParamNames);
        tmpQP = nan(length(qParamNames),1);
        tmpQP(puthere(puthere~=0)) = AUCqParams.AUCMvNM(takethese);
        AUCPerqParam = cat(2,AUCPerqParam,tmpQP);
        toc
    end

    %% AUC distributions
    figure('name','AUC Distr')
    stepsz = 0.05;
    binvec = [0:stepsz:1];
plotvec = stepsz/2:stepsz:1-stepsz/2;
    for qid = 1:length(qParamNames)
        subplot(ceil(sqrt(length(qParamNames))),round(sqrt(length(qParamNames))),qid)
        if nanmedian(AUCPerqParam(qid,:))<0.5
            AUCPerqParam(qid,:) = 1-AUCPerqParam(qid,:);
        end
        nums = histcounts(AUCPerqParam(qid,:),binvec);
        plot(plotvec,nums,'k-');
        hold on
        line([nanmedian(AUCPerqParam(qid,:)) nanmedian(AUCPerqParam(qid,:))],get(gca,'ylim'),'color',[1 0 0])
        [h,p(qid)] = ttest(AUCPerqParam(qid,:),0.5);
        
        title([qParamNames{qid} ', p=' num2str(round(p(qid)*100)./100)])
        xlim([0 1])
        makepretty
        offsetAxes
    end




    %% Figure -- example mouse

    % Number of matched units (matrix)
    for midx = 1:numel(UMFiles)
        figure('Position', [80 700 1700 170],'Name', fileparts(fileparts(UMFiles{midx})));
        s = subplot(1,numel(FPNames)+2,1);
        imagesc(numMatchedUnits{midx})
        xticks(1:size(deltaDays{midx},2))
        xticklabels(days{midx})
        yticks(1:size(deltaDays{midx},1))
        yticklabels(days{midx}(1:end-1))
        axis equal tight
        colorbar
        c = colormap("gray"); colormap(s(1), flipud(c));
        subplot(1,numel(FPNames)+2,2); hold all
        scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}),10,groupColor(groupVector(midx),:),'filled')
        hline(minMatches)
        ylabel('Number of matches')
        xlabel('\Deltadays')
    
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            s = subplot(1,numel(FPNames)+2,fpIdx+2);
            imagesc(squeeze(FPSum.(FPNameCurr).AUC{midx}(1,:,:)))
            xticks(1:size(deltaDays{midx},2))
            xticklabels(days{midx})
            yticks(1:size(deltaDays{midx},1))
            yticklabels(days{midx}(1:end-1))
            axis equal tight
            colormap(s, "RedBlue")
            clim([0 1])
            title(sprintf('Fingerprint %s', FPNameCurr))
        end
    end
    
    %% Figure -- all mice
    % Build matrices across mice -- average in groups according to groupVector
    distMatrix = cell(1,numel(FPNames));
    ROCMatrix = cell(1,numel(FPNames));
    AUCMatrix = cell(1,numel(FPNames));
    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
        distMatrix{fpIdx} = nan(numel(histBinsCenter{fpIdx}), 3, length(groups));
        ROCMatrix{fpIdx} = nan(numel(ROCBins), 2, length(groups));
        AUCMatrix{fpIdx} = nan(2, length(groups));
        for gg = 1:length(groups)
            % Distributions
            tmp = cellfun(@(x) x(:,:,:), FPSum.(FPNameCurr).Distr(groupVector == gg), 'uni', 0);
            distMatrix{fpIdx}(:,:,gg) = nanmean(cat(3,tmp{:}),3);
    
            % ROC
            tmp = cellfun(@(x) x(:,:,:), FPSum.(FPNameCurr).ROC(groupVector == gg), 'uni', 0);
            ROCMatrix{fpIdx}(:,:,gg) = nanmean(cat(3,tmp{:}),3);
    
            % AUC
            tmp = cellfun(@(x) x(:,:), FPSum.(FPNameCurr).AUC(groupVector == gg), 'uni', 0);
            AUCMatrix{fpIdx}(:,gg) = nanmean(cat(2,tmp{:}),2);
        end
    end
    
    % Plot
    distrCols = [0 0.7 0; 1 0 0; 0 0 0.7]; % Within / Match / Non-match
    ROCCols = [1 0 0; 0 0.7 0]; % across Match vs. non-match / within Match vs. non-match
    figure('Position', [400 270 800 700]);
    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
    
        % Plot distribution
        subplot(4,numel(FPNames),0*numel(FPNames)+fpIdx); hold all
        for hid = [3 1 2]
            distr2plt = cumsum(distMatrix{fpIdx}(:,hid,:));
            h = shadedErrorBar(histBinsCenter{fpIdx}, nanmean(distr2plt,3), ...
                nanstd(distr2plt,[],3)./sqrt(sum(~isnan(distr2plt),3)));
            h.mainLine.Color = distrCols(hid,:);
            if ~isempty(h.patch)
                h.patch.FaceColor = distrCols(hid,:);
                h.edge(1).Color = 'none';
                h.edge(2).Color = 'none';
            end
        end
        title(sprintf('%s', FPNameCurr))
        xlabel(whichMetric)
        if strcmp(whichMetric,'Rank'); xticks([1 10 histBinsCenter{midx}(end)]); xticklabels({'1','10',sprintf('>%d',histBinsCenter{midx}(end-1))}); end
        ylabel('Proportion')
        offsetAxes
        makepretty
    
        % Plot ROC
        subplot(4,numel(FPNames),1*numel(FPNames)+fpIdx); hold all
        for hid = 2:-1:1
            h = shadedErrorBar(ROCBins, nanmean(ROCMatrix{fpIdx}(:,hid,:),3), ...
                nanstd(ROCMatrix{fpIdx}(:,hid,:),[],3)./sqrt(sum(~isnan(ROCMatrix{fpIdx}(:,hid,:)),3)));
            h.mainLine.Color = ROCCols(hid,:);
            if ~isempty(h.patch)
                h.patch.FaceColor = ROCCols(hid,:);
                h.edge(1).Color = 'none';
                h.edge(2).Color = 'none';
            end
        end
        plot([0 1], [0 1], 'k--')
        xlim([0 1])
        ylim([0 1])
        xticks([0 1])
        yticks([0 1])
        axis equal tight
        xlabel('False positives')
        ylabel('Hits')
        offsetAxes
        makepretty
    
        % Plot stability of AUC with delta days
        subplot(4,numel(FPNames),2*numel(FPNames)+fpIdx); hold all
        for midx = 1:length(UMFiles)
            xDays = deltaDays{midx}(:);
            xDays(xDays == 0) = 10^(-0.1);
            yVal = mat2vec(FPSum.(FPNameCurr).AUC{midx}(1,:,:));
            nanIdx = isnan(yVal);
            xDays(nanIdx) = [];
            yVal(nanIdx) = [];
            if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                scatter(log10(xDays),yVal,10,ROCCols(1,:),'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
                X = [ones(numel(xDays),1), xDays];
                b = (X\yVal);
                plot(log10(1:max(xDays)), b(1) + b(2)*(1:max(xDays)), 'color',ROCCols(1,:),'LineWidth',1);
            end
            scatter(-0.1,nanmean(mat2vec(FPSum.(FPNameCurr).AUC{midx}(2,:,:))),20,ROCCols(2,:),'filled')
        end
        xlabel('\Delta days')
        ylabel('AUC')
        xticks([-0.1 log10([1 10 100])])
        xticklabels({'within','1','10','100'})
        ylim([0 1])
        hline(0.5)
        offsetAxes
        makepretty

        % Plot stability of AUC with delta days
        subplot(4,numel(FPNames),3*numel(FPNames)+fpIdx); hold all
        for midx = 1:length(UMFiles)
            xDays = deltaDays{midx}(:);
            xDays(xDays == 0) = 10^(-0.1);
            yVal = mat2vec(FPSum.(FPNameCurr).AUC{midx}(1,:,:));
            nanIdx = isnan(yVal);
            xDays(nanIdx) = [];
            yVal(nanIdx) = [];
            if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                X = [ones(numel(xDays),1), xDays];
                b = (X\yVal);
                plot(log10(1:max(xDays)), b(1) + b(2)*(1:max(xDays)), 'color',ROCCols(1,:),'LineWidth',1);
            end
            scatter(-0.1,nanmean(mat2vec(FPSum.(FPNameCurr).AUC{midx}(2,:,:))),20,ROCCols(2,:),'filled')
        end
        xlabel('\Delta days')
        ylabel('AUC')
        xticks([-0.1 log10([1 10 100])])
        xticklabels({'within','1','10','100'})
        ylim([0 1])
        hline(0.5)
        offsetAxes
        makepretty
    end
    
    %% Additional figures
    
    % Plot number of matches as a function of delta days
    figure;
    hold all
    for midx = 1:length(UMFiles)
        scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}),10,groupColor(groupVector(midx),:),'filled')
    end
    ylabel('Number of matches')
    xlabel('\Deltadays')
    
    % AUC as a function of delta days
    figure;
    for aucIdx = 1:2
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            subplot(3,numel(FPNames),(aucIdx-1)*numel(FPNames)+fpIdx);
            hold all
            for midx = 1:length(UMFiles)
                if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                    scatter(deltaDays{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(aucIdx,:,:)),10,groupColor(groupVector(midx),:),'filled')
                end
            end
            xlabel('\Delta days')
            ylabel('AUC')
            title(sprintf('Fingerprint %s', FPNameCurr))
            ylim([0 1])
            hline(0.5)
        end
    end
    
    % Dependence of AUC on number of matched units
    figure;
    for aucIdx = 1:2
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            subplot(2,numel(FPNames),(aucIdx-1)*numel(FPNames)+fpIdx);
            hold all
            for midx = 1:length(UMFiles)
                if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                    scatter(numMatchedUnits{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(aucIdx,:,:)),10,groupColor(groupVector(midx),:),'filled')
                end
            end
            xlabel('Number of matched neurons')
            ylabel('AUC')
            title(sprintf('Fingerprint %s', FPNameCurr))
            ylim([0 1])
            hline(0.5)
        end
    end
    
    
    %% Dependence of number matched units on drift
    
%     figure('name','NrUnits versus drift')
%     for midx = 1:length(UMFiles)
%         subplot(1,2,1)
%         scatter(InitialDrift{midx},numMatchedUnits{midx}./MaxAvailableUnits{midx},20,groupColor(groupVector(midx),:),'filled')
%         hold on
%         xlabel('Drift (Eucl Distance)')
%         ylabel('Number of matches')
%         xlim([0 UMparam.NeighbourDist])
%         title('Initial Drift')
%     
%         subplot(1,2,2)
%         scatter(FixedDrift{midx},numMatchedUnits{midx}./ MaxAvailableUnits{midx},20,groupColor(groupVector(midx),:),'filled')
%         hold on
%         xlabel('Drift (Eucl Distance)')
%         ylabel('Proportion of matches')
%         xlim([0 UMparam.NeighbourDist])
%     
%     end
%     InitialDrift(cellfun(@isempty,InitialDrift)) = [];
%     InitialDrift = cellfun(@(X) X(:),InitialDrift,'uni',0);
%     InitialDrift(InitialDrift>UMparam.maxdist) = nan;
%     MaxAvailableUnits(cellfun(@isempty,MaxAvailableUnits)) = [];
%     MaxAvailableUnits = cellfun(@(X) X(:),MaxAvailableUnits,'uni',0);
%     numMatchedUnits(cellfun(@isempty,numMatchedUnits)) = [];
%     numMatchedUnits = cellfun(@(X) X(:),numMatchedUnits,'uni',0);
%     InitialDrift = cat(1,InitialDrift{:});
%     MaxAvailableUnits = cat(1,MaxAvailableUnits{:});
%     numMatchedUnits = cat(1,numMatchedUnits{:});
%     FixedDrift(cellfun(@isempty,FixedDrift)) = [];
%     FixedDrift = cellfun(@(X) X(:),FixedDrift,'uni',0);
%     FixedDrift = cat(1,FixedDrift{:});
%     
%     
%     FixedDrift(FixedDrift>UMparam.maxdist) = nan;
%     InitialDrift(InitialDrift>UMparam.maxdist) = nan;
%     
%     
%     PercNeurons = numMatchedUnits./MaxAvailableUnits;
%     [r,p] = corr(InitialDrift(~isnan(InitialDrift)),PercNeurons(~isnan(InitialDrift)));
%     subplot(1,2,1)
%     title(['Initial Drift, r=' num2str(round(r*100)/100) ', p=' num2str(round(p*100)/100)])
%     
%     subplot(1,2,2)
%     [r,p] = corr(FixedDrift(~isnan(FixedDrift)),PercNeurons(~isnan(FixedDrift)));
%     title(['Corrected Drift, r=' num2str(round(r*100)/100) ', p=' num2str(round(p*100)/100)])
end