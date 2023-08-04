function FigureFlick(UMDir,user,recompute)
    %%% Go through pairs to manually curate matches. Will save a .mat file
    %%% with the label of the pairs (0: uncurated, 1: match, -1: non-match)
    %
    % Inputs:
    %   pathFiles: path to where the .fig files are (mandatory)
    %   user: user's name (for saving results) (default = 'default')
    %   recompute: recompute the labels from scratch (default = 0)
    %
    % Key presses:
    %   Right arrow: next pair
    %   Left arrow: previous pair
    %   Up arrow: label as match
    %   Down arrow: label as non-match
    %   o: label as I don't know (=uncurated)
    %   m: go to first uncurated pair
    %   p: select pair number
    %   s: SAVE

    if ~exist('user','var')
        warning('No user specified. Will use default computername.')
        [ret,name] = system('hostname');%'default';
        name = (name(1:end-1));
    end

    if ~exist('recompute','var')
        recompute = 0;
    end
    
    % Load matchtable
    UMFile = dir(fullfile(UMDir,'UnitMatch.mat'));
    load(fullfile(UMFile.folder,UMFile.name))
    if ~any(ismember(MatchTable.Properties.VariableNames,name))
        eval(['MatchTable.' name ' = zeros(height(MatchTable),1);'])
    end



    % Get list of all figures
    d = dir(fullfile(UMDir,'MatchFigures','*.fig'));
    %     d(cellfun(@(X) datetime(X)<datetime(UMFile.date),{d.date})) = []; % Remove figures that were made prior to this matchtable
    if isempty(d)
        disp(['No figures were created after date of matchtable output: ' UMFile.date])
        return
    end
    UIDs = cell2mat(cellfun(@(y) str2num(y{2}(1:strfind(y{2},'_')-1)), cellfun(@(x) strsplit(x,'UID'), {d.name}, 'uni',0),'uni',0));
    [~,sortIdx] = sort(UIDs,'ascend');
    d = d(sortIdx);
    UIDs = UIDs(sortIdx);
   
    % Fill in GUI data
    guiData = struct;
    guiData.d = d;
    guiData.name = name;
    eval(['tmpMatch = MatchTable.' name ';']);
    try
    tmpMatchIdx = arrayfun(@(X)find(MatchTable.UID1 == X & MatchTable.UID2 == X,1,'first'),UIDs);
    catch ME
        disp(ME)
        disp('Maybe you didn''t yet make the images, or older images are present in this folder, remove these')
    end
    guiData.match = tmpMatch(tmpMatchIdx)';
    guiData.pairIDs = UIDs;
    guiData.curr.pair = guiData.pairIDs(1); % Start at first
    guiData.curr.updateFig = 1;
    guiData.curr.match = guiData.match(guiData.pairIDs == guiData.curr.pair);
    guiData.showFinishBox = 1;

    % Create figure
    blindFlickGUI = figure('color','w','name',name);

    guidata(blindFlickGUI, guiData);
    updatePlot(blindFlickGUI);
end

function updatePlot(blindFlickGUI)
    % Get guidata
    guiData = guidata(blindFlickGUI);

    % Slow -- to change
    if guiData.curr.updateFig == 1
        % Clear previous fig
        clf(blindFlickGUI)

        % Load and copy new one
        currPairPosition = guiData.pairIDs == guiData.curr.pair;
        fig = openfig(fullfile(guiData.d(currPairPosition).folder,guiData.d(currPairPosition).name),'invisible');
        axNew = findobj(fig,'type','axes');
        copyobj(axNew,blindFlickGUI)
        close(fig)
    end

    % Update title with label -- hacky, should be predefined (but clear
    % fig prevents from doing that...)
    if ~isfield(guiData,'titleMain')
        guiData.titleMain = annotation(blindFlickGUI,'textbox', [0, 1, 1, 0], 'string', 'My Text', 'EdgeColor', 'none',...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
    set(guiData.titleMain, 'String', sprintf('pair UID %d: match? %d (%d/%d left)', ...
        guiData.curr.pair,guiData.curr.match,sum(guiData.match == 0),numel(guiData.pairIDs)));

    % Set functions for key presses
    set(blindFlickGUI,'WindowKeyPressFcn',@keyPress,'DeleteFcn',@guiClosure);

end

function keyPress(blindFlickGUI,eventdata)
    % Get guidata
    guiData = guidata(blindFlickGUI);
    
    guiData.curr.updateFig = 0;
    switch eventdata.Key
        case 'rightarrow' % Next pair
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            newPair = guiData.pairIDs(circshift(currPairPosition,1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'leftarrow' % Previous pair
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            newPair = guiData.pairIDs(circshift(currPairPosition,-1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'uparrow' % It's a match
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            guiData.match(currPairPosition) = 1;
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            newPair = guiData.pairIDs(circshift(currPairPosition,1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'downarrow' % It's not a match
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            guiData.match(currPairPosition) = -1;
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            newPair = guiData.pairIDs(circshift(currPairPosition,1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'o' % Back to uncurated
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            guiData.match(currPairPosition) = 0;
        case 'm' % Go to first uncurated
            firstUncuratedPair = find(guiData.match == 0,1);
            if ~isempty(firstUncuratedPair)
                guiData.curr.pair = firstUncuratedPair;
                guiData.curr.updateFig = 1;
            else
                guiData.showFinishBox = 1;
            end
        case 'p' % Go to specific UID
            newPair = str2double(cell2mat(inputdlg('Go to pair:')));
            if ~ismember(newPair,unique(guiData.pairIDs))
                error(['Pair ' num2str(newPair) ' not present'])
            end
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 's' % save
            savedata(guiData)
              
    end
    currPairPosition = guiData.pairIDs == guiData.curr.pair;
    guiData.curr.match = guiData.match(currPairPosition);
  

    if ~any(guiData.match == 0) && guiData.showFinishBox
        msgbox('All pairs have been curated. Don''t forget to save with ''s''')
        guiData.showFinishBox = 0;
    end

    guidata(blindFlickGUI,guiData);
    updatePlot(blindFlickGUI);
end

function savedata(guiData)

SaveDir = strsplit(guiData.d(1).folder,'MatchFigures');
SaveDir = SaveDir{1};
% Load MatchTable
load(fullfile(SaveDir,'UnitMatch.mat'))
if ~any(ismember(MatchTable.Properties.VariableNames,guiData.name))
    eval(['MatchTable.' guiData.name  ' = zeros(height(MatchTable),1);'])
end

match = guiData.match;
tmpMatchIdx = arrayfun(@(X)find(MatchTable.UID1 == X & MatchTable.UID2 == X,1,'first'),guiData.pairIDs);
eval(['MatchTable.' guiData.name '(tmpMatchIdx) = match;']);
disp(['Saving curation to ' fullfile(SaveDir,'UnitMatch.mat')])
save(fullfile(SaveDir,'UnitMatch.mat'),'MatchTable','-append')


end

function guiClosure(blindFlickGUI,eventdata)
answer = questdlg('Saving curation in matchtable?','Save?','Yes','No','Yes');
if strcmpi(answer,'yes')
    % Get guidata
    guiData = guidata(blindFlickGUI);

    savedata(guiData)
end


end