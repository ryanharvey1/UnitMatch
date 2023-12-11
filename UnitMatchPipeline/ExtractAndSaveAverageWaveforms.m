function Path4UnitNPY = ExtractAndSaveAverageWaveforms(clusinfo, param)

%% Called by UnitMatch, but can also be used on its own to save two averaged waveforms per unit per session

%% Read in from param
RedoExtraction = param.RedoExtraction; % Raw waveform and parameter extraction
AllDecompPaths = param.AllDecompPaths;
sampleamount = param.sampleamount; %500; % Nr. waveforms to include
spikeWidth = param.spikeWidth; %83; % in sample space (time)
halfWidth = floor(spikeWidth/2);

%% Extract all cluster info
AllClusterIDs = clusinfo.cluster_id;
if param.GoodUnitsOnly
    Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
else
    Good_Idx = 1:length(clusinfo.Good_ID);
    disp('Use all units including MUA and noise')
end
GoodRecSesID = clusinfo.RecSesID(Good_Idx);

% Define day stucture
nclus = length(Good_Idx);

%% Actual extraction
dataTypeNBytes = numel(typecast(cast(0, 'uint16'), 'uint8')); % Define datatype

% Initialize
Path4UnitNPY = cell(1, nclus);

load(fullfile(param.basepaths{1}, [basenameFromBasepath(param.basepaths{1}), '.session.mat']), 'session')
channel_idx = [session.extracellular.electrodeGroups.channels{:}];


timercounter = tic;
fprintf(1, 'Extracting raw waveforms. Progress: %3d%%', 0)

% load spikes for all sessions
nfiles = length(param.basepaths);
spAll = cell(1, nfiles);
for basepath_i = 1:nfiles
    basepath = param.basepaths{basepath_i};
    load(fullfile(basepath, [basenameFromBasepath(basepath), '.spikes.cellinfo.mat']), 'spikes')
    % spikes,st = importSpikes('basepath',basepath)
    spAll{basepath_i}.st = spikes.spindices(:, 1);
    spAll{basepath_i}.spikeTemplates = spikes.spindices(:, 2);
    spAll{basepath_i}.RecSesID = ones(length(spikes.spindices(:, 1)), 1) * basepath_i;
    spAll{basepath_i}.sample_rate = spikes.sr;
end
% Add all spikedata in one spikes struct - can be used for further analysis
spAll = [spAll{:}];
spnew = struct;
fields = fieldnames(spAll(1));
for fieldid = 1:length(fields)
    try
        spnew.(fields{fieldid}) = cat(1, spAll(:).(fields{fieldid}));
    catch ME
        if strcmp(ME.message, 'Out of memory.')
            spnew.(fields{fieldid}) = cat(1, spAll(1).(fields{fieldid}));
            for tmpid = 2:length(spAll)
                spnew.(fields{fieldid}) = cat(1, spnew.(fields{fieldid}), spAll(tmpid).(fields{fieldid}));
            end
        else
            spnew.(fields{fieldid}) = cat(1, spAll(:).(fields{fieldid}));
            eval(['spnew.', fields{fieldid}, '= cat(2,sp(:).', fields{fieldid}, ');'])
        end
    end
end
sp = spnew;
sp.sample_rate = sp.sample_rate(1);
clear spikes


Currentlyloaded = 0;
for uid = 1:nclus
    fprintf(1, '\b\b\b\b%3.0f%%', uid/nclus*100)
    % if length(param.KSDir)>1
    % tmppath = dir(fullfile(param.KSDir{GoodRecSesID(uid)},'**','RawWaveforms*'));
    % tmppath = dir(fullfile(param.KSDir{GoodRecSesID(uid)}, 'RawWaveforms'));
    if ~exist(fullfile(param.KSDir{GoodRecSesID(uid)}, 'RawWaveforms'), 'dir')
        mkdir(fullfile(param.KSDir{GoodRecSesID(uid)}, 'RawWaveforms'))
    end

    % else %Stitched KS
    %     tmppath = dir(fullfile(param.KSDir{1},'**','RawWaveforms*'));
    % end
    % if length(tmppath) > 1
    %     % Probably stitched:
    %     tmppath = tmppath(GoodRecSesID(uid));
    % end
    Path4UnitNPY{uid} = fullfile(fullfile(param.KSDir{GoodRecSesID(uid)}), 'RawWaveforms', ['Unit', num2str(AllClusterIDs(Good_Idx(uid))), '_RawSpikes.npy']); %0-indexed

    if exist(Path4UnitNPY{uid}) && ~RedoExtraction
        continue
    else
        if ~(GoodRecSesID(uid) == Currentlyloaded) % Only load new memmap if not already loaded

            % Map the data
            % clear memMapData
            spikeFile = dir(AllDecompPaths{GoodRecSesID(uid)});

            n_samples = spikeFile.bytes / (param.nChannels * dataTypeNBytes);
            ap_data = memmapfile(AllDecompPaths{GoodRecSesID(uid)}, 'Format', {'int16', [param.nChannels, n_samples], 'data'});

            % memMapData = ap_data.Data.data;
            Currentlyloaded = GoodRecSesID(uid);
        end

        % Spike samples
        idx1 = (sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(uid)) & ...
            sp.RecSesID == clusinfo.RecSesID(Good_Idx(uid))) .* round(sp.sample_rate)); % Spike times in samples;

        %Extract raw waveforms on the fly - % Unit uid
        if sampleamount < length(idx1)
            spikeIndicestmp = sort(datasample(idx1, sampleamount, 'replace', false));
        else
            spikeIndicestmp = sort(idx1);
        end
        spikeMap = nan(spikeWidth, param.nChannels, length(spikeIndicestmp));
        for iSpike = 1:length(spikeIndicestmp)
            thisSpikeIdx = int32(spikeIndicestmp(iSpike));
            if thisSpikeIdx > halfWidth && (thisSpikeIdx + halfWidth) < size(ap_data.Data.data, 2) % check that it's not out of bounds

                tmp = smoothdata(double(ap_data.Data.data(channel_idx, thisSpikeIdx-halfWidth:thisSpikeIdx+halfWidth)), ...
                    floor(0.000066*sp.sample_rate), 'gaussian', floor(0.00016*sp.sample_rate));

                tmp = (tmp - mean(tmp(:, 1:floor(0.00066*sp.sample_rate)), 2))';

                tmp(:, end+1:param.nChannels) = nan(size(tmp, 1), param.nChannels-size(tmp, 2));
                % Subtract first 10 samples to level spikes
                spikeMap(:, :, iSpike) = tmp(1:spikeWidth, :);
            end
        end
        %Actual number of wavefroms
        nwavs = sum(sum(~isnan(nanmean(spikeMap, 2)), 1) == spikeWidth); % Actual number of waves
        for cv = 1:2
            if cv == 1
                wavidx = floor(1:nwavs/2);
            else
                wavidx = floor(nwavs/2+1:nwavs);
            end
            spikeMapAvg(:, :, cv) = nanmedian(spikeMap(:, :, wavidx), 3);
        end
        spikeMap = spikeMapAvg;
        clear spikeMapAvg
        writeNPY(spikeMap, Path4UnitNPY{uid})

    end
end

fprintf('\n')
disp(['Extracting raw waveforms took ', num2str(round(toc(timercounter)./60)), ' minutes for ', num2str(nclus), ' units'])
