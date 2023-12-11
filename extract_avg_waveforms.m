function save_path = extract_avg_waveforms(param)

%% Read in from param
sampleamount = param.sampleamount; %500; % Nr. waveforms to include
spikeWidth = param.spikeWidth; %83; % in sample space (time)
halfWidth = floor(spikeWidth/2);

for session_i = 1:length(param.basepaths)
    basepath = param.basepaths{session_i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath, [basename, '.session.mat']), 'session')

    d = dir(fullfile(basepath,[basename,'.dat']));
    nSamp = d.bytes/2/session.extracellular.nChannels;
    ap_data = memmapfile(fullfile(basepath,[basename,'.dat']), 'Format', {'int16', [session.extracellular.nChannels, nSamp], 'data'});

    channel_idx = [session.extracellular.electrodeGroups.channels{:}];
    load(fullfile(basepath, [basename, '.spikes.cellinfo.mat']), 'spikes')

    if ~exist(fullfile(basepath, 'RawWaveforms'),'dir')
        mkdir(fullfile(basepath, 'RawWaveforms'))
    end
    % for iCluster = 1:length(spikes.rawWaveform_all)
    %     tmpspkmap(:,:,1) = spikes.rawWaveform_all{iCluster}(channel_idx,:)';
    %     tmpspkmap(:,:,2) = spikes.rawWaveform_all{iCluster}(channel_idx,:)';
    %     writeNPY(tmpspkmap, fullfile(savePath,'RawWaveforms',['Unit' num2str(spikes.UID(iCluster)) '_RawSpikes.npy']))
    % end
    WaitMessage = parfor_wait(length(spikes.times));
    parfor unit_i = 1:length(spikes.times)
        save_path{unit_i,session_i} = extract_waveform(unit_i,spikes,ap_data,channel_idx,sampleamount,spikeWidth,halfWidth,param,basepath);
        WaitMessage.Send;
    end
    WaitMessage.Destroy;
    clear ap_data
end

save_path = vertcat(save_path(:,1),save_path(:,2));

save_path = save_path(~cellfun('isempty',save_path));

end

function save_path = extract_waveform(unit_i,spikes,ap_data,channel_idx,sampleamount,spikeWidth,halfWidth,param,basepath)

save_path = fullfile(basepath,'RawWaveforms',['Unit' num2str(unit_i) '_RawSpikes.npy']);
if exist(save_path,'file')
    return
end

st = spikes.times{unit_i} .* round(spikes.sr);
if sampleamount < length(st)
    st_ind = sort(datasample(st, sampleamount, 'replace', false));
else
    st_ind = sort(st);
end

for iSpike = 1:length(st_ind)

    thisSpikeIdx = int32(st_ind(iSpike));

    if thisSpikeIdx > halfWidth && (thisSpikeIdx + halfWidth) < size(ap_data.Data.data, 2) % check that it's not out of bounds

        tmp = smoothdata(...
            double(ap_data.Data.data(channel_idx, thisSpikeIdx-halfWidth:thisSpikeIdx+halfWidth)), ...
            floor(0.000066*spikes.sr),...
            'gaussian', floor(0.00016*spikes.sr)...
            );

        tmp = (tmp - mean(tmp(:, 1:floor(0.00066*spikes.sr)), 2))';

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

writeNPY(spikeMapAvg, save_path)
end
