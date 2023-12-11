%% DEMO UNIT MATCH 


%% READ ME
% If you do not use the suggested pipeline to extract raw waveforms (e.g. you don't use Neuropixels/SpikeGLX), make
% sure your 'KiloSortPaths' contains a subfolder called 'RawWaveforms'. There
% should be a NPY file for every cluster with the dimensions of
% UMparam.spikeWidth X nRecordingChannels X 2 (1 for each half of a
% recording). This should contain the average waveform (recommended of at
% least 500 spikes) for every recording channel for every half of a
% recording for that cluster.
UMparam.basepaths = {'Z:\Data\HMC1\day8', 'Z:\Data\HMC1\day9'};
%% User input: 
UMparam.SaveDir = 'Z:\Data\HMC1\test_unit_match'; % Recommended to use end this path with \Probe0\IMRO_1\ if more probes/IMRO tables were used or \AllProbes\AllIMRO\ otherwise
% UMparam.KSDir = {'Z:\Data\HMC3\day4\Kilosort_2023-07-28_220442', 'Z:\Data\HMC3\day5\Kilosort_2023-07-29_100711'};  % This is a cell array with a path, in the path there should be a subfolder called 'RawWaveforms'. 
UMparam.KSDir = UMparam.basepaths; % ksdir to basepath, because we are using cell explorer to pull data
% N.B. if you want to use the functional score evaluation of UnitMatch, 'KSDir' should also contain typical 'Kilosort output', (e.g. spike times etc.)

%% N.B. the following user input can also be automatically extracted and prepared/cleaned up using UMparam = ExtractKilosortData(KiloSortPaths, UMparam) for Kilosorted data of SpikeGLX recorded data (see next section);
for basepath_i = 1:length(UMparam.basepaths)
    UMparam.RawDataPaths{basepath_i} = fullfile(UMparam.basepaths{basepath_i},[basenameFromBasepath(UMparam.basepaths{basepath_i}),'.dat']);
end
UMparam.AllDecompPaths = UMparam.RawDataPaths;
% UMparam.RawDataPaths = {'Z:\Data\HMC3\day4\day4.dat', 'Z:\Data\HMC3\day5\day5.dat'};  % This is a cell array with info on where to find the decompressed recording (.cbin files) --> Necessary when you want UnitMatch to do waveform extraction
% UMparam.AllDecompPaths = {'Z:\Data\HMC3\day4\day4.dat', 'Z:\Data\HMC3\day5\day5.dat'};  % This is a cell array with info on where to find the decompressed recording (.bin files) --> Necessary when you want UnitMatch to do waveform extraction

for basepath_i = 1:length(UMparam.basepaths)
    load(fullfile(UMparam.basepaths{basepath_i},[basenameFromBasepath(UMparam.basepaths{basepath_i}),'.session.mat']),'session')
    channel_idx = [session.extracellular.electrodeGroups.channels{:}];  
    UMparam.AllChannelPos{basepath_i} = [session.extracellular.chanCoords.x(channel_idx),session.extracellular.chanCoords.y(channel_idx)];
end
% UMparam.AllChannelPos = {[RecordingSites_Recording1],[RecordingSites_Recording2]}; % These are coordinates of every recording channel on the probe (e.g. nRecordingChannels x 2)
clear clusinfo
clusinfo = struct;
clusinfo.cluster_id = [];
clusinfo.Good_ID = [];
clusinfo.RecSesID = [];
clusinfo.Probe = [];
clusinfo.ProbeID = [];
clusinfo.Depth = [];
clusinfo.Shank = [];

for basepath_i = 1:length(UMparam.basepaths)
    basepath = UMparam.basepaths{basepath_i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'spikes')
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
    % Note, this can be kilosort input, 
    % - clusinfo (this is a struct that contains per unit the following information):
    clusinfo.cluster_id = [clusinfo.cluster_id,spikes.UID];% * cluster_id (e.g. kilosort output clus_id)
    clusinfo.Good_ID = [clusinfo.Good_ID,ones(1,length(spikes.UID))];% * Good_ID: ones for units that should be included in the analysis
    clusinfo.RecSesID = [clusinfo.RecSesID,ones(1,length(spikes.UID))*basepath_i];% * RecSesID: Recording Session ID
    clusinfo.Probe = [clusinfo.Probe,ones(1,length(spikes.UID))];% * Probe: Which probe (if just 1, ones of numel cluster_id)
    clusinfo.ProbeID = [clusinfo.ProbeID,ones(1,length(spikes.UID))];
    clusinfo.Depth = [clusinfo.Depth,cell_metrics.trilat_y];% * Depth: depth on probe (optional)
    clusinfo.Shank = [clusinfo.Shank,spikes.shankID];% * Shank: Which shank (optional)
        % * Coordinates: Typically 3D Allen Common Coordinate framework coordinates per unit (optional)
end
clusinfo.cluster_id = clusinfo.cluster_id';
clusinfo.Good_ID = clusinfo.Good_ID';
clusinfo.RecSesID = clusinfo.RecSesID';
clusinfo.Probe = clusinfo.Probe';
clusinfo.ProbeID = clusinfo.ProbeID';
clusinfo.Depth = clusinfo.Depth';
clusinfo.Shank = clusinfo.Shank';


% N.B. clusinfo can also be automatically extracted using clusinfo =
% getClusinfo

% basepath = 'Z:\Data\HMC3\day4'
% make avg waveforms


%% Add paths and subpaths
% mfilePath = mfilename('fullpath');
% if contains(mfilePath,'LiveEditorEvaluationHelper')
%     mfilePath = matlab.desktop.editor.getActiveFilename;
% end
% Components = strsplit(mfilePath,filesep);
% addpath(genpath(fullfile(Components{1:end-1})));

%% Optional (for Kilosort + SpikeGLX users) --- see ExampleAnalysisPipelines for more detail!!
% UMparam = ExtractKilosortData(UMparam.KSDir, UMparam); % Extract KS data and do some noise removal, optionally decompresses cbin to bin data and uses BOMBCELL quality metric to define good single units
% clusinfo = getClusinfo(UMparam.KSDir); % prepare clusinfo struct

%% Load default parameters
UMparam = DefaultParametersUnitMatch(UMparam,session);


%% UnitMatch algorithm:
[UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);
if UMparam.AssignUniqueID
    AssignUniqueID(UMparam.SaveDir);
end

%%% N.B. From here it is all evaluation, you don't need this to use UnitMatch
%%% results in your further analysis
%% Automatic evaluation:
EvaluatingUnitMatch(UMparam.SaveDir); % Within session cross-validation
% QualityMetricsROCs(UMparam.SaveDir); % Only works in combination with BOMBCELL
ComputeFunctionalScores(UMparam.SaveDir) % Only works when having access to Kilosort output (e.g. spike times etc.) 

%% Curation:
if UMparam.MakePlotsOfPairs
    DrawPairsUnitMatch(UMparam.SaveDir);
    if UMparam.GUI
        FigureFlick(UMparam.SaveDir)
        pause
    end
end
load(fullfile(UMparam.SaveDir,'MatchTable.mat'));

writetable(MatchTable,fullfile(UMparam.SaveDir,'MatchTable.csv'))

% idx = MatchTable.RecSes1 ~= MatchTable.RecSes2 & MatchTable.MatchProb > .9;
% MatchTable(idx,:)