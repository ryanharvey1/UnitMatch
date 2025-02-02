%% DEMO UNIT MATCH 
% This is a demo using the FigShareData on 10.6084/m9.figshare.24305758

%% READ ME
% If you do not use the suggested pipeline to extract raw waveforms (e.g. you don't use Neuropixels/SpikeGLX), make
% sure your 'KiloSortPaths' contains a subfolder called 'RawWaveforms'. There
% should be a NPY file for every cluster with the dimensions of
% UMparam.spikeWidth X nRecordingChannels X 2 (1 for each half of a
% recording). This should contain the average waveform (recommended of at
% least 500 spikes) for every recording channel for every half of a
% recording for that cluster.

%% User input: 
UMparam.SaveDir = 'H:\FigShare_UnitMatch\ExampleOutput'; % Recommended to use end this path with \Probe0\IMRO_1\ if more probes/IMRO tables were used or \AllProbes\AllIMRO\ otherwise
UMparam.KSDir = {'H:\FigShare_UnitMatch\Mouse1\2019-11-21\Probe0\1','H:\FigShare_UnitMatch\Mouse1\2019-11-22\Probe0\1'};  % This is a cell array with a path, in the path there should be a subfolder called 'RawWaveforms'. 
% N.B. if you want to use the functional score evaluation of UnitMatch, 'KSDir' should also contain typical 'Kilosort output', (e.g. spike times etc.)

%% N.B. the following user input can also be automatically extracted and prepared/cleaned up using UMparam = ExtractKilosortData(KiloSortPaths, UMparam) for Kilosorted data of SpikeGLX recorded data (see next section);
% UMparam.RawDataPaths = {'\\path\to\firstrecording','\\path\to\secondrecording','\\path\to\nthrecording'};  % This is a cell array with info on where to find the decompressed recording (.cbin files) --> Necessary when you want UnitMatch to do waveform extraction
% UMparam.AllDecompPaths = {'\\path\to\firstrecording','\\path\to\secondrecording','\\path\to\nthrecording'};  % This is a cell array with info on where to find the decompressed recording (.bin files) --> Necessary when you want UnitMatch to do waveform extraction
% UMparam.AllChannelPos = {[RecordingSites_Recording1],[RecordingSites_Recording2]}; % These are coordinates of every recording channel on the probe (e.g. nRecordingChannels x 2)
% clusinfo = struct; % Note, this can be kilosort input, 
% - clusinfo (this is a struct that contains per unit the following information):
% * cluster_id (e.g. kilosort output clus_id)
% * Good_ID: ones for units that should be included in the analysis
% * RecSesID: Recording Session ID
% * Probe: Which probe (if just 1, ones of numel cluster_id)
% * Depth: depth on probe (optional)
% * Shank: Which shank (optional)
% * Coordinates: Typically 3D Allen Common Coordinate framework coordinates per unit (optional)


% N.B. clusinfo can also be automatically extracted using clusinfo =
% getClusinfo

%% Add paths and subpaths - this will only work when running script in sequence, otherwise just manually add all UnitMatch subfolders
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
Components = strsplit(mfilePath,filesep);
addpath(genpath(fullfile(Components{1:end-1})));

%% Optional (for Kilosort + SpikeGLX users) --- see ExampleAnalysisPipelines for more detail!!
UMparam = ExtractKilosortData(UMparam.KSDir, UMparam); % Extract KS data and do some noise removal, optionally decompresses cbin to bin data and uses BOMBCELL quality metric to define good single units
clusinfo = getClusinfo(UMparam.KSDir); % prepare clusinfo struct

%% Load default parameters
UMparam = DefaultParametersUnitMatch(UMparam);

%% UnitMatch algorithm:
[UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);
if UMparam.AssignUniqueID
   [UniqueIDConversion, MatchTable] = AssignUniqueID(UMparam.SaveDir);
end

%%% N.B. From here it is all evaluation, you don't need this to use UnitMatch
%%% results in your further analysis
%% Automatic evaluation:
EvaluatingUnitMatch(UMparam.SaveDir); % Within session cross-validation
QualityMetricsROCs(UMparam.SaveDir); % Only works in combination with BOMBCELL (and is added to path!!)
ComputeFunctionalScores(UMparam.SaveDir) % Only works when having access to Kilosort output (e.g. spike times etc.) 

%% Curation:
if UMparam.MakePlotsOfPairs
    DrawPairsUnitMatch(UMparam.SaveDir);
    if UMparam.GUI
        FigureFlick(UMparam.SaveDir)
        pause
    end
end

%% Visualization
PlotUnitsOnProbe(clusinfo,UMparam,UniqueIDConversion,WaveformInfo)
