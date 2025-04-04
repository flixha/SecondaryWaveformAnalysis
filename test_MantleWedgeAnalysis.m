clear all
close all

% set parameters
loadFromPreviousState = false;              % activate for quick-loading after saving state
saveState = true;                           % whether to save state at some points in the code
loadCorrelationsFromPreviousState = false;  % whether to quick-load correlation objects
reloadFM3Darrivals = false;                 % whether to reload fm3d arrivals from FM3D output text
                                            %   files (can be slow)
deconvolveSTFs = false;                     % Whether to try empirical source time function
                                            %   deconvolution to make signals sharper
doProcessStationGathers = true;             % whether to process gathers for individual stations
doProcessEventGathers = false;              % whether wo process gathers for individual events
loadWaveforms = true;                       % whether to reload waveforms
targetSamplingRateForProcessing = 100;      % sampling rate for processing (can be high for good
                                            %   better alignment
targetSamplingRateAfterFilter = 30;         % target sampling rate after filtering (can be lower to
                                            %   speed up plotting

% Path to FM3D output folders
fm3Dpath = ['/Volumes/nasdata2/Documents2/Greece_MWcluster/',...
    'FastMarching/Runs/EventRuns_0.02deg_DepthCorrectlyInverted/'];

% earthquake catalog file once saved from Matlab:
catalogFile = 'greekevents_629_2018-10-24.mat';
catalogFile_withoutWav = 'greekevents_629_withoutWaveforms_2018-10-24.mat';
corEventObjFile = 'InterfaceEvents_EventCorrelationObject_232_2018-10-24.mat';
corEventObjFileForSTF = 'InterfaceEvents_EventCorrelationObjectForSTFs_232_2018-10-24.mat';
corObjFile = 'events_CorrelationObject_468_2018-11-07.mat';
corObjFileProcessed = 'events_CorrelationObject_PreProcessed468_2018-11-21.mat';

% Stations to read in from Seisan files:
requestStations = {'FDF','ANWB'}';
%requestStations = {'OS1', 'OS2', 'OS3', 'OS4','OS5B'}'; %200706
%requestStations = {'F00A','F01A','F02A','F03A','F04A'}'; %200701

% List of stations for correlation analysis
cstation = {'FDF','ANWB'}';
%cstation = {'OS1','OS2','OS3','OS4','OS5B'}'; %200706
%cstation = {'F00A','F01A','F02A','F03A','F04A'}'; %200701

corcomp = {'Z','E','N'};  % Components to be correlated
%corcomp = {'Z','1','2'};
minFrequency = 1.5;
maxFrequency = 10;

tic
secPerDay = 60*60*24;
dayPerSec = 1/secPerDay;

%% read seisan databases with event locations and waveforms etc.
SEISAN_TOP='MWS/Seisan-linux';
seisanREAdir{1} = fullfile(SEISAN_TOP, 'REA', 'ANTIP');
%simulTOP = '/Volumes/nasdata/Greece_MWcluster/Relocation3D/';
%simulrRELOCdir{1} = fullfile(simulTOP, 'MEDUT_Reloc');
seisanWAVdir = fullfile('MWS/Seisan-linux/WAV/ANTIP');
%seisanWAVdir = fullfile('corrected_waveforms');
detTimeWindow=seconds(4);
startTime = '2007/01/09 00:00:00'; 
endTime = '2007/01/09 59:59:59';
% endTime = '2007/01/31 00:00:00';
% endTime = '2007/12/31 00:00:00';

if loadFromPreviousState && ~loadCorrelationsFromPreviousState
    load(catalogFile, 'events');
elseif loadFromPreviousState && loadCorrelationsFromPreviousState
    %load waveform-less catalog object
    load(catalogFile_withoutWav, 'events');
elseif ~loadFromPreviousState
    tic
    % events = readSeisanDatabases(SEISAN_TOP, seisanREAdir,...
    %     loadWaveforms, seisanWAVdir, startTime, endTime, requestStations);
    events = readSeisanDatabases(SEISAN_TOP, seisanREAdir,...
        loadWaveforms, seisanWAVdir, startTime, endTime);
    toc
    if saveState
        save(['events_', datestr(datetime('today'),'yyyy-mm-dd'),...
            '.mat'],'events','-v7.3');
    end
end

% Filter out FDF station from the events structure.
if isfield(events, 'traces')
    for i = 1:length(events)
        if isfield(events(i).traces, 'st')
            for j = length(events(i).traces.st):-1:1 % Iterate backwards to safely remove elements
                if strcmp(events(i).traces.st(j).meta.station, 'FDF')
                    events(i).traces.st(j) = []; % remove element.
                end
            end
        end
    end
end