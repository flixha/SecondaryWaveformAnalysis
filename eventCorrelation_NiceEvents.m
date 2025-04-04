clear all
close all

%set parameters
loadWaveforms = true;
loadFromPreviousState = false;
saveState = true;
loadCorrelationsFromPreviousState = false;
reloadFM3Darrivals = false;
deconvolveSTFs = false;
doProcessEventGathers = true;
stackTraces = false;
plotSpectrograms = false;
targetSamplingRateForProcessing = 100;
targetSamplingRateAfterFilter = 30;

fm3Dpath = ['/Volumes/nasdata2/Documents2/Greece_MWcluster/FastMarchin',...
    'g/Runs/EventRuns_0.02deg_DepthCorrectlyInverted_2172Events/'];

fm3DarrivalMatFile = 'arrivals_0.02degGrid_2172events_relevantStations.mat';

savePrefix = 'niceEvents_';
catalogFile = 'niceEvents__2018-11-14.mat';
%catalogFile = eclogiteEvents_withoutWaveforms189_2018-11-09.mat;
catalogFile_withoutWav = 'niceEvents_withoutWaveforms4_2018-11-14.mat';
% catalogFile_withoutWav = 'eclogiteEvents_withoutWaveforms325_2018-10-22.mat';
corObjFile = '';
corEventObjFile = '';
corEventObjFileForSTF = '';

 
cstation = {'S001','S002','PE06','S003','PE03','S004','S005','S006','S007','S008','S009',...
    'S010','PE05','S011','S012','S013','PE02','S014','S015','S016','S017','S018','S019',...
    'S020','S021','AT03','S022','S023','S024','S025','S026','S027','AT02','S028','S029',...
    'S030','S031','S032','S033','S034','S035','S036','S037','S038','S039',...
    'S040','S041','S042','S104','S124','S126',...
    'PE01','PE07','PE08','PE09','PE10','PE11',...
    'TRI','ERM','LKR','ITM'}';

% 'KLV','LTK','VLX','VLY','AIOA','DRO','PTL','LAKA','VLI','DION','KALE','ANX'

corcomp = {'Z','E','N'};
minFrequency = 1.0;
maxFrequency = 10;
tic

secPerDay = 60*60*24;
dayPerSec = 1/secPerDay;

%% read seisan databases with SIMUL & DD-relocations, waveforms etc.

SEISAN_TOP = '/Volumes/Stilton/SEI';
seisanREAdir{1} = fullfile(SEISAN_TOP, 'REA', 'MEDUT');
seisanREAdir{2} = fullfile(SEISAN_TOP, 'REA', 'MEDET');
% seisanREAdir{3} = fullfile(SEISAN_TOP, 'REA', 'MEDAU');
% seisanREAdir{4} = fullfile(SEISAN_TOP, 'REA', 'GRNW7');
% seisanREAdir{4} = fullfile(SEISAN_TOP, 'REA', 'MEDNC');
seisanREAdir{3} = fullfile(SEISAN_TOP, 'REA', 'MEDSC');

simulTOP = '/Volumes/nasdata/Greece_MWcluster/Relocation3D/';
simulrRELOCdir{1} = fullfile(simulTOP, 'MEDUT_Reloc');
simulrRELOCdir{2} = fullfile(simulTOP, 'MEDET_Reloc');
% simulrRELOCdir{3} = fullfile(simulTOP, 'MEDAU_Reloc');
% simulrRELOCdir{4} = fullfile(simulTOP, 'GRNW_Reloc_new');
% simulrRELOCdir{4} = fullfile(simulTOP, 'MEDNC_Reloc');
simulrRELOCdir{3} = fullfile(simulTOP, 'MEDSC_Reloc');
seisanWAVdir = fullfile('/Volumes/nasdata/SEISAN/WAV');

detTimeWindow = seconds(4);

loadDD = true;
DDfolder = ['/Volumes/nasdata/Greece_MWcluster/Relocation_hypoDD/',...
    '2174Events/StatisticRuns/20kmCT_15kmCC_Damp600/hypoDD_setup'];
loadJackknifeErrorsFromFile = true;

startTime = '2006/06/15 00:00:00'; 
endTime = '2007/04/24 00:00:00';

if loadFromPreviousState && ~loadCorrelationsFromPreviousState
    load(catalogFile, 'greekevents');
elseif loadFromPreviousState && loadCorrelationsFromPreviousState
    %load waveform-less catalog object
    load(catalogFile_withoutWav, 'greekevents');
elseif ~loadFromPreviousState
    tic
    greekevents = readMEDUSAseisanDatabases(SEISAN_TOP, seisanREAdir,...
        loadWaveforms, seisanWAVdir, startTime, endTime, simulrRELOCdir,...
        loadDD, DDfolder, loadJackknifeErrorsFromFile, detTimeWindow,...
        cstation);
    toc
end
%% 
if reloadFM3Darrivals
    arrivals = loadFM3Darrivals(fm3Dpath, [1:1:greekevents.numberOfEvents]);
else
    load(fm3DarrivalMatFile);
end

if greekevents.numberOfEvents < numel(arrivals)
    % then I need to correct the EventIDs of all events in greekevents
    gev2 = load('eclogiteEvents_withoutWaveforms189_2018-11-09.mat',...
        'greekevents');
    gev2 = gev2.greekevents;
    for j=1:1:greekevents.numberOfEvents
        [minTimeDiff, minIdx] = min(abs(greekevents.table.otime(j) -...
            gev2.table.otime));
        greekevents.table.EventID(j) = gev2.'.EventID(minIdx);
    end
    clear gev2;
end


if ~loadFromPreviousState && saveState
    save([savePrefix, '_', datestr(datetime('today'),...
        'yyyy-mm-dd'),'.mat'],'greekevents','-v7.3'); %,'-v7.3'
end
    
%% find earthquakes in subclusters

% selected high magnitude/quality-events:
% eventI = [276, 313, 10, ...
%     4, 12, 144,...
%     232, 468, 387,...
%     337, 120 , 225];

%maybe MW: 514, 481, IS: 566, 538, deeper: 164, 582
% MW, IF, IC, ISM selection:
% eventI = [10, 4, 127, 164];

% MW: 2006_06_19    11:31:49.2 or 20061007195802 20061120210228
evDates(1) = datenum(2006,06,19,11,31,49.2);
evDates(2) = datenum(2006,10,07,19,58,02);
% IF: 2006_06_15    21:47:45.0 or 20060907235122 20061026152249 
% MEDET: 20070407143255.5
%evDates(3) = datenum(2006,09,07,23,51,22); 
evDates(3) = datenum(2007,04,07,14,32,55.5); 
%evDates(3) = datenum(2006,06,15,21,47,45.0); % this would be a great event, but the FM3D run failed!
evDates(4) = datenum(2006,06,20,18,12,02);
% IC: 2006_08_28    04:22:46.4
% uhm not really.. what about 20061029231610.2 20061203222312.9 
% or 20070104051225.7 or 2006,08,18,04,22,46.4
evDates(5) = datenum(2006,10,29,23,16,10.2);
evDates(6) = datenum(2006,12,03,22,23,12.9);
% ISM: 2006_09_25    04:49:12.0 or 20060617051214.8
evDates(7) = datenum(2006,09,25,04,49,12.0);
evDates(8) = datenum(2006,06,17,05,12,14.8);
% lower plane of the DSZ seismicity: only ONE event:
% 2007-04-23 MEDSC
evDates(9) = datenum(2007,04,23,22,22,33);

nSelectedEvents = length(evDates);
eventI = [];
for j=1:1:nSelectedEvents
    [minTimeDiff, evIDX] = min(abs(greekevents.otime - evDates(j)));
    if minTimeDiff < 5 * dayPerSec
        eventI = [eventI;evIDX];
    end
end



%% Event-based correlation objects
eventStations = {'S002','PE06','S003','PE03','S004','S005','S006','S007','S008','S009',...
    'S010','PE05','S011','S012','S013','PE02','S014','S015','S016','S017','S018','S019',...
    'S020','S021','AT03','S022','S023','S024','S025','S026','S027','S028','S029',...
    'S040','S041','S042','S104','S124','S126',...
    'TRI','ERM','LKR','ITM'}';
% 'S001', 'PE09', 'S031','S032','S033','S030',
% 'S035','S037','S038','S039','S036','AT02',
% 'PE01','PE07','PE08','PE10','PE11'
% 'KLV','LTK','VLX','VLY','AIOA','DRO','PTL','LAKA','VLI','DION','KALE','ANX

stfStackStations = {'S009',...
    'S010','PE05','S011','S012','S013','PE02','S014','S015','S016','S018','S019',...
    'S020','S021','AT03','S022','S023','S024','S025','S026','S027','AT02','S028','S029',...
    'S030','S031','S032','S033','S034','S035','S036','S037','S038','S039',...
    'S040','S041','S042','S104','S124','S126',...
    'PE01','PE07'}';
    
%'KLV','LTK','VLX','VLY','AIOA','DRO','PTL','LAKA','VLI','DION','KALE','ANX'

tic

% Find the source time function of each event and deconvolce all waveforms of that event by it
% 1st strategy:
%     Do a k-means clustering analysis on the first arrival to sort into 2
%     clusters: the positive and the negative first arivals. Then flip the
%     negative ones, and stack them together with the positive ones into
%     one STF estimate based on the first 0... - x seconds.
if loadCorrelationsFromPreviousState
    load(corEventObjFile, 'eventC');
end

if deconvolveSTFs || doProcessEventGathers
    if loadCorrelationsFromPreviousState
        %load(corEventObjFile, 'eventC');
        load(corEventObjFileForSTF, 'eventCforSTF');
        cstation = fieldnames(c.Z)';
    else
        % CHeck whether greekevents still contains waveforms - if there are
        % none, reload the whole database from .mat-file
        if all(cellfun(@isempty,greekevents.waveforms))
            clear greekevents;
            load(catalogFile,'greekevents');
        end
        % Build the event-correlation object catalog!
        eventC = buildEventCorrelationCatalog(greekevents, eventI,...
            eventStations,corcomp, targetSamplingRateForProcessing);
        eventCforSTF = buildEventCorrelationCatalog(greekevents, eventI,...
            stfStackStations,corcomp, targetSamplingRateForProcessing);
    end
    toc
end

% Load waveform-less catalog / empty all waveforms from the catalog to save
% memory
if loadFromPreviousState
    clear greekevents
    load(catalogFile_withoutWav, 'greekevents');
else
    % remove all waveforms from greekevents
    for j=1:1:greekevents.numberOfEvents
        greekevents.waveforms{j} = {};
    end
    if saveState
        if deconvolveSTFs
            save([savePrefix,'EventCorrelationObject_',...
                num2str(length(eventI)),...
                '_', datestr(datetime('today'),'yyyy-mm-dd'), '.mat'],...
                'eventC');
            save([savePrefix,'EventCorrelationObjectForSTFs_',...
                num2str(length(eventI)),...
                '_', datestr(datetime('today'),'yyyy-mm-dd'), '.mat'],...
                'eventCforSTF');
        end
        save([savePrefix, 'withoutWaveforms', num2str(length(eventI)),...
            '_', datestr(datetime('today'),'yyyy-mm-dd'), '.mat'],...
            'greekevents');
    end
end


%% Deconvolve the waveforms by a source time function estimate. First, 
% create that estimate from sorting the waveforms into 2 bins, for positive
% and negative first arrivals (however, at first we don't know which bin is
% which).

if deconvolveSTFs
    %make copies of the event-structures
    % eventC2 = eventCforSTF;
    % eventC3 = eventCforSTF;

    workcatalog = subset(greekevents,eventI);
    deconvI = find(workcatalog.mag >= 2.5);
    minimumChannelsForSTFStack = 10;
    minSNR = 1.5;
    weightBeforeStacking = true;
    doplot = false;
    %j = 28;
    % deconvI = 186;

    STFestimate = estimate_source_time_function(workcatalog, eventCforSTF,...
        deconvI, minimumChannelsForSTFStack, minSNR, weightBeforeStacking,...
        doplot);
end

%% Deconvolve

if deconvolveSTFs
    eventC2 = eventCforSTF;
    % A = get(eventC2.Z{j,1},'waveforms');
    % % deconvolve(eventC2.Z{j,1}, STFestimate)
    % % B = deconvolve(A, STFestimate, 'method', 'iterative','iterations',100,'minimalErrorChange',0.0001,'gaussianWidth',0.2)
    % % B = deconvolve(A, STFestimate, 'method', 'spectraldivision','regularization','wat','tshift',0.2);
    % B = deconvolve(A, STFestimate, 'method', 'timedomain','tshift',0.2);
    % eventC4 = eventC2;
    % eventC4.Z{j,1} = set(eventC4.Z{j,1},'waveform',B);
    % plot(eventC4.Z{j,1},'wig')
end

%% Process Event Objects
if doProcessEventGathers
    % process the event-based correlation objects
    minFreq = 1.5; maxFreq = 10; minSNR = 5;
    %minFreq = 1.5; maxFreq = 10; minSNR = 2;
    eventC3 = processEventCorrelObjects(eventC, minFreq, maxFreq, minSNR);

    % three-components processing!
    n = 0.5; J = 1; K = 2; dt = 1/targetSamplingRateForProcessing;
    width = 0.5; adjustAlignment = true;
    eventC4 = threeComponentProcessingForEvObj(eventC3, n, J, K, dt,...
        width, adjustAlignment);
    % throw out events with less than XXX usable waveforms

    % interesting events:
    % plot(eventC3.Z{179})
    % especialy interesting: 179
    % 1 3 208
    % for event with magnitude above ca. 2, losing the low frequencies (below
    % 1.5) is maybe not great when band-pass filtering. Causes more acausal
    % ringing, sometimes worse alignment and less clear secondary arrivals

    %plotWavesOfEventSortedBy(eventC4, 1, 'Z', 'X', 'wig', 'scale', 1.2)
    %plotWavesOfEventSortedBy(eventC4, 1, 'Zp', 'X', 'wig')
    
%     plotWavesOfEventSortedBy(eventC4, 1, 'Zp', 'X', 'bwig',...
%         'markArrivalTimes',true, 'labelarrivals',true,'arrivals',...
%         arrivals, 'plotEnvelope', true)
    % plotWavesOfEventSortedBy(eventC4, 1, 'Z', 'X', 'wigbyy','scale',5)
end

%% Plot selected events with three channels etc.
eventIDs = greekevents.table.EventID(eventI);
% eventIDs = 69;
plotThreeChannelEvent = true;
printFigure = true;
plotEnvelope = true;
baseScale = 1;
fileNameAddition0 = '';
%
p = [1 2 3 5];

if ~exist('eventCOut','var')
    eventCOut = resampleEventCorrObject(eventC4, targetSamplingRateAfterFilter);
end
if plotThreeChannelEvent
    plot_and_combine_ChannelsPerEvent(eventC4, eventIDs, p,...
        printFigure, plotEnvelope, baseScale, fileNameAddition0, arrivals)
end


