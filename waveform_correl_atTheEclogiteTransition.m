clear all
close all

%set parameters
loadWaveforms = true;
loadFromPreviousState = false;
saveState = false;
loadCorrelationsFromPreviousState = false;
reloadFM3Darrivals = false;
deconvolveSTFs = false;
doProcessStationGathers = true;
doProcessEventGathers = true;
stackTraces = false;
plotSpectrograms = false;
targetSamplingRateForProcessing = 100;
targetSamplingRateAfterFilter = 30;

fm3Dpath = ['/Volumes/nasdata2/Documents2/Greece_MWcluster/FastMarchin',...
    'g/Runs/EventRuns_0.02deg_DepthCorrectlyInverted_2172Events/'];

% fm3DarrivalMatFile = 'arrivals_0.02degGrid_2172events_forEclogiteEvents.mat';
fm3DarrivalMatFile = 'arrivals_0.02degGrid_2172events_relevantStations.mat';
% cstation = {'IDHR','S013','S014','S015','S016','ERM',...
%      'LTK','DID','VLY','LKR','DSF'};
% central/sputhern section:

catalogFile = 'eclogiteEvents__2018-11-09.mat';
catalogFile_withoutWav = '';
corObjFile = 'eclogiteEvents_CorrelationObject_189_2018-11-09.mat';
corEventObjFile = 'eclogiteEvents_withoutWaveforms189_2018-11-09.mat';
corEventObjFileForSTF = '';

%south-central part (Attica, E Peloponnese)
cstation = {'S013','S014','S015','S016','S040','KRND',...
    'DID','VLY','ATHU','DION','PTL','PSAR','VILL','PROD',...
    'AXAR','EPID','LOUT','THAL','EREA','LTK','GUR','KLV','IDHR'};
% 'DIDY','ERM','ACOR',
%northern part:
% cstation = {'AGG','MAKR','AXAR','ATAL','LKR','DSF','MRKA','THL','KALE',...
%     'XOR','NEO','EVR'};
 
cstation = {'S002','PE06','S003','PE03','S004','S005','S006','S007','S008','S009',...
    'S010','PE05','S011','S012','S013','PE02','S014','S015','S016','S017','S018','S019',...
    'S020','S021','AT03','S022','S023','S024','S025','S026','S027','AT02','S028','S029',...
    'S030','S031','S032','S033','S034','S035','S036','S037','S038','S039',...
    'S040','S041','S042','S104','S124','S126',...
    'PE01','PE07','PE08','PE09','PE10','PE11',...
    'TRI','ERM','LKR','ITM',...
    'KLV','LTK','VLX','VLY','AIOA','DRO','PTL','LAKA','VLI','DION','KALE','ANX'}';

corcomp = {'Z','E','N'};
minFrequency = 1.5;
maxFrequency = 10;
tic

secPerDay = 60*60*24;
dayPerSec = 1/secPerDay;

%% read seisan databases with SIMUL & DD-relocations, waveforms etc.

SEISAN_TOP = '/Volumes/Stilton/SEI';
seisanREAdir{1} = fullfile(SEISAN_TOP, 'REA', 'MEDUT');
seisanREAdir{2} = fullfile(SEISAN_TOP, 'REA', 'MEDET');
seisanREAdir{3} = fullfile(SEISAN_TOP, 'REA', 'MEDAU');
seisanREAdir{4} = fullfile(SEISAN_TOP, 'REA', 'GRNW7');
seisanREAdir{5} = fullfile(SEISAN_TOP, 'REA', 'MEDNC');
seisanREAdir{6} = fullfile(SEISAN_TOP, 'REA', 'MEDSC');

simulTOP = '/Volumes/nasdata/Greece_MWcluster/Relocation3D/';
simulrRELOCdir{1} = fullfile(simulTOP, 'MEDUT_Reloc');
simulrRELOCdir{2} = fullfile(simulTOP, 'MEDET_Reloc');
simulrRELOCdir{3} = fullfile(simulTOP, 'MEDAU_Reloc');
simulrRELOCdir{4} = fullfile(simulTOP, 'GRNW_Reloc_new');
simulrRELOCdir{5} = fullfile(simulTOP, 'MEDNC_Reloc');
simulrRELOCdir{6} = fullfile(simulTOP, 'MEDSC_Reloc');
seisanWAVdir = fullfile('/Volumes/nasdata/SEISAN/WAV');

detTimeWindow = seconds(4);

loadDD = true;
DDfolder = ['/Volumes/nasdata/Greece_MWcluster/Relocation_hypoDD/',...
    '2174Events/StatisticRuns/20kmCT_15kmCC_Damp600/hypoDD_setup'];
loadJackknifeErrorsFromFile = true;

startTime = '2007/04/20 00:00:00'; 
endTime = '2007/04/30 00:00:00';

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
    if saveState
        save(['eclogiteEvents_', '_', datestr(datetime('today'),...
            'yyyy-mm-dd'),'.mat'],'greekevents','-v7.3'); %,'-v7.3'
    end
end
%% 
if reloadFM3Darrivals
    arrivals = loadFM3Darrivals(fm3Dpath, [1:1:greekevents.numberOfEvents]);
else
    load(fm3DarrivalMatFile);
end

%% find earthquakes in subclusters
% minX = 80;
% minY = -130;
% maxY = 130;

maxMag = 4;

minX = 80;
minY = -90;
%maxY = 130;
maxY = 70;


% if contains(cstation,'AGG')
%     minY = 70;
%     maxY = 200;
% end

minDepth = 85;
deepMinDistFromSlabTop = -1;

minDepthForDFST = 60;
shallowMinDistFromSlabTop = -10;

% earthquakes within the blueschist-eclogite transition domain, i.e. where
% the slab is deeper than ~70 km
eventI = find(...
    greekevents.table.mag <= maxMag &...
    greekevents.table.x >= minX &...
    greekevents.table.y >= minY &...
    greekevents.table.y < maxY &...
    ((greekevents.table.depth >= minDepth &...
    greekevents.table.DistFromSlabTop <= deepMinDistFromSlabTop) |...
    (greekevents.table.depth >= minDepthForDFST &...
    greekevents.table.DistFromSlabTop <= shallowMinDistFromSlabTop)));

%all:
eventI = sort(eventI);

nSelectedEvents = length(eventI);

%% now get me all waveforms for the events that I want to compare into one
% structure of correlatable events
%that's the structure c: c.(component).(stationname) is a correlation object
tic
if loadCorrelationsFromPreviousState
    load(corObjFile, 'c');
    cstation = fieldnames(c.Z)';
else
    % CHeck whether greekevents still contains waveforms - if there are
    % none, reload the whole database from .mat-file
    if all(cellfun(@isempty,greekevents.waveforms))
        clear greekevents;
        load(catalogFile, 'greekevents');
    end
    c = buildCorrelationCatalog(greekevents, eventI, cstation, corcomp,...
        targetSamplingRateForProcessing);
    if saveState
        save(['eclogiteEvents_CorrelationObject_', num2str(length(eventI)),...
            '_', datestr(datetime('today'),'yyyy-mm-dd'), '.mat'], 'c');
    end
end
toc

%% MAKE some FIGUREs 
% Show 3D model:
% [model3DFig] = plotEarthquakeCatalog3D(greekevents,eventI);

% Make another figure on top of GRT images:
% plotGRTfigureWithEarthquakes(greekevents,-130,70,-200,200,false,false)
% plotGRTfigureWithEarthquakes(subset(greekevents,eventI),minX,300,minY,...
%     maxY,false,false)

%% Event-based correlation objects
eventStations = {'S002','PE06','S003','PE03','S004','S005','S006','S007','S008','S009',...
    'S010','PE05','S011','S012','S013','PE02','S014','S015','S016','S017','S018','S019',...
    'S020','S021','AT03','S022','S023','S024','S025','S026','S027','AT02','S028','S029',...
    'S030','S031','S032','S033','S034','S035','S036','S037','S038','S039',...
    'S040','S041','S042','S104','S124','S126',...
    'PE01','PE07','PE08','PE09','PE10','PE11',...
    'TRI','ERM','LKR',...
    'KLV','LTK','VLX','VLY','AIOA','DRO','PTL','LAKA','VLI','DION','KALE','ANX'}';

stfStackStations = {'S009',...
    'S010','PE05','S011','S012','S013','PE02','S014','S015','S016','S018','S019',...
    'S020','S021','AT03','S022','S023','S024','S025','S026','S027','AT02','S028','S029',...
    'S030','S031','S032','S033','S034','S035','S036','S037','S038','S039',...
    'S040','S041','S042','S104','S124','S126',...
    'PE01','PE07',...
    'KLV','LTK','VLX','VLY','AIOA','DRO','PTL','LAKA','VLI','DION','KALE','ANX'}';

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
            save(['eclogiteEvents_EventCorrelationObject_',...
                num2str(length(eventI)),...
                '_', datestr(datetime('today'),'yyyy-mm-dd'), '.mat'],...
                'eventC');
            save(['eclogiteEvents_EventCorrelationObjectForSTFs_',...
                num2str(length(eventI)),...
                '_', datestr(datetime('today'),'yyyy-mm-dd'), '.mat'],...
                'eventCforSTF');
        end
        save(['eclogiteEvents_withoutWaveforms', num2str(length(eventI)),...
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
    minFreq = 1.5; maxFreq = 10; minSNR = 2;
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
    plotWavesOfEventSortedBy(eventC4, 1, 'Zp', 'X', 'bwig',...
        'markArrivalTimes',true, 'labelarrivals',true,'arrivals',...
        arrivals, 'plotEnvelope', true)
    % plotWavesOfEventSortedBy(eventC4, 1, 'Z', 'X', 'wigbyy','scale',5)
end

%% Plot selected events with three channels etc.
eventNums = 1;
plotThreeChannelEvent = true;
printFigure = true;
plotEnvelope = false;
baseScale = 1.2;
fileNameAddition0 = '';
p = [1 2 3 5];

if ~exist('eventCOut','var')
    eventCOut = resampleEventCorrObject(eventC4, targetSamplingRateAfterFilter);
end
if plotThreeChannelEvent
    plot_and_combine_ChannelsPerEvent(eventCOut, eventNums, p,...
        printFigure, plotEnvelope, baseScale, fileNameAddition0, arrivals)
end

%% now process the waveforms per station and channel

tic
comp = {'Z','N','E'};

if loadCorrelationsFromPreviousState
    load(corObjFileProcessed,'c2');
else
    c2 = processStationGathers(c, cstation, comp, minFrequency,...
        maxFrequency);
    if saveState
        save(['eclogiteEvents_CorrelationObject_PreProcessed', num2str(...
            length(eventI)), '_', datestr(datetime('today'),...
            'yyyy-mm-dd'), '.mat'], 'c2');
    end
end
toc

%% Create threecomponent-objects, rotate and apply polarizationfilter and 
% First, throw out waveforms with bad Signal-to-noise-ratio

if doProcessStationGathers
    minimumSNR = 2.5;
    comp = {'Z','N','E'};
    c3 = requestMinimumSNR(c2, 2.5);
    %c3 = resampleNetworkCorrObject(c3, targetSamplingRateAfterFilter);
    % adjust the trigger times of all traces

    % load('c3_CorrelationObject_PolarizationFiltered_aligned483_2018-08-02.mat');

    c2station = cstation;
    % c2station = {'S011'};
    tic
    % polarizationfilter...: (zrt, 1, 1, 2, 0.02, 0.6) n, J, K, dt, width:
    % these values seem reasonable
    % c4 = threeComponentProcessing(c3, c2station, 1, 1, 2, 0.02, 0.6, true);
    % n =0.5 or 1 changes very little
    % for rather low frequency filters, 1.5 -6 Hz:
    %c4 = threeComponentProcessing(c3, c2station, 0.5, 1, 2, 0.01, 1.0, true);
    % the width is crucial! chose e.g. 0.5 s width for 1.5 - 10 Hz filter
    % let's try 0.35 s window length for 1.5 - 15 Hz filter
    windowLength = 10/maxFrequency * 0.5;
    
    dt = 1/targetSamplingRateForProcessing;
    width = 0.5;
    c4 = threeComponentProcessing(c3, c2station, 0.5, 1, 2, dt, width, true);
    % c4 = threeComponentProcessing(c3, c2station, 0.5, 1, 2, 0.01,...
    %     windowLength, true);

    toc

    % plotWavesAtStationSortedBy(c3,'Zg','S011','distfromslabtop','wig')
    % plotGRTfigureWithEarthquakes(subset(greekevents,eventI),-100,100,true)
    if saveState
        save(['eclogiteEvents_CorrelationObject_3compProcessed', num2str(...
            length(eventI)), '_', datestr(datetime('today'),...
            'yyyy-mm-dd'), '.mat'], 'c3');
    end
end

%% make nice printable plots for correlation objects

conversionStations = cstation;
plotConversionStations = true;
printFigure = true;
plotEnvelope = true;
baseScale = 1.2;
fileNameAddition0 = '';

%resample to reduce filesize
if ~exist('cOut','var')
    cOut = resampleNetworkCorrObject(c4, targetSamplingRateAfterFilter);
end

if plotEnvelope
    fileNameAddition0 = [fileNameAddition0, 'Envelope_'];
end

if plotConversionStations
    for p=1:1:5
        % Plot the three-component figures 3 times: 
        % 1. with raw ZRT channels,
        % 2. with polarization-filtered ZRT channels,
        % 3. with polarization-filtered ZRT channels plus theoretical arrivals
        % 4. with polarization-filtered ZRT channels, plus arrivals and labels
        switch p
            case 1
                plotArrivals = false;
                labelArrivals = false;
                plotComp = {'Z','R','T'};
                fileNameAddition0 = [fileNameAddition0, ''];
                plotEnv = false;
                plotType = 'wig'; linewidth = 0.5;
            case 2
                plotArrivals = true; labelArrivals = false;
                plotComp = {'Zp','Rp','Tp'};
                fileNameAddition0 = [fileNameAddition0, ''];
                plotEnv = false;
                plotType = 'wigbyy'; linewidth = 0.5; scale = 1;
            case 3
                plotArrivals = false; labelArrivals = false;
                plotComp = {'Zp','Rp','Tp'};
                fileNameAddition0 = [fileNameAddition0, ''];
                plotEnv = plotEnvelope;
                plotType = 'bwig'; linewidth = 0.1;
            case 4
                plotArrivals = true; labelArrivals = false;
                plotComp = {'Zp','Rp','Tp'};
                fileNameAddition0 = [fileNameAddition0, ''];
                plotEnv = plotEnvelope;
                plotType = 'bwig'; linewidth = 0.1;
            case 5
                plotArrivals = true;
                labelArrivals = true;
                plotComp = {'Zp','Rp','Tp'};
                fileNameAddition0 = [fileNameAddition0, ''];
                plotEnv = plotEnvelope;
                plotType = 'bwig'; linewidth = 0.1;
        end
        if plotEnv
            scale = baseScale * 1.65;
        else
            scale = baseScale;
        end


        for j=1:1:length(conversionStations)
            for k=1:1:numel(plotComp)
                %plot(c2.(plotComp).(conversionStations{j}).corr,'wig',1.5)
                plotWavesAtStationSortedBy(cOut, plotComp{k},...
                    conversionStations{j}, 'depth', plotType,...
                    'scale', scale, 'linewidth', linewidth,...
                    'markarrivaltimes', plotArrivals, 'arrivals',...
                    arrivals, 'maxconversions', 1, 'maxreflections', 1,...
                    'maxPhaseLength', 4, 'labelarrivals',labelArrivals,...
                    'plotEnvelope', plotEnv)

                xlim([-2 28]);
                oldPos=get(gcf,'Position');
                set(gcf,'Position',[1+(j-1)*230, oldPos(2), 350, oldPos(4)])
                ax(k) = gca();
            end

            if printFigure
                formatThreeComponentWaveformFigures(ax,conversionStations{j},...
                    plotComp, plotType, plotArrivals, labelArrivals, ...
                    fileNameAddition0);
                close all;
            end
        end
    end
end
