clear all
close all

% Check which user is running the script
user_name = char(java.lang.System.getProperty('user.name'));

% Set up python environment
%if strcmp(user_name, "felix")  % Felix, Mac
%    pyenv(Version="/Volumes/MacHD1/Users/felix/Software/miniconda/miniconda3/envs/GMTpy312/bin/python")
%end

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
catalogFile = 'events_2025-04-07.mat';
catalogFile_withoutWav = 'greekevents_629_withoutWaveforms_2018-10-24.mat';

corEventObjFile = 'events_CorrelationObject_7_2025-04-07.mat';
entObjFileForSTF = 'InterfaceEvents_EventCorrelationObjectForSTFs_232_2018-10-24.mat';

corObjFile = 'events_CorrelationObject_468_2018-11-07.mat';
corObjFileProcessed = 'events_CorrelationObject_PreProcessed468_2018-11-21.mat';

% Stations to read in from Seisan files:
requestStations = {'FDF','ANWB'}';
%requestStations = {'OS1', 'OS2', 'OS3', 'OS4','OS5B'}'; %200706
%requestStations = {'F00A','F01A','F02A','F03A','F04A'}'; %200701


% List of stations for correlation analysis
% cstation = {'FDF','ANWB'}';
cstation = {'OS1','OS2','OS3','OS4','OS5B'}'; %200706
%cstation = {'F00A','F01A','F02A','F03A','F04A'}'; %200701
    

% corcomp = {'Z','E','N'};  % Components to be correlated
corcomp = {'Z','one','two'};
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
%seisanWAVdir = fullfile('converted_mseed_test/corrected_waveforms');

detTimeWindow=seconds(4);

startTime = '2007/06/01 00:00:00'; 
endTime = '2007/07/01 00:00:00';

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


%% 
if reloadFM3Darrivals
    arrivals = loadFM3Darrivals(fm3Dpath, [1:1:events.numberOfEvents]);
else
    load('arrivals_0.02degGrid.mat', 'arrivals');
end

%% select earthquakes in specific regions / clusters

% !!Test Lesser Antilles!!
%get earthquakes within a region around the GRT cross section image
%interface fault cluster
I_MWCluster = find(events.table.lon > -62 & events.table.lon < -60 &...
                   events.table.lat > 14 & events.table.lat < 15.5 &...
                   events.table.depth > 35 & events.table.depth < 60 );

eventI = [I_MWCluster]';
eventI = sort(eventI);

% all:
eventI=[1:1:events.numberOfEvents];

nSelectedEvents = length(eventI);

%% now get me all waveforms for the events that I want to compare into one
% structure of correlatable events
% that's the structure c: c.(component).(stationname) is a correlation object
tic
if loadCorrelationsFromPreviousState
    load(corObjFile, 'c');
    cstation = fieldnames(c.Z)';
else
    % CHeck whether events still contains waveforms - if there are
    % none, reload the whole database from .mat-file
    if all(cellfun(@isempty, events.waveforms))
        clear events;
        load(catalogFile, 'events');
    end
    c = buildCorrelationCatalog(events, eventI, cstation, corcomp, ...
                                targetSamplingRateForProcessing);
    if saveState
        %save(['events_CorrelationObject_', num2str(length(eventI)),...
        %    '_', datestr(datetime('today'),'yyyy-mm-dd'), '.mat'],...
        %    'c','-v7.3');
        save(['events_CorrelationObject_', num2str(length(eventI)),...
            '_', datestr(datetime('today'),'yyyy-mm-dd'), '.mat'], 'c');
    end
end
toc

%% MAKE some FIGUREs 
% Show 3D model:
% [model3DFig] = plotEarthquakeCatalog3D(events,eventI);

% Make another figure on top of GRT images:
% plotGRTfigureWithEarthquakes(events,65,125,minY,maxY,false,true)

%% Event-based correlation objects
eventStations = {'FDF','ANWB'}';


% reduce to those that are actually in the dataset
eventStations = eventStations(contains(eventStations, requestStations));

notSTFStackStations = {'FDF','ANWB'};

% reduce to those that are actually in the dataset
% stfStackStations = stfStackStations(contains(stfStackStations, requestStations));
stfStackStations = eventStations(~contains(eventStations, notSTFStackStations));

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
        % CHeck whether events still contains waveforms - if there are
        % none, reload the whole database from .mat-file
        if all(cellfun(@isempty,events.waveforms))
            clear events;
            load(catalogFile,'events');
        end
        % Build the event-correlation object catalog!
        eventC = buildEventCorrelationCatalog(events, eventI, ...
            eventStations,corcomp, targetSamplingRateForProcessing);
        eventCforSTF = buildEventCorrelationCatalog(events, eventI, ...
            stfStackStations,corcomp, targetSamplingRateForProcessing);
    end
    toc
end

% Load waveform-less catalog / empty all waveforms from the catalog to save memory
if loadFromPreviousState
    clear events
    load(catalogFile_withoutWav, 'events');
else
    % remove all waveforms from events
    for j=1:1:events.numberOfEvents
        events.waveforms{j} = {};
    end
    if saveState
        if deconvolveSTFs || doProcessEventGathers
            save(['events_EventCorrelationObject_',...
                  num2str(length(eventI)), '_', datestr(datetime('today'),...
                'yyyy-mm-dd'), '.mat'], 'eventC', '-v7.3');
            save(['events_EventCorrelationObjectForSTFs_',...
                  num2str(length(eventI)), '_', datestr(datetime('today'),...
                'yyyy-mm-dd'), '.mat'], 'eventCforSTF');
        end
        save(['events_', num2str(length(...
            events.numberOfEvents)), '_withoutWaveforms_',...
            datestr(datetime('today'),'yyyy-mm-dd'), '.mat'],...
            'events');
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

    workcatalog = subset(events,eventI);
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
        width, adjustAlignment, false);
    % throw out events with less than XXX usable waveforms

% interesting events:
% plot(eventC3.Z{179})
% for events with magnitude above ca. 2, losing the low frequencies (below
% 1.5) is maybe not great when band-pass filtering. Causes more acausal
% ringing, sometimes worse alignment and less clear secondary arrivals

    plotWavesOfEventSortedBy(eventC4, 179, 'Z', 'X', 'wig', 'scale', 1.2)
    plotWavesOfEventSortedBy(eventC4, 179, 'Z', 'X', 'wigbyy')
end

%% Plot selected events with three channels etc.
eventNums = 1;
plotThreeChannelEvent = true;
printFigure = true;
plotEnvelope = false;
baseScale = 1.2;
fileNameAddition0 = '';

if plotThreeChannelEvent && doProcessEventGathers
    plot_and_combine_ChannelsPerEvent(eventC4, eventNums, printFigure,...
        plotEnvelope, baseScale, fileNameAddition0, arrivals);
end

%% now process the waveforms per station and channel

tic
% comp = {'Z','N','E'};
% comp = {'Z','one','two'};
comp = corcomp;

if loadCorrelationsFromPreviousState
    load(corObjFileProcessed, 'c2');
else
    c2 = processStationGathers(c, cstation, comp, minFrequency,...
        maxFrequency);
    if saveState
        save(['events_CorrelationObject_PreProcessed', num2str(...
            length(eventI)), '_', datestr(datetime('today'),...
            'yyyy-mm-dd'), '.mat'], 'c2');
    end
end
toc

%% Create threecomponent-objects, rotate and apply polarizationfilter and 
% First, throw out waveforms with bad Signal-to-noise-ratio

if doProcessStationGathers
    minimumSNR = 2.5;
    % comp = {'Z','N','E'};
    comp = corcomp;
    c3 = requestMinimumSNR(c2, 2.5);
    %c3 = resampleNetworkCorrObject(c3, targetSamplingRateAfterFilter);
    % adjust the trigger times of all traces

    % load('c3_CorrelationObject_PolarizationFiltered_aligned483_2018-08-02.mat');

    c2station = cstation;
    tic
    % polarizationfilter...: (zrt, 1, 1, 2, 0.02, 0.6) n, J, K, dt, width:
    % these values seem reasonable
    % c4 = threeComponentProcessing(c3, c2station, 1, 1, 2, 0.02, 0.6, true, false);
    % n =0.5 or 1 changes very little
    % for rather low frequency filters, 1.5 -6 Hz:
    %c4 = threeComponentProcessing(c3, c2station, 0.5, 1, 2, 0.01, 1.0, true, false);
    % the width is crucial! chose e.g. 0.5 s width for 1.5 - 10 Hz filter
    % let's try 0.35 s window length for 1.5 - 15 Hz filter
    windowLength = 10/maxFrequency * 0.5;
    
    dt = 1/targetSamplingRateForProcessing;
    width = 0.5;
    c4 = threeComponentProcessing(c3, c2station, 0.5, 1, 2, dt, width,...
        true, false, false);
    c4 = removeEmptyTraces(c4);
    % c4 = threeComponentProcessing(c3, c2station, 0.5, 1, 2, 0.01,...
    %     windowLength, true, false);
    toc
    
    c5 = componentAgc(c4, {'Zp','Rp','Tp'}, 2); % or 1.5 s length?
    %c5 = agcEachTrace(c4, 1);

    % plotWavesAtStationSortedBy(c3,'Zg','S011','distfromslabtop','wig')
    % plotGRTfigureWithEarthquakes(subset(events,eventI),-100,100,true)
    if saveState
        save(['events_CorrelationObject_3compProcessed', num2str(...
            length(eventI)), '_', datestr(datetime('today'),...
            'yyyy-mm-dd'), '.mat'], 'c3');
    end
end



%% Figure equal to stacked figure, but without stacking 
%         baseScale = 1.2;
%         [plotArrivals, labelArrivals, plotComp, fileNameAddition0,plotEnv,...
%             plotType, linewidth, scale, doPlotHistograms] =...
%                 getThreeCompPlotOptions(1,'nonstack_', false, baseScale);
%         %for nc=1:1:numel(plotComp)
%         plotWavesAtStationSortedBy(c5, 'Z', sstation{k}, 'DFST',...
%             plotType, 'scale', scale, 'linewidth', linewidth, 'markarrivaltimes',...
%             plotArrivals, 'arrivals', arrivals, 'maxconversions', 1, 'maxreflections',...
%             1, 'maxPhaseLength', 5, 'avoidPtoS', false, 'labelarrivals',labelArrivals,...
%             'plotEnvelope', plotEnv)
%         ax = gca;
%         newax = formatThreeComponentWaveformFigures(ax,...
%                 sstation{k}, {'Z'}, plotType,...
%                 plotArrivals, plotEnvelope, labelArrivals,...
%                 fileNameAddition0);
%         % end
%         
%         
%         [plotArrivals, labelArrivals, plotComp, fileNameAddition0,plotEnv,...
%             plotType, linewidth, scale, doPlotHistograms] =...
%                 getThreeCompPlotOptions(5,'nonstack_', true, baseScale);
%         % for nc=1:1:numel(plotComp)
%         plotWavesAtStationSortedBy(c5, 'Zp', sstation{k}, 'DFST',...
%             plotType, 'scale', scale, 'linewidth', linewidth, 'markarrivaltimes',...
%             plotArrivals, 'arrivals', arrivals, 'maxconversions', 1, 'maxreflections',...
%             1, 'maxPhaseLength', 5, 'avoidPtoS', false, 'labelarrivals',labelArrivals,...
%             'plotEnvelope', plotEnv, 'colorful', true)
%         
%         ax = gca;
%         newax = formatThreeComponentWaveformFigures(ax,...
%                 sstation{k}, {'Zp'}, plotType,...
%                 plotArrivals, plotEnvelope, labelArrivals,...
%                 fileNameAddition0);

%% 
baseScale = 1.2;
[plotArrivals, labelArrivals, plotComp, fileNameAddition0,plotEnv,...
    plotType, linewidth, scale, doPlotHistograms] =...
        getThreeCompPlotOptions(1,'proc_raw_', false, baseScale);
plotComp = {'Z', 'N', 'E'};
clear axes
for nc=1:1:numel(plotComp)
    plotWavesAtStationSortedBy(c, plotComp{nc}, 'S010', 'DFST',...
        'wig', 'scale', scale, 'linewidth', linewidth, 'markarrivaltimes',...
        plotArrivals, 'arrivals', arrivals, 'maxconversions', 1, 'maxreflections',...
        1, 'maxPhaseLength', 5, 'avoidPtoS', false, 'labelarrivals',labelArrivals,...
        'plotEnvelope', plotEnv)
    axes(nc) = gca;
    axes(nc).YLim = [345 360];
    axes(nc).XLim = [-2 16];
end
newax = formatThreeComponentWaveformFigures(axes, 'S010', plotComp,...
        plotType, false, false, false, fileNameAddition0);


%% plot waveforms with 
% plotWavesAtStationSortedBy(c4,'Z','S011','method','Z','plottype','wig',...
%    'markArrivalTimes',true)

%% Plot waveforms at a station sorted by some parameter of the event (depth, lat, x-)
% Stations S011 and S012 and S013 ( this has maybe lowest SNR, and weaker 
% conv/refl compared to P- and S) seem to be best suitable for researching the
% reflections/conversions between P- and S- arrival. Also 
% comp = 'Zp';
% sta = 'S010';

%sort by lat

%sort by distance from trench events.table.x
%- this sorting might be preferrable, conversion waveforms align quite well

%sort by distance along strike (along trench) events.table.x

%sort by depth events.table.z

%sort by distance from slab Top
% plotWavesAtStationSortedBy(c4,'Zg','S011','distfromslabtop','wig')
% plot(c2.(comp).(sta).corr,'wig',0.05,traceOrder,...
%    c2.(comp).(sta).cat.table.x,'linewidth',0.1);


%% make nice printable plots for correlation objects

% plotComp = {'Zp'};

conversionStations = {'S009','S010','S011','PE05','PE07','S012',...
    'PE02','S013','S014','S015','S016'};
% conversionStations = {'S009','S010','S011','PE05','PE07','S012',...
%     'PE02','S013','S014','S015','S016','S019','S020','S022'};
% conversionStations = {'S015'}; 
% conversionStations = {'S009','S010','S011','S014','S015','PE02','PE07'}; 
% conversionStations = {'S009','S010','S011','PE05','PE07','S012',...
%     'PE02','S013','S014','S015','S016','S018','S019','S020','S022',...
%     'S027','S040'};
plotConversionStations = true;
printFigure = true;
plotEnvelope = true;
baseScale = 1.2;
fileNameAddition0 = 'Whist_';
clear cOut;

if plotEnvelope
    fileNameAddition0 = [fileNameAddition0, 'Envelope_'];
    %resample to reduce filesize
    if ~exist('cOut','var')
        cOut = resampleNetworkCorrObject(c5, targetSamplingRateAfterFilter);
    end
end

if plotConversionStations
    for p=[3,4,5] %1:1:4 1,2:5 1,3,4,5
        % Plot the three-component figures 3 times: 
        % 1. with raw ZRT channels,
        % 2. plotted against proper y-axis
        % 3. with polarization-filtered ZRT channels,
        % 4. with polarization-filtered ZRT channels plus theoretical arrivals
        % 5. with polarization-filtered ZRT channels, plus arrivals and labels
        [plotArrivals, labelArrivals, plotComp, fileNameAddition0,plotEnv,...
            plotType, linewidth, scale, doPlotHistograms] =...
            getThreeCompPlotOptions(p,fileNameAddition0, plotEnvelope, baseScale);
        if plotEnv
            scale = baseScale * 1.65;
            
            %cOut = SimplifyWaveformsOfNetworkCorrObject(cOut, 2);
        else
            scale = baseScale;
        end

        for j=1:1:length(conversionStations)
            clear ax;
            % give each trace higher relative scale when there are more
            % traces (just for the visualization)
            %wavs = get(c4.(plotComp{1}).(conversionStations{j}).corr,...
            %    'waveform');
            %if numel(wavs) < 50
            %    scale = baseScale;
            %else
            %    scale = baseScale * (numel(wavs)/50) ^ (1/8);
            %end
            
            if doPlotHistograms
                [axx, SonRp_vs_PonZp, SonRp_vs_SonTp, SonTp_vs_PonZp] =...
                    plotWavComparisonHistogram(c4, conversionStations{j},...
                    'method','DFST','arrivals',arrivals,'windowLength',...
                    1,'scaling','lin', 'RunMedianWidth', 4);
                if mean(SonRp_vs_PonZp) < 1
                    plotPMS = true;
                    avoidPtoS = false;
                end
            end
            
            
            for k=1:1:numel(plotComp)
                %plot(c2.(plotComp).(conversionStations{j}).corr,'wig',1.5)
                %try
                plotWavesAtStationSortedBy(cOut, plotComp{k},...
                    conversionStations{j}, 'DFST', plotType,...
                    'scale', scale, 'linewidth', linewidth,...
                    'markarrivaltimes', plotArrivals, 'arrivals',...
                    arrivals,'maxconversions',1,'maxreflections', 1,...
                    'maxConvAndRefl', 2, 'maxPhaseLength', 5,...
                    'labelarrivals',labelArrivals, 'plotEnvelope',...
                    plotEnv, 'avoidPtoS', false, 'plotPMS',true,...
                    'colorful',true)
                    %'SV_vs_P',SonRp_vs_PonZp, 'SH_vs_P', SonTp_vs_PonZp,...
                %catch
                %    continue
                %end

                if contains(conversionStations{j},{'S019';'S02';'S040'})
                    xlim([-2 22]);
                elseif contains(conversionStations{j},{'S03','S04'})
                    xlim([-2 28]);
                else
                    xlim([-2 16]);
                end
                oldPos=get(gcf,'Position');
                set(gcf,'Position',[1+(j-1)*230, oldPos(2), 350, oldPos(4)])
                ax(k) = gca();
            end

            if doPlotHistograms
                ax = [ax, axx];
            end

            if printFigure
                %try
                newax = formatThreeComponentWaveformFigures(ax,...
                        conversionStations{j}, plotComp, plotType,...
                        plotArrivals, plotEnvelope, labelArrivals,...
                        fileNameAddition0);
                %catch
                %    continue
                %end
            end
            if length(conversionStations) > 1
                close all;
            end
        end
    end
end

% 
pos = [240   475   374   205];
set(gcf,'Position',pos)

set(gca,'Position',[0.10 0.1 0.89 0.85])
% set(gca,'Position',[0.10 0.05 0.89 0.8])
xlim([-2 16]);
%yl = get(gca,'ylim');
%ylim([round(yl(1)) round(yl(1))+30])
ylim([339 352])
ylim([287 300])
% ylim([277 290])

set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3),pos(4)])
set(gcf,'renderer','painters')


%%
% dataformat=;
% arrivalObject = ARRIVAL.RETRIEVE('seisan', 'param1', _value1_, ...

% test for seisan s-file loading:
% Sfile(fileread(fullfile(sfiles(j).dir, sfiles(j).name)))



%% print nice figure
printFigure=false;
if printFigure
    box on
    hAll = findall(gcf);
    for idx = 1 : length(hAll)
      try
        set(hAll(idx),'LineWidth',0.02);
      catch
        % never mind...
      end
      try
        set(hAll(idx),'fontsize',4);
      catch
        % never mind...
      end
      try
        set(hAll(idx),'Color','k');
      catch
        % never mind...
      end
    end
    set(gcf, 'Color', 'w');
    set(gca, 'Color','w')
    set(gcf,'Units','centimeters');
    set(gcf,'PaperUnits','centimeters');
    pos = get(gcf,'Position');
%     set(gcf,'Position',[pos(1),pos(2),20,pos(4)*20/pos(3)]);
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3),pos(4)])
%     set(gcf,'PaperOrientation','portrait');
%     set(gcf,'PaperType','a4');
	print(gcf,'PE07','-dpdf','-painters','-r300')
%     export_fig('lala.png','-CMYK','-R600','-depsc');
end

