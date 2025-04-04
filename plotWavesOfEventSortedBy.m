function [pickOffset] = plotWavesOfEventSortedBy(eventC, j, comp, varargin)

% varargin: method, plottype, linewidth, scale, 
% plottype is 'sha' or 'wig' 
SECPERDAY = 24*60*60;


p = inputParser;

defaultMethod = 'DISTANCEFROMTRENCH';
validMethods = {'LAT', 'LATITUDE', 'LON', 'LONGITUDE', 'X',...
    'DISTFROMTRENCH', 'DISTANCEFROMTRENCH', 'Y', 'DISTALONGSTRIKE',...
    'DISTANCEFROMEVENT','DISTFROMEVENT','DFE','OFFSET'};
checkMethod = @(x) any(validatestring(x,validMethods));

defaultPlotType = 'WIG';
validPlotType = {'SHA', 'WIG', 'WIGBYY', 'CWIG', 'BWIG','RAW'};
checkPlotType = @(x) any(validatestring(x,validPlotType));

defaultLineWidth = 0.5;
defaultScale = 1;

defaultMarkArrivalTimes = false;

addRequired(p,'eventC',@isstruct);
addRequired(p,'j',@isnumeric);
addRequired(p,'comp',@ischar);
addOptional(p,'method', defaultMethod, checkMethod);
addOptional(p,'plottype', defaultPlotType, checkPlotType);
addParameter(p,'linewidth', defaultLineWidth, @isnumeric)
addParameter(p,'scale', defaultScale, @isnumeric)

defaultArrivals = struct();
addParameter(p,'arrivals', defaultArrivals, @isstruct)
addParameter(p,'markArrivalTimes', defaultMarkArrivalTimes, @islogical)

defaultLabelArrivals = false;
addOptional(p,'labelArrivals', defaultLabelArrivals, @islogical);

defaultMaxConversions = 1;
addOptional(p,'maxConversions', defaultMaxConversions, @isnumeric);

defaultMaxReflections = 1;
addOptional(p,'maxReflections', defaultMaxReflections, @isnumeric);

defaultMaxInteractions = 1;
addOptional(p,'maxInteractions', defaultMaxInteractions, @isnumeric);

defaultMaxPhaseLength = 4;
addParameter(p,'maxPhaseLength', defaultMaxPhaseLength, @isnumeric)

defaultPlotEnvelope = false;
addParameter(p,'plotEnvelope', defaultPlotEnvelope, @islogical)

parse(p, eventC, j, comp, varargin{:})
% tic
% THIS part EATS time. 7 second execution for these simple operations below
% method = p.Results.method;
% plottype = p.Results.plottype;
% linewidth = p.Results.linewidth;
% markArrivalTimes = p.Results.markArrivalTimes;
% arrivals = p.Results.arrivals;
% maxConversions = p.Results.maxConversions;
% labelArrivals = p.Results.labelArrivals;
% maxReflections  = p.Results.maxReflections;
% maxInteractions = p.Results.maxInteractions;
% maxPhaseLength = p.Results.maxPhaseLength;
% toc

% this takes only 0.4 s compared to the above block
% tic
res = p.Results;
method = res.method;
plottype = res.plottype;
plotEnvelope = res.plotEnvelope;
linewidth = res.linewidth;
markArrivalTimes = res.markArrivalTimes;
arrivals = res.arrivals;
maxConversions = res.maxConversions;
labelArrivals = res.labelArrivals;
maxReflections = res.maxReflections;
maxInteractions = res.maxInteractions;
maxPhaseLength = res.maxPhaseLength;
% toc


isStacked = false;

% if ismember('stddepth', eventC.cat.table.Properties.VariableNames)
%     isStacked = true;
% end

%read station coordinates
stationNames = get(get(eventC.Z{j},'waveforms'),'station');
if exist('stations_LLEXYZ.mat', 'file') == 2
    load('stations_LLEXYZ.mat','stations');
else
    stations = readtable('stations_LatLonElev_XYZ.dat');
    stations.Properties.VariableNames = {'name' 'lat' 'lon' 'elev' 'x'...
        'y' 'z'};
end
[requStationID, stationSort] = ismember(stations.name, stationNames);

reqStations = stations(requStationID,:);
stationSortOrder = stationSort(requStationID);

[~, inversSortOrder] = ismember([1:1:numel(stationSortOrder)],...
    stationSortOrder);
reqStations = reqStations(inversSortOrder,:);


% Plot not the waveforms but their envelopes/Hilbert transforms
plusEnv = 0;
if plotEnvelope
    eventC.(comp){j} = hilbert(eventC.(comp){j});
    %w = get(eventC.(comp){j},'waveform');
    % w = power(w,3/4);
    %eventC.(comp){j} = set(eventC.(comp){j},'waveform',w);
    if strcmpi(plottype,'wig') || strcmpi(plottype,'bwig') ||...
            strcmpi(plottype,'cwig') || strcmpi(plottype,'raw')
        plusEnv = -0.5;
    end
end

% if raw, then scale the records across channels
if strcmpi(plottype, 'raw')
    compNames = fieldnames(eventC);
    if contains(comp,'p')
        ci = find(contains(compNames,'p'));
    elseif contains(comp,'R')
        ci = find(contains(compNames,{'Z';'R';'T'}) & ~contains(...
            compNames,'p') & ~contains(compNames,'C'));
    else
        ci = find(contains(compNames,{'Z';'N';'E'}) & ~contains(...
            compNames,'p') & ~contains(compNames,'C'));
    end
    AmpPerComp = zeros(length(get(eventC.(compNames{ci(1)}){j},'waveform')),0);
    for c=1:1:length(ci)
        w = get(eventC.(compNames{ci(c)}){j},'waveform');
        MaxAmpPerComp = [AmpPerComp, max(w)];
    end
    MeanMaxAmp = mean(abs(MaxAmpPerComp),2);
    
    for c=1:1:length(ci)
        w = get(eventC.(compNames{ci(c)}){j},'waveform');
        w = times(w,1./MeanMaxAmp);
        eventC.(comp){j} = set(eventC.(compNames{ci(c)}){j},'waveform',w);
    end
end



switch upper(method)
    case {'LAT','LATITUDE'}
        %sort by lat
        [sortedLats, sortLatIndices] = sort(reqStations.lat);
        traceOrder = sortLatIndices;
        sortedValues = sortedLats;
        titlename = 'latitude';
        ylabelname = titlename;
    case {'LON','LONGITUDE'}
        %sort by lat
        [sortedLons,sortLonIndices] = sort(reqStations.lon);
        traceOrder = sortLonIndices;
        sortedValues = sortedLons;
        titlename = 'longitude';
        ylabelname = titlename;
    case {'X','DISTFROMTRENCH','DISTANCEFROMTRENCH'}
        %sort by distance from trench greekevents.table.x
        %- this sorting might be preferrable, conversion waveforms align quite well
        [sortedXs,sortXIndices] = sort(reqStations.x);
        traceOrder = sortXIndices;
        sortedValues = sortedXs;
        titlename = 'distance along dip';
        ylabelname = [titlename, ' (km), station'];
    case {'Y','DISTALONGSTRIKE'}
        %sort by distance along strike (along trench) greekevents.table.x
        [sortedYs,sortYIndices] = sort(reqStations.y);
        traceOrder = sortYIndices;
        sortedValues = sortedYs;
        titlename = 'distance along strike';
        ylabelname = [titlename, ' (km), station'];
    case {'DISTANCEFROMEVENT','DISTFROMEVENT','DFE','OFFSET'}
        %calculate distance to event
        p1 =[reqStations.x; reqStations.y; reqStations.z];
        p2 = [eventC.cat.table.x(j); eventC.cat.table.x(j);...
            eventC.cat.table.x(j);];
        offsets = (p1 - p2).^2;
        
        %sort by distance from event
        [sortedOffsets, sortOffsetIndices] = sort(offsets);
        traceOrder = sortOffsetIndices;
        sortedValues = sortedOffsets;
        titlename = 'Offset';
        ylabelname = [titlename, ' (km), station'];
end

if ~any(contains(p.UsingDefaults,'scale'))
    scale = p.Results.scale;
    if strcmpi(method,'WIGBYY')
        markScale = abs(sortedValues(end) - sortedValues(1))/numel(traceOrder)*1;
    else
        markScale = 1;
    end
end



switch upper(plottype)
    case {'WIG', 'SHA', 'CWIG', 'BWIG', 'RAW'}
        % Find a suitable value for the scale of wiggles if none provided:
        if ~exist('scale', 'var')
            scale = 1;
            markScale = 1;
        end
        plot(eventC.(comp){j}, plottype, scale, traceOrder);
        
        ylabel(ylabelname,'fontsize',12)
        labelValues = sortedValues(get(gca,'YTick'))';
        if size(labelValues,1) < size(labelValues,2)
            labelValues = labelValues';
        end
        
        % If this is a stacked object, then display number of traces in
        % each stack and standard deviations
        nlabels = numel(labelValues);
        if isStacked
            newLabels = [char(stationNames),...
                num2str(labelValues, '%4.0f'),...
                repmat('\pm: ', nlabels, 1),...
                num2str(stdValues, '%3.1f'),...
                repmat(' n: ',nlabels,1),...
                num2str(eventC.cat.table.nStackedTraces(traceOrder),...
                '%2.0f')];
        else
            newLabels = [num2str(labelValues, '%4.0f '),...
                char(pad(stationNames(sortXIndices),6,'left'))];
        end
        
        set(gca,'YTickLabel', newLabels);
    case 'WIGBYY'
        if ~exist('scale', 'var')
           scale = abs(sortedValues(end) - sortedValues(1))/numel(traceOrder)*1;
           markScale = scale;
        end
        if scale == 0
            scale = abs(sortedValues(end) - sortedValues(1))/numel(traceOrder)*1;
            markScale = scale;
        end
            
        % Set y-axis to the sorted parameter, i.e. plot wiggle against that
        % sort-parameter
        plot(eventC.(comp){j}, plottype, scale, traceOrder, sortedValues,...
           'linewidth', linewidth);
        
        addsome = abs(sortedValues(end) - sortedValues(1))*0.02;      
        ylabel(ylabelname,'fontsize',12)
        
        ylim([min(sortedValues)-addsome, max(sortedValues)+addsome]);
        set(gca,'YTickMode','auto')
        set(gca,'YTickLabelMode','auto')
        set(gca,'YTickLabelMode','auto')
        set(gca,'YMinorTick','on')
end

%logical array to store whether arrivals were labelled on the last few
%lines or no
nLines = numel(traceOrder);
minLinesWithoutLabel = round(nLines/1);
lastXLines = false(minLinesWithoutLabel,1);

ax1 = gca;
hold(ax1, 'on')

if markArrivalTimes
    % check if arrivals were sent to function, else load from file
    if isempty(fields(arrivals))
        load('arrivals_0.02degGrid.mat','arrivals');
    end

    
    % chose only the arrivals for this event
    evID = eventC.cat.table.EventID(j);
    arrivals = arrivals(evID);
    if arrivals.EventID ~= evID
        warning('event IDs do not match')
        return
    end
    % plot arrivals 
    % check that theere are arrivals for that event
    if isempty(arrivals.arrivals)
        return
    end
    

    % initialize some helper variables
    prevEventArrivalPhases = cell(0,0);   
    pickOffset = zeros(nLines,1);

    thisEventArrivalPhases = cell(0,0);
    thisEventArrivalTimes = zeros(0,0);
    wavs = get(eventC.(comp){j},'waveform');
    triggers = get(eventC.(comp){j},'trig');
    
    for k = traceOrder'
        sta = get(wavs(k),'station');
        %Check that there
%         if traceOrder(k) > length(arrivals)
%             continue
%         end
                
        % get the right line for plotting
        %[arrivalLine, arrIdx] = find(traceOrder==k);
        arrivalLine = find(traceOrder==k);
        % us some absolute value if plotting by a value on the y-axis
        if strcmpi(plottype,'WIGBYY')
            arrivalLine = sortedValues(arrivalLine);
        end

        % check that there are arrivals for that station
        if ~isfield(arrivals.arrivals, sta)
            continue
        end
        if isempty(arrivals.arrivals.(sta).time)
            continue
        end
        narrivals = numel(arrivals.arrivals.(sta).time);

        % Two ways to find the first arrival: from the FM3D paths
        % or from the picks (triggers)
        firstArrival = getFirstArrivalTimeForEventFromFM3D(eventC,...
            arrivals, sta, traceOrder, arrivalLine);
        firstArrivalFromPicks = (triggers(k) - eventC.cat.table.otime(j))...
            * SECPERDAY;
        pickOffset(k) = firstArrival - firstArrivalFromPicks;

        q = 0;
        for p=1:1:narrivals
            % Check the number of conversions along a path, and only allow 
            % plotting when it is less than maxConversions
            thisPhase = char(arrivals.arrivals.(sta).phase(p));
            plotPhase = true;

            % do not plot when phase contains too many conversions
            if arrivals.arrivals.(sta).nConv(p) > maxConversions
                plotPhase = false;
            % don't plot if too many reflections
            elseif arrivals.arrivals.(sta).nReflections(p) >...
                    maxReflections
                plotPhase = false;
            % don't plot if it is a P-to -S conversion (except for events
            % from the lower plane of the DSZ)
            elseif strcmp(thisPhase(end),'S') && arrivals.arrivals.(sta...
                    ).nConv(p)>=1 %&& eventC.cat.table.DistFromSlabTop > -15 
                plotPhase = false;
            % if polarizationfiltered: don't plot the arrivals that
            % shouldn't be visible on that channel
            % elseif strcmp(comp,'Zp') && strcmp(thisPhase(end),'S')
            %     plotPhase = false;
            % elseif (contains(comp,'R') || contains(comp,'T')) &&...
            %         strcmp(thisPhase(end),'P')
            %     plotPhase = false;
            elseif isThisPhaseTooLong(thisPhase, maxPhaseLength,...
                    eventC.cat.table.DistFromSlabTop(j))
                plotPhase = false;
            % don't plot if too many interactions
            elseif  arrivals.arrivals.(sta).nConv(p) +...
                    arrivals.arrivals.(sta).nReflections(p)...
                    > maxInteractions
                plotPhase = false;
                % still plot if reflection and conversion was the
                % same interaction
                % need to be no more than x interactions and x
                % reflections
                if arrivals.arrivals.(sta).nConv(p)==...
                        maxConversions && arrivals.arrivals.(...
                        sta).nReflections(p) == maxReflections
                    if are_reflec_and_convers_oneInteract(thisPhase)
                        plotPhase = true;
                    end
                end                        
            end

            if plotPhase
                q = q + 1;
                thisEventArrivalPhases{q,1} = char(arrivals.arrivals.(...
                    sta).phase(p));
                thisEventArrivalTimes(q,1) = arrivals.arrivals.(...
                    sta).time(p) - firstArrival;
            end


            if plotPhase
                x = arrivals.arrivals.(sta).time(p) - firstArrival;
                if x < 0
                    warning(['Arrival earlier than first arrival,',...
                        'something went wrong']);
                end
                y = arrivalLine;
                if strcmpi(plottype,'CWIG')
                    pickColor = [1 0.64 0];
                elseif strcmpi(plottype, 'WIG') ||...
                    strcmpi(plottype, 'WIGBYY') || strcmpi(plottype,'BWIG')
                    pickColor = [1 0 0];
                else
                    pickColor = [0 0 0];
                end
                plot(ax1, [x, x], [y-markScale/2 y+markScale/2]+plusEnv,...
                    'color', pickColor, 'LineWidth', 2);
            end
        end
        
        if labelArrivals
            % Mark the phases with text only if the
            % previous event did not have the same arrivals
            plotText = true;
            if length(thisEventArrivalPhases) == length(...
                    prevEventArrivalPhases)
                phasesEqual = strcmp(prevEventArrivalPhases,...
                    thisEventArrivalPhases);
                if all(phasesEqual)
                    plotText = false;
                end
            end
            %mark arrivals only for the first MW, first IF, first crustal,
            %first slabmantle EQ
            %don't label when the last 10 lines weren't empty
            if any(lastXLines)
                plotText = false;
            end
                
            if plotText
                [sortedEventArrivalTimes, sortI] = sort(thisEventArrivalTimes);
                sortedEventPhases = thisEventArrivalPhases(sortI,:);
                for p=1:1:length(thisEventArrivalTimes)
                    %get a shortened phase string
                    shortPhase = shorten_phase_descriptor(...
                        subset(eventC.cat,j),...
                        sortedEventPhases{p}, sortedEventArrivalTimes(p));
                    %if there is long enough between arrivals, then print
                    % with defautl alignment
                    if strcmp(shortPhase,'')
                        continue
                    end
                    %y - scale
                    ylimits = ylim;
                    text(ax1, sortedEventArrivalTimes(p), ylimits(1),...
                        shortPhase, 'color', pickColor,'Rotation', 90,...
                        'FontSize', 9, 'HorizontalAlignment','left',...
                        'VerticalAlignment','middle');
                    %if mod(p,2)
                    %    text(ax1, sortedEventArrivalTimes(p), y + scale,...
                    %        shortPhase, 'color', pickColor,'Rotation', 90,...
                    %        'FontSize', 10, 'HorizontalAlignment','left');
                    %else
                    %    text(ax1, sortedEventArrivalTimes(p), y - scale,...
                    %        shortPhase, 'color', pickColor,'Rotation', 90,...
                    %        'FontSize', 10, 'HorizontalAlignment','right');
                    %end
                end
            end
            lastXLines(1:end-1,1) = lastXLines(2:end,1);
            lastXLines(end,1) = plotText;
        end        
        
        prevEventArrivalPhases = thisEventArrivalPhases;
    end
end

if strcmpi(plottype,'SHA') && scale~=1
    caxis(caxis./scale)
end


ax1 = gca;
eventTime = [strrep(eventC.cat.table.yyyy_mm_dd(j,:),'_','-'),'T',...
    eventC.cat.table.hh_mm_ss(j,:)];
title({[eventTime, ', channel: ', comp];['sorted by ', titlename, ]},...
    'fontsize',14)
%ax1.Title.Position(2) = ax1.Title.Position(2) * 1.05;
ybounds = get(gca,'YLim');
ax1.Title.Position(2) = ybounds(1) - diff(ybounds)*0.047;


nYTicks = length(ax1.YTickLabel);
maxNYTicks = 50;
if nYTicks > maxNYTicks
    yTickStep = round(nYTicks/maxNYTicks);
    newYTickLabels = ax1.YTickLabel(1:yTickStep:end,:);
    ax1.YTick = ax1.YTick(1:yTickStep:end);
    ax1.YTickLabel = newYTickLabels;
end
    
ax1.Position = [0.1 0.03 0.89 0.9];
ax1.XMinorTick = 'on';
xlim([-4 30])

fig1 = gcf;
oldpos = fig1.Position;
fig1.Position(4) = oldpos(4);
