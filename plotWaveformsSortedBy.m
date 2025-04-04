function [pickOffset] = plotWaveformsSortedBy(c, comp, sta, varargin)

% varargin: method, plottype, linewidth, scale, 
% plottype is 'sha' or 'wig' 
SECPERDAY = 24*60*60;


p = inputParser;

defaultMethod = 'DISTANCEFROMTRENCH';
validMethods = {'LAT', 'LATITUDE', 'LON', 'LONGITUDE', 'X',...
    'DISTFROMTRENCH', 'DISTANCEFROMTRENCH', 'Y', 'DISTALONGSTRIKE',...
    'Z', 'DEPTH', 'DISTFROMSLABTOP', 'DISTANCEFROMSLABTOP', 'SLABTOP',...
    'DFST', 'S-P','SMP','SMINUSP'};
checkMethod = @(x) any(validatestring(x,validMethods));

defaultPlotType = 'WIG';
validPlotType = {'SHA', 'WIG', 'WIGBYY', 'CWIG', 'BWIG'};
checkPlotType = @(x) any(validatestring(x,validPlotType));

defaultLineWidth = 0.5;
defaultScale = 1;

defaultMarkArrivalTimes = false;

addRequired(p,'c',@isstruct);
addRequired(p,'comp',@ischar);
addRequired(p,'sta',@ischar);
addOptional(p,'method', defaultMethod, checkMethod);
addOptional(p,'plottype', defaultPlotType, checkPlotType);
addParameter(p,'linewidth', defaultLineWidth, @isnumeric)
addParameter(p,'scale', defaultScale, @isnumeric)
addParameter(p,'arrivals', @isstruct)
addParameter(p,'markArrivalTimes', defaultMarkArrivalTimes, @islogical)

defaultLabelArrivals = false;
addOptional(p,'labelArrivals', defaultLabelArrivals, @islogical);

defaultMaxConversions = 1;
addOptional(p,'maxConversions', defaultMaxConversions, @isnumeric);

defaultMaxReflections = 1;
addOptional(p,'maxReflections', defaultMaxReflections, @isnumeric);

defaultMaxInteractions = 1;
addOptional(p,'maxInteractions', defaultMaxInteractions, @isnumeric);

defaultMaxPhaseLength = 5;
addParameter(p,'maxPhaseLength', defaultMaxPhaseLength, @isnumeric)

parse(p, c, comp, sta, varargin{:})
method = p.Results.method;
plottype = p.Results.plottype;
linewidth = p.Results.linewidth;
markArrivalTimes = p.Results.markArrivalTimes;
maxConversions = p.Results.maxConversions;
labelArrivals = p.Results.labelArrivals;
maxReflections  = p.Results.maxReflections;
maxInteractions = p.Results.maxInteractions;
maxPhaseLength = p.Results.maxPhaseLength;



    

% if length(varargin) >= 1
%     method = varargin{1};
% else
%     method = 'DISTANCEFROMTRENCH';
% end
% 
% if length(varargin) >= 2
%     plottype = varargin{2};
% else
%     plottype = 'SHA';
% end
% 
% if length(varargin) >= 3
%     linewidth = varargin{3};
% else
%     linewidth = 0.5;
% end
% 
% if length(varargin) >= 4
%     scale = varargin{4};
%     % evaluate a suitable scale later if not given here
% end

invertplot = false;
isStackedCatalog = false;

if ismember('stddepth', c.(comp).(sta).cat.table.Properties.VariableNames)
    isStackedCatalog = true;
end

switch upper(method)
    case {'LAT','LATITUDE'}
        %sort by lat
        [sortedLats,sortLatIndices] = sort(c.(comp).(sta).cat.lat);
        traceOrder = sortLatIndices;
        sortedValues = sortedLats;
        titlename = 'latitude';
        ylabelname = titlename;
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stdlat(traceOrder);
        end
    case {'LON','LONGITUDE'}
        %sort by lat
        [sortedLons,sortLonIndices] = sort(c.(comp).(sta).cat.lon);
        traceOrder = sortLonIndices;
        sortedValues = sortedLons;
        titlename = 'longitude';
        ylabelname = titlename;
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stdlon(traceOrder);
        end        
    case {'X','DISTFROMTRENCH','DISTANCEFROMTRENCH'}
        %sort by distance from trench greekevents.table.x
        %- this sorting might be preferrable, conversion waveforms align quite well
        [sortedXs,sortXIndices] = sort(c.(comp).(sta).cat.table.x);
        traceOrder = sortXIndices;
        sortedValues = sortedXs;
        titlename = 'distance from trench';
        ylabelname = [titlename, ' (km)'];
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stdx(traceOrder);
        end
    case {'Y','DISTALONGSTRIKE'}
        %sort by distance along strike (along trench) greekevents.table.x
        [sortedYs,sortYIndices] = sort(c.(comp).(sta).cat.table.y);
        traceOrder = sortYIndices;
        sortedValues = sortedYs;
        titlename = 'distance along strike';
        ylabelname = [titlename, ' (km)'];
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stdy(traceOrder);
        end
    case {'Z','DEPTH'}
        %sort by depth greekevents.table.z
        [sortedDepths, sortDepthIndices] = sort(c.(comp).(sta).cat.table.depth);
        traceOrder = sortDepthIndices;
        sortedValues = sortedDepths;
        titlename = 'depth';
        ylabelname = [titlename, ' (km)'];
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stddepth(traceOrder);
        end
    case {'DISTFROMSLABTOP','DISTANCEFROMSLABTOP','SLABTOP','DFST'}
        %sort by distance from slab Top
        [sortedDFSTs,sortDFSTIndices] = sort(...
            c.(comp).(sta).cat.table.DistFromSlabTop,'ascend');
        traceOrder = sortDFSTIndices;
        sortedValues = sortedDFSTs;
        titlename = 'distance from slab top';
        ylabelname = [titlename, ' (km)'];
        invertplot = true;
        if isStackedCatalog
            stdValues = c.(comp).(sta).cat.table.stdDistFromSlabTop(traceOrder);
        end
    case {'S-P','SMP','SMINUSP'}
        %sort by S-P
        nsmp = size(c.(comp).(sta).cat.arrivals,1);
        SmPtime = zeros(nsmp,1);
        for j=1:1:nsmp
            %if there are any arrivals at all for that station and event
            if size((c.(comp).(sta).cat.arrivals{j}),2) > 0
                %if there are both p- and s- arrivals
                if ~isempty(c.(comp).(sta).cat.arrivals{j}.s_time) &&...
                        ~isempty(c.(comp).(sta).cat.arrivals{j}.p_time)
                    SmPtime(j) = c.(comp).(sta).cat.arrivals{j}.s_time-...
                        c.(comp).(sta).cat.arrivals{j}.p_time;
                end
            else
                SmPtime.(sta)(j) = NaN;
            end
        end
        SmPtime = SmPtime * SECPERDAY;
        [sortedSmP,sortSmPIndices] = sort(SmPtime);
        traceOrder = sortSmPIndices;
        sortedValues = sortedSmP;
        titlename = 'S-minus-P arrival time';
        ylabelname = [titlename, ' (s)'];
    case {'MAG','M','MAGNITUDE'}
        %sort by distance from slab Top
        [sortedMags,sortMagIndices] = sort(...
            c.(comp).(sta).cat.mag,'descend');
        traceOrder = sortMagIndices;
        sortedValues = sortedMags;
        titlename = 'Magnitude';
        ylabelname = titlename;
        %invertplot = true;
end

if ~any(contains(p.UsingDefaults,'scale'))
    scale = p.Results.scale;
    if strcmpi(method,'WIGBYY')
        markScale = abs(sortedValues(end) - sortedValues(1))/numel(traceOrder)*1;
    else
        markScale = 1;
    end
end

% plot(c.(comp).(sta).corr,'wig',1,traceOrder);
% plot(c.(comp).(sta).corr,'sha',1,traceOrder);
% plot(c.(comp).(sta).corr, plottype, 1, traceOrder);
% To plot .....
% plot(c.(comp).(sta).corr,'wig',0.05,traceOrder,c.(comp).(sta).cat.table.x,'linewidth',0.1);

% Switch around waveforms and plot when the chosen sort value should go
% from large at the top to small at the bottom (that applies only to the 
% distance from slab top for now)
if invertplot
    if strcmp(plottype,'wig') || strcmp(plottype,'wigbyy')
        w = get(c.(comp).(sta).corr,'waveform');
        w = times(w,-1);
        c.(comp).(sta).corr = set(c.(comp).(sta).corr,'waveform',w);
    end
    ydirection = 'normal';
else
    ydirection = 'reverse';
end
            
switch upper(plottype)
    case {'WIG', 'SHA', 'CWIG', 'BWIG'}
        % Find a suitable value for the scale of wiggles if none provided:
        if ~exist('scale', 'var')
            scale = 1;
            markScale = 1;
        end
        plot(c.(comp).(sta).corr, plottype, scale, traceOrder);
        
        ylabel(ylabelname,'fontsize',12)
        labelValues = sortedValues(get(gca,'YTick'))';
        if size(labelValues,1) < size(labelValues,2)
            labelValues = labelValues';
        end
        
        % If this is a stacked object, then display number of traces in
        % each stack and standard deviations
        nlabels = numel(labelValues);
        if isStackedCatalog
            newLabels = [ num2str(labelValues, '%6.2f'),...
                repmat('\pm: ', nlabels, 1),...
                num2str(stdValues, '%3.1f'),...
                repmat(' n: ',nlabels,1),...
                num2str(c.(comp).(sta).cat.table.nStackedTraces(traceOrder),...
                '%2.0f')];
        else
            newLabels = num2str(labelValues, '%6.2f');
        end
        
        set(gca,'YTickLabel', newLabels);
    case 'WIGBYY'
        if ~exist('scale', 'var')
           scale = abs(sortedValues(end) - sortedValues(1))/numel(traceOrder)*1;
           markScale = scale;
        end
        % Set y-axis to the sorted parameter, i.e. plot wiggle against that
        % sort-parameter
        plot(c.(comp).(sta).corr, 'wig', scale, traceOrder, sortedValues,...
           'linewidth', linewidth);
%         plot(c.(comp).(sta).corr, 'wig', scale, sortedValues,...
%             'linewidth', linewidth);
        
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
minLinesWithoutLabel = round(nLines/8);
lastXLines = false(minLinesWithoutLabel,1);

ax1 = gca;
hold(ax1, 'on')

if markArrivalTimes
%     load('arrivals.mat');
%     load('arrivals_VgridsDepthCorrectlyInverted.mat');
%     load('arrivals_wNconv.mat');
    load('arrivals_0.02degGrid.mat');
    
    prevEventArrivalPhases = cell(0,0);
    
    % Go from the bottom up if the plot is inverted at the end so that the
    % phase-text appears at the top
    if invertplot
        traceOrder = flip(traceOrder);
        sortedValues = flip(sortedValues);
    end
    
    pickOffset = zeros(nLines,1);
    
    for j=traceOrder'
        evID = c.(comp).(sta).cat.table.EventID(j);
        
        if arrivals(evID).EventID ~= evID
            warning('event IDs do not match')
            continue
        end
            
            
        thisEventArrivalPhases = cell(0,0);
        thisEventArrivalTimes = zeros(0,0);
        %Check that there
%         if traceOrder(j) > length(arrivals)
%             continue
%         end
                
        % get the right line for plotting
        [arrivalLine, arrIdx] = find(traceOrder==j);
        if invertplot
            arrivalLine = length(traceOrder) - arrivalLine + 1;
        end
        
        if strcmpi(method,'WIGBYY')
            arrivalLine = sortedValues(arrIdx);
        end

        % plot arrivals 
        % check that theere are arrivals for that event
        if ~isempty(arrivals(evID).arrivals)
            % check that there are arrivals for that station 
            if ~isempty(arrivals(evID).arrivals.(sta).time)
                narrivals = numel(arrivals(evID).arrivals.(sta).time);
                
                % Two ways to find the first arrival: from the FM3D paths
                % or from the picks (triggers)
                firstArrival = getFirstArrivalTimeFromFM3D(c, arrivals,...
                    sta, comp, traceOrder, evID, arrivalLine);

                triggers = get(c.(comp).(sta).corr,'trig');
                firstArrivalFromPicks = (triggers(j) -...
                    c.(comp).(sta).cat.table.otime(j)) * SECPERDAY;
                
                pickOffset(j) = firstArrival - firstArrivalFromPicks;

                q = 0;
                for k=1:1:narrivals
                    % Check the number of conversions along a path, and
                    % only allow plotting when it is less than
                    % maxConversions
                    thisPhase = char(...
                        arrivals(evID).arrivals.(sta).phase(k));
%                     plotPhase = true;
%                     nConv = 0;
%                     pIdx = 1:3:length(thisPhase);
%                     for p=2:1:length(pIdx)
%                         if ~strcmp(thisPhase(pIdx(p)),...
%                                 thisPhase(pIdx(p-1)))
%                             nConv = nConv + 1;
%                         end
%                     end
                    plotPhase = true;

                    % do not plot when phase contains too many conversions
                    if arrivals(evID).arrivals.(sta).nConv(k) > maxConversions
                        plotPhase = false;
                    % don't plot if too many reflections
                    elseif arrivals(evID).arrivals.(sta).nReflections(k) >...
                            maxReflections
                        plotPhase = false;
                    % don't plot if it is a P-to -S conversion
                    elseif strcmp(thisPhase(end),'S') &&...
                            arrivals(evID).arrivals.(sta).nConv(k)>=1
                        plotPhase = false;
                        
                    % if polarizationfiltered: don't plot the arrivals that
                    % shouldn't be visible on that channel
%                     elseif strcmp(comp,'Zp') && strcmp(thisPhase(end),'S')
%                         plotPhase = false;
%                     elseif (contains(comp,'R') || contains(comp,'T')) &&...
%                             strcmp(thisPhase(end),'P')
%                         plotPhase = false;
                    elseif isThisPhaseTooLong(thisPhase, maxPhaseLength,...
                            c.(comp).(sta).cat.table.DistFromSlabTop(j))
                        plotPhase = false;
                    % don't plot if too many interactions
                    elseif  arrivals(evID).arrivals.(sta).nConv(k) +...
                            arrivals(evID).arrivals.(sta).nReflections(k)...
                            > maxInteractions
                        plotPhase = false;
                        % still plot if reflection and conversion was the
                        % same interaction
                        % need to be no more than x interactions and x
                        % reflections
                        if arrivals(evID).arrivals.(sta).nConv(k)==...
                                maxConversions &&...
                                arrivals(evID).arrivals.(...
                                sta).nReflections(k) == maxReflections
                            if are_reflec_and_convers_oneInteract(...
                                thisPhase)
                                plotPhase = true;
                            end
                        end                        
                    end
                    
                    if plotPhase
                        q = q + 1;
                        thisEventArrivalPhases{q,1} = char(...
                            arrivals(evID).arrivals.(sta).phase(k));
                        thisEventArrivalTimes(q,1) = arrivals(...
                            evID).arrivals.(sta).time(k) - firstArrival;
                    end
                    

                    if plotPhase    
                        %plot(ax1, arrivals(evID).arrivals.(sta).time(k) - firstArrival,...
                        %    arrivalLine, '+k');
                        x = arrivals(evID).arrivals.(sta).time(k) - firstArrival;
                        if x < 0
                            warning(['Arrival earlier than first arrival,',...
                                'something went wrong']);
                        end
                        y = arrivalLine;
                        if strcmpi(plottype,'CWIG')
                            pickColor = [1 0.64 0];
                        elseif strcmpi(plottype, 'WIG') ||...
                            strcmpi(plottype, 'WIGBYY') ||...
                                strcmpi(plottype,'BWIG')
                            pickColor = [1 0 0];
                        else
                            pickColor = [0 0 0];
                        end

                        plot(ax1, [x, x], [y-markScale/2 y+markScale/2],...
                            'color', pickColor, 'LineWidth', 2);

                    end
                end
            end
        end
        
        if labelArrivals
            % Mark the phases with text only if the
            % previous event did not have the same arrivals
            plotText = true;
%                             thisEventArrivalPhases =...
%                                 arrivals(evID).arrivals.(sta).phase;
            if length(thisEventArrivalPhases) ==...
                    length(prevEventArrivalPhases)
                phasesEqual = strcmp(prevEventArrivalPhases,...
                    thisEventArrivalPhases);
                if all(phasesEqual)
                    plotText = false;
                end
            end
            %mark arrivals only for the first MW, first IF, first crustal,
            %first slabmantle EQ
            %if c.(comp).(sta).cat.table.
            %don't label when the last 10 lines weren't empty
            if any(lastXLines)
                plotText = false;
            end
                
            if plotText
                %Y = repmat(y, length(thisEventArrivalTimes),1);
                %text(ax1, thisEventArrivalTimes, Y,...
                %            thisEventArrivalPhases,...
                %            'color', pickColor,'Rotation', 90,...
                %            'FontSize', 10);
                [sortedEventArrivalTimes, sortI] = sort(thisEventArrivalTimes);
                sortedEventPhases = thisEventArrivalPhases(sortI,:);
                for k=1:1:length(thisEventArrivalTimes)
                    %if k>1
                    %    pickDist = sortedEventArrivalTimes(k) -...
                    %        sortedEventArrivalTimes(k-1);
                    %else
                    %    pickDist = 99;
                    %end
                    %get a shortened phase string
                    shortPhase = shorten_phase_descriptor(...
                        sortedEventPhases{k}, sortedEventArrivalTimes(k));
                    %if there is long enough between arrivals, then print
                    % with defautl alignment
                    if ~strcmp(shortPhase,'')
                        if mod(k,2)
                            text(ax1, sortedEventArrivalTimes(k),...
                                y + scale, shortPhase,...
                                'color', pickColor,'Rotation', 90,...
                                'FontSize', 10, 'HorizontalAlignment','left');
                        else
                            text(ax1, sortedEventArrivalTimes(k),...
                                y - scale, shortPhase,...
                                'color', pickColor,'Rotation', 90,...
                                'FontSize', 10, 'HorizontalAlignment','right');
                        end
                    end
                        
                        
                        
                        
                        
%                         if pickDist > 1
%                             text(ax1, sortedEventArrivalTimes(k), y,...
%                                     sortedEventPhases(k),...
%                                     'color', pickColor,'Rotation', 90,...
%                                     'FontSize', 10, 'HorizontalAlignment','right');
%                         else
%                             % print offset to the left/bottom
%                     
%                     if mod(k,2)
%                         text(ax1, sortedEventArrivalTimes(k), y,...
%                                 sortedEventPhases(k),...
%                                 'color', pickColor,'Rotation', 90,...
%                                 'FontSize', 10, 'HorizontalAlignment','left');
%                     else
%                         text(ax1, sortedEventArrivalTimes(k), y,...
%                             sortedEventPhases(k),...
%                             'color', pickColor,'Rotation', 90,...
%                             'FontSize', 10, 'HorizontalAlignment','right');
%                     end
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


set(gca, 'YDir', ydirection);
ax1 = gca;

title([sta, '.', comp, ' sorted by ', titlename],'fontsize',14)
ax1.Title.Position(2) = ax1.Title.Position(2) * 1.05;

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
xlim([-4 16])

fig1 = gcf;
oldpos = fig1.Position;
fig1.Position(4) = oldpos(4);
