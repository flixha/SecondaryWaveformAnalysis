function eventCout = processEventCorrelObjects(eventC, minFreq, maxFreq,...
    minSNR)

    % This function processes correlation objects that each make up
    % seismograms of one event at different stations.

    channels = fields(eventC);
    channels = channels([2:end]);
    
    % first remove events with empty correlation objects
    acceptCorr = ones(eventC.cat.numberOfEvents, 1);
    for j=1:1:eventC.cat.numberOfEvents
        if isempty(get(eventC.Z{j},'waveform'))
            acceptCorr(j) = 0;
        end
    end
    if any(~acceptCorr)
        selectAccepted = find(acceptCorr);
        eventC.cat = subset(eventC.cat, selectAccepted);
        eventC.Z = eventC.Z(selectAccepted);
        eventC.N = eventC.N(selectAccepted);
        eventC.E = eventC.E(selectAccepted);
    end
    
    % now process
    for j=1:1:eventC.cat.numberOfEvents      
        for nc=1:1:numel(channels)
            eventC.(channels{nc}){j} = demean(eventC.(channels{nc}){j});
            eventC.(channels{nc}){j} = detrend(eventC.(channels{nc}){j});
            
             % Fill NaN-value gaps
            waves = get(eventC.(channels{nc}){j},'waveform');
            waves = fillgaps(waves,0);
            eventC.(channels{nc}){j} = set(eventC.(channels{nc}){j},...
                'waveform', waves);

            eventC.(channels{nc}){j} = crop(eventC.(channels{nc}){j},...
                -10,30);
            eventC.(channels{nc}){j} = taper(eventC.(channels{nc}){j});
            if minFreq > 0 && maxFreq < 100
                eventC.(channels{nc}){j} = butter(eventC.(channels{nc}){j},...
                    [minFreq maxFreq]);
            end
        end
        disp(['Event Object processing completed to ',...
            num2str(j/eventC.cat.numberOfEvents * 100, '%3.0f'), ' %'])
    end
    
    eventC = requestMinimumSNRofEventObj(eventC, minSNR);
    
    
    eventCout = eventC;