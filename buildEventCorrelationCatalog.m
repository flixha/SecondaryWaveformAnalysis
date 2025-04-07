function [eventCorrelationCatalog] = buildEventCorrelationCatalog(seisCatalog,...
    eventSelection,correlationStations,correlationChannels,targetSamplingRate)

    %function to build a complete correlatable catalog for many stations, many
    %channels, and many events. those can then be individually accessed within
    %that structure through a "." (dot)-notation.
    %seisanCatalog is a catalog object of earthquakes with waveforms 
    % correlationStations is a cell array with the names of the desired
    % stations
    % correlationChannels is a cell array with the names of the desired
    % channels
    % eventSelection are the indices of the events to be processed
%%

    SECPERDAY = 60*60*24;

    cstation=correlationStations;
    corcomp=correlationChannels;
    eventI = eventSelection;
    if size(eventI,1) > size(eventI,2)
        eventI = eventI';
    end
    
    eventC = struct();
    eventC.cat = subset(seisCatalog, eventI);
    
    for ncc=1:1:numel(corcomp)
        eventC.(corcomp{ncc}) = cell(0,0);
    end
    

    
    %for each event
    m = 0;
    for j=[eventI]
        m = m+1;
        w = seisCatalog.waveforms{j};
        if isempty(w)
            warning(['No waveforms for event ', num2str(j), ' ',...
                datestr(seisCatalog.otime(j),'yyyy-mm-dd HH:MM:SS')])
            continue
        end
            
        wstation = get(w,'station');
        wchannel = get(w,'channel');
        
        % I don't understand why this happens, but sometimes the above two
        % get statements return a X*X cell of station names for a
        % waveform object with only X waveforms
        if size(wstation,1) == size(wstation,2) && size(wstation,1) > 1
            wstation = wstation(:,1);
        end
        if size(wchannel,1) == size(wchannel,2) && size(wchannel,1) > 1
            wchannel = wchannel(:,1);
        end
        wchannel = char(wchannel);
        wcomponent = cellstr(wchannel(:,end));
        %for each channel
        for ncc=1:1:numel(corcomp)
            eventC.(corcomp{ncc}){m,1} = correlation();
            %for each station
            for k=1:1:length(cstation)
                staCompWavI = find(strcmp(wstation,cstation{k}) &...
                                   strcmp(wcomponent, getCorCompValue(corcomp{ncc})));
                
                if ~isempty(staCompWavI)
                    clear prevWav prevTrig ptriggers;
                    prevWav = get(eventC.(corcomp{ncc}){m,1},'waveform');
                    prevTrig = get(eventC.(corcomp{ncc}){m,1},'trig');
                
                    %find the trigger
                    trigSta = {seisCatalog.arrivals{j}.stacode};
                    trigStaI = find(strcmp(trigSta,cstation{k}));
                
                    %only add waveform and trigger to the structure if
                    %there is a pick
                    if ~isempty(trigStaI)
                        
                        ptriggers={seisCatalog.arrivals{j}.p_time};
                        striggers={seisCatalog.arrivals{j}.s_time};
                        
                        if ~isempty(ptriggers{trigStaI})

                            testWav = seisCatalog.waveforms{j}(staCompWavI);
                            testWav = testWav(1);
                            testfreq = get(testWav,'freq');
                            testdata = get(testWav,'data');
                            if testfreq > 9.99 & sum(testdata)~=0
                                %resample waveform to target Sampling rate
                                crunchFactor = get(testWav,'freq') /...
                                    targetSamplingRate;
                                crunchFactor = round(crunchFactor,4);
                                [Q, P] = numden(sym(crunchFactor));
                                Q = double(Q);
                                P = double(P);

                                % assume W is an existing waveform
                                D = double(testWav);
                                ResampleD = resample(D,P,Q); 
                                % put back into waveform, but don't forget 
                                % to update the frequency
                                testWav = set(testWav,'data',ResampleD,...
                                    'Freq', targetSamplingRate);
                                % and for good measure... update the 
                                % waveform's history
                                % testWav = testWav.addHistory(...
                                % 'Resampled data outside of waveform object');

                                %set the waveform for the correlation event
                                setWav = testWav;
                                if exist('prevWav')
                                    setWav=[prevWav;setWav];
                                end
                                %setWav = extract(setWav, 'TIME',...
                                %    seisCatalog.otime(j),...
                                %    seisCatalog.otime(j)+90/SECPERDAY);
                                eventC.(corcomp{ncc}){m,1} = set( ...
                                    eventC.(corcomp{ncc}){m,1},'waveform',setWav);


                                % set the trigger to either the pick or an
                                % assumed value
                                if ~isempty(trigStaI) && ~isempty(...
                                        ptriggers{trigStaI})
                                    setPTrig=ptriggers{trigStaI};
                                else
                                    setPTrig=[];
                    %                 setPTrig=NaN;
    %                                 setPTrig=greekevents.otime(j)+10/(60*60*24);
                                end

                                if ~isempty(trigStaI) && ~isempty(...
                                        ptriggers{trigStaI})
                                    setSTrig=striggers{trigStaI};
                                else
                                    setSTrig=NaN;
                                end

                                if exist('prevTrig')
                                    setPTrig=[prevTrig;setPTrig];
                                end
                                eventC.(corcomp{ncc}){m,1} = set(eventC.(...
                                    corcomp{ncc}){m,1},'trig',setPTrig);

                                % sort the arrivals for that very station for this event to
                                % the event-at-station-catalog
                                % eventC.(corcomp{ncc}).(cstation{...
                                % k}).cat.table.arrivals{end,1} =...
                                % seisCatalog.table.arrivals{j}(trigStaI);
                            end
                        end
                    end
                end
            end
        end
    end
    
    eventCorrelationCatalog = eventC;
end
