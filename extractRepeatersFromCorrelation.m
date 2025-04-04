function corrStruc = extractRepeatersFromCorrelation(c,repeaters,comp,cstation)
    %this function extract events with a specific index from a correlation
    %object, which allows me to easily extract the repeaters and plot them.
    %
%     EventIDs=unique(EventIDs);
%     EventIDs=sort(EventIDs);
    EventIDs = repeaters.sequenceSortedIDs;

    corrStruc = struct();
    for nc=1:1:length(comp)
        corrStruc.(comp{nc}) = struct();
        for k=1:1:length(cstation)
            corrStruc.(comp{nc}).(cstation{k}) = struct();
            % only if the networkCorr structure has this station and
            % channel
            if isfield(c.(comp{nc}),cstation{k})
                corrStruc.(comp{nc}).(cstation{k}).corr = correlation();
                corrStruc.(comp{nc}).(cstation{k}).cat = Catalog();
                %find all the waveforms in the correlation object for the
                %repeaters
                [corObjrow, repeaterIDrow] = find(c.(comp{nc}).(cstation{k}).cat.table.EventID == EventIDs');
                if ~isempty(corObjrow)
                    corrStruc.(comp{nc}).(cstation{k}).corr = subset(c.(comp{nc}).(cstation{k}).corr,corObjrow);
                    corrStruc.(comp{nc}).(cstation{k}).cat = subset(c.(comp{nc}).(cstation{k}).cat,corObjrow);
                    
%                     setWav = get(c.(comp{nc}).(cstation{k}).corr,'waveforms');
%                     setWav = setWav(corObjrow,:);
%                     corrStruc.(comp{nc}).(cstation{k}).corr =...
%                         set(corrStruc.(comp{nc}).(cstation{k}).corr,'waveforms',setWav);
% 
%                     setTrig = get(c.(comp{nc}).(cstation{k}).corr,'TRIG');
%                     setTrig = setTrig(corObjrow,:);
%                     corrStruc.(comp{nc}).(cstation{k}).corr =...
%                         set(corrStruc.(comp{nc}).(cstation{k}).corr,'TRIG',setTrig);
% 
%                     corrStruc.(comp{nc}).(cstation{k}).corr = demean(corrStruc.(comp{nc}).(cstation{k}).corr);
%                     corrStruc.(comp{nc}).(cstation{k}).corr = detrend(corrStruc.(comp{nc}).(cstation{k}).corr);
%                     corrStruc.(comp{nc}).(cstation{k}).corr = crop(corrStruc.(comp{nc}).(cstation{k}).corr,-4,14);
%                     corrStruc.(comp{nc}).(cstation{k}).corr = taper(corrStruc.(comp{nc}).(cstation{k}).corr);    
%                     corrStruc.(comp{nc}).(cstation{k}).corr = butter(corrStruc.(comp{nc}).(cstation{k}).corr,[2.5 8]);
%                     corrStruc.(comp{nc}).(cstation{k}).corr = xcorr(corrStruc.(comp{nc}).(cstation{k}).corr,[-0.1 12]);
%                     corrStruc.(comp{nc}).(cstation{k}).corr = adjusttrig(corrStruc.(comp{nc}).(cstation{k}).corr,'MIN');
%                     corrStruc.(comp{nc}).(cstation{k}).corr = sort(corrStruc.(comp{nc}).(cstation{k}).corr);
% 
                    %now set linkage and cluster IDs in the
                    %networkCorrelation structure
                    if length(get(corrStruc.(comp{nc}).(cstation{k}).corr,'trig'))>1
                        corrStruc.(comp{nc}).(cstation{k}).corr = linkage(corrStruc.(comp{nc}).(cstation{k}).corr);
                        clusterIDs = repeaters.sequenceNumbers(repeaterIDrow);
                        corrStruc.(comp{nc}).(cstation{k}).corr =...
                            set(corrStruc.(comp{nc}).(cstation{k}).corr,'clust', clusterIDs);
                    end
                end

            end
        end
    end

end