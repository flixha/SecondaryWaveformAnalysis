function c = subsetStationGather(c, EventIDs)


if size(EventIDs,2) > size(EventIDs,1)
    EventIDs = EventIDs';
end

comp = fieldnames(c);
sta = fieldnames(c.(comp{1}));

for nc=1:1:length(comp)
    for s=1:1:length(sta)
        
        if isa(c.(comp{nc}).(sta{s}),'threecomp')
            trIdx = ismember(c.(comp{1}).(sta{s}).cat.table.EventID, EventIDs);
            c.(comp{nc}).(sta{s}) = c.(comp{nc}).(sta{s})(trIdx);
        else
            trIdx = ismember(c.(comp{nc}).(sta{s}).cat.table.EventID, EventIDs);
            c.(comp{nc}).(sta{s}).corr = subset(c.(comp{nc}).(sta{s}...
                ).corr, trIdx);
            c.(comp{nc}).(sta{s}).cat = subset(c.(comp{nc}).(sta{s}...
                ).cat, trIdx);
        end
    end
end

