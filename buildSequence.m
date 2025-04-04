function outSequence = buildSequence(inSequence,pairIDs,pairRow)
    %this is a recursive function that builds a repeater sequence based on
    % -a sequence (list of eventIDs within that sequence)
    % -a list of repeater pairs
    % -the current row of the pair-array of which the partner (2nd ID in
    % row) is checked/added to the sequence
    outSequence = inSequence;
    %if the sequence contains the EventID of the event's partner
    partner = pairIDs(pairRow,2);
    partnerCol = find(inSequence==partner);
    if ~isempty(partnerCol)
        return
    else
        %find the end of the sequence (the first NaN in the line)
        col = min(find(isnan(inSequence)));
        outSequence(col) = partner;
        partnersPartners=find(pairIDs(:,1)==partner);
        for p=1:1:length(partnersPartners)
            outSequence = buildSequence(outSequence,pairIDs,partnersPartners(p));
        end
    end
end