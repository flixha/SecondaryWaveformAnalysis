function repeaters = findRepeaters(Cat, c, comp, csta, minCCC,...
    minPercHighCCC, minNHighCCC, maxDist)
% this function analyses a network structure and checks which events
% fulfill the criteria of repeaters (minimum correlation at a minimum
% number of minimum percentage of stations. Returns a repeater structure.
% The networkCorrelation object that this function requires already need
% correlation matrices in each correlation object.
%(Cat,c,comp,minCCC,minPercHighCCC,minNHighCCC,maxDist);

    
    %number of channels per station, stations and allround n of channels to
    %compare
    nChannels = length(comp);
    nStations = length(csta);
    nStationChannels = nChannels * nStations;

    corr_matrices=struct();
    %initialize empty 3-D matrices
    helpCorrMatrix = NaN(Cat.numberOfEvents,Cat.numberOfEvents,nStationChannels);
    corrNanOrZeroMatrices = NaN(Cat.numberOfEvents,Cat.numberOfEvents,nStationChannels);
    %fill 3-D matrices with correlation values and logic statement whether zero
    %or NaN, respectively
    for nc=1:1:length(comp) 
        for k=1:1:length(csta)
            if isfield(c.(comp{nc}),csta{k})
                % fill the global correlation matrices with NaNs, theninput the 
                % correlation coefficient at the right position in the global 
                % correlation matrices
                corr_matrices.(comp{nc}).(csta{k}) = NaN(Cat.numberOfEvents);
                corr_matrix = get(c.(comp{nc}).(csta{k}).corr,'CORR');
                eventIDs = c.(comp{nc}).(csta{k}).cat.table.EventID;
                corr_matrices.(comp{nc}).(csta{k})(eventIDs,eventIDs) = corr_matrix;

                %save all correlation matrices into one 3-dimensional array: then
                %they can be summed
                if isfield(corr_matrices.(comp{nc}),csta{k})
                    helpCorrMatrix(:,:, nChannels * k - (nChannels - nc)) =...
                        corr_matrices.(comp{nc}).(csta{k});
                    %save an equal-sized matrix with boolean whether CCC is NaN or zero
                    corrNanOrZeroMatrices(:,:, nChannels * k - (nChannels - nc)) =...
                        (isnan(corr_matrices.(comp{nc}).(csta{k})) | ...
                        corr_matrices.(comp{nc}).(csta{k})==0);
                else
                    %if that channel/station has no cross correlation matrix, say
                    %that everything is NaN
                    corrNanOrZeroMatrices(:,:, nChannels * k - (nChannels - nc)) =...
                         true(Cat.numberOfEvents,Cat.numberOfEvents);
                end
            end
        end
    end         

    % number of CCCs per event that are NaN or Zero (or not)
    nNanOrZero = sum(corrNanOrZeroMatrices,3);
    nNotNanOrZero = nStationChannels - nNanOrZero;
    % number of CCCs greather the set threshold
    % nCCCgreaterThreshold = sum( (helpCorrMatrix > minCCC ) == true,3);
    nCCCgreaterThreshold = sum( (helpCorrMatrix > minCCC & helpCorrMatrix ~=1) == true,3);

    %find events with high CCCs according to criteria
    isHighCorrelationEvent = (nCCCgreaterThreshold > nNotNanOrZero * minPercHighCCC |...
            nCCCgreaterThreshold >= minNHighCCC);
    % isHighCorrelationEvent = ( nCCCgreaterThreshold >= minNHighCCC);    
    %seperate similar events that are not themselves
    [similarEventRow, similarEventColumn] = find(isHighCorrelationEvent == true);
    %find the autocorrelations (event compared to itself at a station), and
    %take only events that are not themselves further
    similarEventsThatAreNotThemselves = find( similarEventRow ~= similarEventColumn);
    similarEventIDs = similarEventRow(similarEventsThatAreNotThemselves);
    similarPairIDs =  [similarEventRow(similarEventsThatAreNotThemselves), ...
        similarEventColumn(similarEventsThatAreNotThemselves)];

    %Check the distance between events with similar waveforms
    similarEventDistance = sqrt(...
        (Cat.table.x(similarPairIDs(:,1)) - Cat.table.x(similarPairIDs(:,2))).^2 +...
        (Cat.table.y(similarPairIDs(:,1)) - Cat.table.y(similarPairIDs(:,2))).^2 +...
        (Cat.table.z(similarPairIDs(:,1)) - Cat.table.z(similarPairIDs(:,2))).^2);

    %Define similar events as repeaters when they are close enough 
    repeatingEventRow = find(similarEventDistance < maxDist);
    repeaters.EventIDs = similarPairIDs(similarEventDistance < maxDist);
    repeaters.cat = Cat.table(repeaters.EventIDs,:);
    repeaters.PairIDs =  [similarPairIDs(repeatingEventRow,1), similarPairIDs(repeatingEventRow,2)];
    repeaters.PairDistance = similarEventDistance(repeatingEventRow);

    %% Now go through all detected repeater pairs and build a list of repeater sequences:
    % If event 1 highly correlates with event 2 and event 2 with event 3, then
    % they all make up one repeater sequence. Each sequence is stored in its
    % own row in repeaterSequences, and it is made up of EventIDs in the
    % row's columns.
    nRepeaters = length(unique(repeaters.EventIDs));
    repeaterSequences = NaN(nRepeaters);
    for j=1:1:length(repeaters.PairIDs)
        currentEventID = repeaters.PairIDs(j,1);
        [foundSeqRow, foundSeqCol] = find(repeaterSequences==currentEventID);
        if ~isempty(foundSeqRow)
            currentSequence = repeaterSequences(foundSeqRow,:);
            repeaterSequences(foundSeqRow,:) = buildSequence(currentSequence,repeaters.PairIDs,j);
        else
            putSequenceRow = min(find(isnan(repeaterSequences(:,1))));
            repeaterSequences(putSequenceRow,1) = currentEventID;
            currentSequence = repeaterSequences(putSequenceRow,:);
            repeaterSequences(putSequenceRow,:) = buildSequence(currentSequence,repeaters.PairIDs,j);
        end
    end
    repeaterSequences(:,~any(~isnan(repeaterSequences), 1))=[];
    repeaterSequences(~any(~isnan(repeaterSequences), 2),:)=[];
    for j=1:1:size(repeaterSequences,1)
        repeaterSequences(j,:) = sort(repeaterSequences(j,:));
    end
    repeaters.Sequences = repeaterSequences;

    sequNumMatrix = ones(size(repeaterSequences)).*[1:1:size(repeaterSequences,1)]';

    %now get indices of repeaters in the squence-order, needed for plotting
    sequenceSortedIDs = reshape(repeaterSequences',numel(repeaterSequences),1);
    sequNumbers = reshape(sequNumMatrix',numel(sequNumMatrix),1);
    sequNumbers(~any(~isnan(sequenceSortedIDs), 2),:)=[];
    sequenceSortedIDs(~any(~isnan(sequenceSortedIDs), 2),:)=[];
    repeaters.sequenceSortedIDs = sequenceSortedIDs;
    repeaters.sequenceNumbers = sequNumbers;

end