function eventCout = requestMinimumSNRofEventObj(eventC, minimumSNR)

% takes a network-object and removes all 3-c traces where no component has
% a signal-to-noise ratio below minimumSNR. 

comp = fieldnames(eventC);
comp = comp([2:end]);

eventCout = eventC;

number_removed = 0;
number_kept = 0;

for j=1:1:eventC.cat.numberOfEvents
    
     
    if length(get(eventC.Z{j},'waveform')) > 0
        noiseWZ = get(crop(eventC.Z{j}, -2.5, -0.5), 'waveform');
        signalWZ = get(crop(eventC.Z{j}, 0, 2), 'waveform');
        
        noiseWN = get(crop(eventC.N{j}, -2.5, -0.5), 'waveform');
        signalWN = get(crop(eventC.N{j}, 0, 2), 'waveform');     

        noiseWE = get(crop(eventC.E{j},-2.5,-0.5), 'waveform');
        signalWE = get(crop(eventC.E{j},0,2), 'waveform');      
        
        % select only those events where all channels are available
        Zstations = get(signalWZ,'station');
        Nstations = get(signalWN,'station');
        Estations = get(signalWE,'station');
        threeCompStations = Zstations(contains(Zstations, Nstations));
        threeCompStations = threeCompStations(contains(...
            threeCompStations, Estations));
        if numel(threeCompStations) ~= mean(...
                [numel(Zstations);numel(Nstations);numel(Estations)])
            % then we need to select only those where there are three!
            selectZsta = find(contains(Zstations, threeCompStations));
            selectNsta = find(contains(Nstations, threeCompStations));
            selectEsta = find(contains(Estations, threeCompStations));
            eventC.Z{j} = subset(eventC.Z{j}, selectZsta);
            eventC.N{j} = subset(eventC.N{j}, selectNsta);
            eventC.E{j} = subset(eventC.E{j}, selectEsta);
            
            % and get noise and signal again
            noiseWZ = get(crop(eventC.Z{j}, -2.5, -0.5), 'waveform');
            signalWZ = get(crop(eventC.Z{j}, 0, 2), 'waveform');
            noiseWN = get(crop(eventC.N{j}, -2.5, -0.5), 'waveform');
            signalWN = get(crop(eventC.N{j}, 0, 2), 'waveform');     
            noiseWE = get(crop(eventC.E{j},-2.5,-0.5), 'waveform');
            signalWE = get(crop(eventC.E{j},0,2), 'waveform');
        end
        
        
        SNR_Z = mean(abs(signalWZ))./mean(abs(noiseWZ));
        SNR_N = mean(abs(signalWN))./mean(abs(noiseWN));
        SNR_E = mean(abs(signalWE))./mean(abs(noiseWE));
            
        goodSNR = find(SNR_Z > minimumSNR | SNR_N > minimumSNR |...
            SNR_E > minimumSNR);
    else
        goodSNR = [];
    end

    for nc=1:1:length(comp)
        eventCout.(comp{nc}){j} = subset(eventC.(comp{nc}){j}, goodSNR);
%         c2.(comp{nc}).(cstation{k}).cat = subset(...
%             c.(comp{nc}).(cstation{k}).cat, goodSNR);
    end

    number_removed = number_removed + length(SNR_Z) - length(goodSNR);
    number_kept = number_kept + length(goodSNR);

    disp([num2str(round(j/eventC.cat.numberOfEvents*100)), ' % done, kept: ',...
        num2str(number_kept),', removed: ',num2str(number_removed)])
end
