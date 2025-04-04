function eventCOut = resampleEventCorrObject(eventC, targetSamplingRate)
    %resample all waveforms in a network-correlation object to a selected
    %target sampling rate
    comp = fieldnames(eventC)';
    
    eventCOut = eventC;
    
    for k=2:1:numel(comp)
        for j=1:1:numel(eventC.(comp{k}))
            %check that we're not resampling the three component objects
            if ~isa(eventC.(comp{k}){j},'correlation')
                continue
            end
            
            wavs = get(eventC.(comp{k}){j}, 'waveform');

            
            %resample waveform to target Sampling rate
            crunchFactor = get(wavs,'freq') ./ targetSamplingRate;
            crunchFactor = round(crunchFactor, 4);
            [Q, P] = numden(sym(crunchFactor));
            Q = double(Q);
            P = double(P);
            if all(Q==Q(1)) && all(P==P(1))
                Q = Q(1); P = P(1);
                D = double(wavs);
                ResampleD = resample(D,P,Q);
            end

            parfor w=1:1:numel(wavs)
                wav = wavs(w);
                % put back into waveform, but don't forget to 
                % update the frequency               
                wav = set(wav,'data',ResampleD(:,w), 'Freq',...
                    targetSamplingRate);
                wavs(w) = wav;
            end
            eventCOut.(comp{k}){j} = set(eventCOut.(comp{k}){j},...
                'waveform', wavs);
        end
    end