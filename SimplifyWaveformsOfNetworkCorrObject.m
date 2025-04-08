function cOut = SimplifyWaveformsOfNetworkCorrObject(c, reductionFactor)
    %resample all waveforms in a network-correlation object to a selected
    %target sampling rate
    comp = fieldnames(c)';
    cstation = fieldnames(c.(comp{1}))';
    
    cOut = c;
    
    for k=1:1:numel(comp)
        for s=1:1:numel(cstation)
            %check that we're not resampling the three component objects
            if isa(c.(comp{k}).(cstation{s}),'threecomp')
                continue
            end
            if ~isa(c.(comp{k}).(cstation{s}).corr,'correlation')
                continue
            end          
            wavs = get(c.(comp{k}).(cstation{s}).corr, 'waveform');
            
            % parfor w=1:1:numel(wavs)
            for w=1:1:numel(wavs)
                wav = wavs(w);
                
                data = get(wav,'data');
                ResampleD = data(1:reductionFactor:end);
                SamplingRate = get(wav,'freq') ./ reductionFactor;
                % put back into waveform, but don't forget to 
                % update the frequency               
                wav = set(wav,'data',ResampleD, 'Freq', SamplingRate);
                wavs(w) = wav;
            end
            cOut.(comp{k}).(cstation{s}).corr = set(...
                    cOut.(comp{k}).(cstation{s}).corr, 'waveform', wavs);
        end
    end