function cOut = resampleNetworkCorrObject(c, targetSamplingRate)
    %resample all waveforms in a network-correlation object to a selected
    %target sampling rate
    comp = fieldnames(c)';
    cstation = fieldnames(c.(comp{1}))';
    cstation = cstation(~strcmpi(cstation, 'cat'));
    
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

            
            %resample waveform to target Sampling rate
            crunchFactor = get(wavs,'freq') ./ targetSamplingRate;
            crunchFactor = round(crunchFactor, 4);
            if length(sym(crunchFactor)) > 20
                crunchFactor = round(crunchFactor, 2);
            end
            
            [Q, P] = numden(sym(crunchFactor));
            Q = double(Q);
            P = double(P);
            if all(Q==Q(1)) && all(P==P(1))
                Q = Q(1); P = P(1);
                D = double(wavs);
                
                % Q and P can't be arbritrarily long; check if they need to
                % be rounded / "truncated":
                max_fact = Q .* P;
                if any(max_fact > 2^31)
                    cor_fact = max(Q, P);
                    % Q = round(Q ./ 2^32);
                    % P = round(P ./ 2^32);
                    Q = round(Q ./ cor_fact * 2^15);
                    P = round(P ./ cor_fact * 2^15);
                elseif any(Q == 0) || any(P == 0)
                    error('Resampling-coefficient cannot be zero')
                end
                ResampleD = resample(D,P,Q);
            end

            % To run this in parallel:
            % parfor w=1:1:numel(wavs)
            for w=1:1:numel(wavs)
                wav = wavs(w);
                % put back into waveform, but don't forget to 
                % update the frequency               
                wav = set(wav,'data', ResampleD(:,w), 'Freq',...
                    targetSamplingRate);
                wavs(w) = wav;
            end
            cOut.(comp{k}).(cstation{s}).corr = set(...
                    cOut.(comp{k}).(cstation{s}).corr, 'waveform', wavs);
        end
    end