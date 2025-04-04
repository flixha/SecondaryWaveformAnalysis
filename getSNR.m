function SNR = getSNR(c, cstation, comp)

% returns Signal-to-noise-ratio of all records of a station

comp = fieldnames(c);
cstation = fieldnames(c.Z)';

SNR = c;

for k=1:1:length(cstation)  
    if isfield(c.(comp{1}),cstation{k})
        XXXXXX need to finish tthis function
        SNR.(comp{1}).(cstation{k}).corr
        
        noiseWZ = get(crop(c.(comp{1}).(cstation{k}).corr,-2.5,-0.5),...
            'waveform');
        signalWZ = get(crop(c.(comp{1}).(cstation{k}).corr,0,2),...
            'waveform');
        SNR_Z = mean(abs(signalWZ))./mean(abs(noiseWZ));
        
        noiseWN = get(crop(c.(comp{2}).(cstation{k}).corr,-2.5,-0.5),...
            'waveform');
        signalWN = get(crop(c.(comp{2}).(cstation{k}).corr,0,2),...
            'waveform');
        SNR_N = mean(abs(signalWN))./mean(abs(noiseWN));
        
        noiseWE = get(crop(c.(comp{3}).(cstation{k}).corr,-2.5,-0.5),...
            'waveform');
        signalWE = get(crop(c.(comp{3}).(cstation{k}).corr,0,2),...
            'waveform');
        SNR_E = mean(abs(signalWE))./mean(abs(noiseWE));
        
    end
end