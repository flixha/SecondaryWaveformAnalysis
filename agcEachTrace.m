function c2 = agcEachTrace(c, winLength)

    c2 = c;
    
    comp = fieldnames(c2);
    cstation = fieldnames(c2.(comp{1}));
    
    for k=1:1:length(cstation)
        for nc = 1:1:length(comp)-1
            c2.(comp{nc}).(cstation{k}).corr = agc(c2.(comp{nc}).(...
                cstation{k}).corr,winLength);
        end
    end
    