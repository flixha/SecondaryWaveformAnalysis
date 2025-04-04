function STFestimate = estimate_source_time_function(workcatalog, eventC,...
    index, minimumChannelsForSTFStack, minSNR, weightBeforeStacking,...
    doplot)

if size(index,1) > size(index,2)
    index = index';
end
p = 1;
STFestimate = waveform();

%make a copy 
eventC2 = eventC;

for j=[index]
   % Conditions when not to attempt deconvolution:
    % 1. when total number of channels too low
    if numel(get(eventC.Z{j,1},'waveforms')) < minimumChannelsForSTFStack
        continue
    end
    
    % define boundaries of the expected STF limits - dependent on the
    % magnitude of the event
    mymag = workcatalog.mag(j);
    M0 = 10^(1.5 * mymag + 9.1);
    
    %corner frequency estimate
    % Prieto et al. 2004
    x = 0;
    %﻿Hiramatsu et al. 2002, for different stress drops.
    % low stress drop:
    x = 2;
    fc_est = 10^(-1/3 * (log10(M0) +x) + 5.67);
    
    %
    
    %STF duration, relation from Beresnev (2001)
    td = 0.6/fc_est;
    
    STFcorrBoudaries = [-td/2 td];
    %STFcorrBoudaries = [-0.2 0.3];
    
    % some preprocessing
    eventC.Z{j,1} = demean(eventC.Z{j,1});
    eventC.Z{j,1} = detrend(eventC.Z{j,1});
    eventC.Z{j,1} = crop(eventC.Z{j,1},-6,30);
    eventC.Z{j,1} = detrend(eventC.Z{j,1});
    %eventC.Z{j,1} = demean(eventC.Z{j,1});
    eventC.Z{j,1} = taper(eventC.Z{j,1});
    % eventC.Z{j,1} = butter(eventC.Z{j,1},[0.1 24.9]);
    % eventC.Z{j,1} = butter(eventC.Z{j,1},[fc_est/10 2*fc_est]);

    %calculate trigger-adjustment on envelopes
    eventC2.Z{j,1} = hilbert(eventC.Z{j,1});
    %eventC2.Z{j,1} = butter(eventC2.Z{j,1},[0.1 4]);
    eventC2.Z{j,1} = xcorr(eventC2.Z{j,1}, 0.5*STFcorrBoudaries);
    
    if doplot
        plotSTF = crop( eventC2.Z{j,1}, STFcorrBoudaries*4);
        plot(plotSTF, 'wig', 2);
        ax1 = gca;
        figSTF = gcf;
        for k=1:1:numel(ax1.Children)
            set(ax1.Children(k),'color','r')
            set(ax1.Children(k),'LineWidth',0.5)
        end
        xlim( STFcorrBoudaries * 4);
        grid on
        set(figSTF,'Position',[50 44  464  1061]);
        %eventC2.Z{j,1} = adjusttrig(eventC2.Z{j,1},'MIN');
        eventC2.Z{j,1} = adjusttrig(eventC2.Z{j,1},'MMM', td/2);

        plotSTF = crop( eventC2.Z{j,1}, STFcorrBoudaries*4);
        plot(plotSTF, 'wig', 2);
        ax2 = gca;
        for k=1:1:numel(ax2.Children)
            %set(ax2.Children(k),'color','r')
            set(ax2.Children(k),'LineWidth',1)
        end
        copyobj(ax2.Children,ax1);
        close(ax2.Parent);
    end
    

    %adjust triggers on actual waveforms, and cross-corr again
    eventC.Z{j,1} = set(eventC.Z{j,1}, 'trig', get(eventC2.Z{j,1},'trig'));
    eventC.Z{j,1} = xcorr(eventC.Z{j,1}, STFcorrBoudaries);
    eventC.Z{j,1} = linkage(eventC.Z{j,1});
    %eventC.Z{j,1} = norm(eventC.Z{j,1});


    % %find me those stations that I have deemed suitable for stacking in general
    % getWave = get(eventC.Z{j,1},'waveform');
    % recStations = get(getWave,'station');
    % stackI = find(ismember(recStations, stfStackStations));
    % setWave = getWave(stackI);
    % eventC.Z{j,1} = set(eventC.Z{j,1},'waveform',setWave);
    % 
    % %select traces for stack based on SNR
    % 
    % eventC.Z{j,1} = crop(eventC.Z{j,1},-0.5,1.5);
    % eventC.Z{j,1} = xcorr(eventC.Z{j,1}, STFcorrBoudaries);
    % 
    % % plot(eventC.Z{j,1},'dend')
    % eventC.Z{j,1} = stack(eventC.Z{j,1},stackI);
    % eventC.Z{j,1} = norm(eventC.Z{j,1});
    % % eventC.Z{j,1} = adjusttrig(eventC.Z{j,1},'MIN',1);
    % plot(eventC.Z{j,1},'wig')

    %NOW:
    % 1. start at station S011, S010, S014 or S016
    % eventC2.Z{j,1} = crop(eventC.Z{j,1},-0.2,0.6);
    eventC2.Z{j,1} = crop(eventC.Z{j,1}, -6.0, 10);
    eventC2.Z{j,1} = detrend(eventC2.Z{j,1});
    tempWaveforms = get(eventC2.Z{j,1},'waveform');
    tempTriggers = get(eventC2.Z{j,1},'trig');
    channelTags = get(tempWaveforms,'ChannelTag');
        
    % calculate the maximum relative amplitude of the first motion vs. 
    % double/XX the whole STF
    %  This is to chose a station as standard for comparison where there is
    %  a very clear first motion waveform compared to the following ringing
    cForSTFAmplitudeComp = crop(eventC2.Z{j,1}, 2*STFcorrBoudaries);
    cForSTFAmplitudeComp = detrend(cForSTFAmplitudeComp);
    cForSTFAmplitudeComp = norm(cForSTFAmplitudeComp);
    cForSTFAmplitudeComp = crop(cForSTFAmplitudeComp, 0.5*STFcorrBoudaries);
    waveformsForSTFAmplitudeComp = get(cForSTFAmplitudeComp,'waveforms');
    maxSTFAmp = max(abs(waveformsForSTFAmplitudeComp));

    
    %tempWaveforms_short = get(crop(eventC2.Z{j,1},...
    %    0.5*STFcorrBoudaries),'waveforms');
    tempWaveforms_short = get(crop(butter(eventC.Z{j,1},...
        [fc_est/10 2*fc_est]), 0.5*STFcorrBoudaries),'waveforms');

    %find the channel for a preferred standard waveform
    % XXX but check which of these has a very clear peaky first polarity -
    % problems arise for stations that record waveforms from closeby the
    % nodal planes
    k=0;
    stationPriorities = {'S011','S010','S014','S016','S013','PE02','PE05'};
    %stationPriorities = {'PE05','S011','S010','S014','S016','S013','PE02',};
    m = [];
    %while isempty(m)
    %   k = k+1;
    m = zeros(numel(stationPriorities),0);
    mSTFAmp = zeros(numel(stationPriorities),0);
    for k=1:1:numel(stationPriorities)
        f = find(strcmp(get(channelTags,'station'),stationPriorities{k}) & ...
            (strcmp(get(channelTags,'channel'),'BH0Z') |...
            strcmp(get(channelTags,'channel'),'BHZ') |...
            strcmp(get(channelTags,'channel'),'HHZ')) );
        if f~=0
            mSTFAmp(k) = maxSTFAmp(f);
        else
            mSTFAmp(k) = 0;
        end
    end
    %get the priority-channel with highest maxSTFAmp
    [bestSTFAmp, bestSTFAmpI] = max(mSTFAmp);
    
    standardWaveform_short = tempWaveforms_short(bestSTFAmpI);
    standardWaveform = tempWaveforms(bestSTFAmpI);
    standardTrigger = tempTriggers(bestSTFAmpI);
    % 2. Check that SNR (pre vs post-trigger) is good enough.

    preSignal = crop(eventC2.Z{j,1}, [-2.5 -0.5]);
    preSignalWav = get(preSignal,'waveforms');
    preSignalMeanAmp = mean(abs(preSignalWav));
    postSignal = crop(eventC2.Z{j,1}, [0 2.0]);
    postSignalWav = get(postSignal,'waveforms');
    postSignalMeanAmp = mean(abs(postSignalWav));
    mySNR = postSignalMeanAmp./preSignalMeanAmp;
    
    badSNRindex = find(mySNR < minSNR);
    badSNRchannelTags = get(tempWaveforms(badSNRindex),'channeltag');
    %set to empty if none with bad SNR is found so later comparisons work
    if isempty(badSNRchannelTags)
        badSNRchannelTags = ChannelTag();
    end
    % cancel deconvolution when number of SNR-acceptable channels too low
    if numel(tempWaveforms) - numel(badSNRchannelTags) < 10
        continue
    end

    % 3. compare the source waveform. If the correlation is positive, put into
    wavData1 = get(standardWaveform_short,'data');
    stackBin1 = correlation();
    stackBin1 = set(stackBin1,'waveforms',standardWaveform);
    stackBin1 = set(stackBin1,'trig',standardTrigger);

    stackBin2 = correlation();
    for k=1:1:numel(channelTags)
        compareWaveform = tempWaveforms_short(k);
        if get(standardWaveform,'channeltag') ~= get(compareWaveform,'channeltag')
            wavData2 = get(compareWaveform,'data');
            
            zeroLagCorr = xcorr(wavData1, wavData2, 0);
            % Put in same polarity-bin if corr positive. If correlation is negative, other bin.
            if zeroLagCorr >= 0 
                prevWaveforms1 = get(stackBin1,'waveforms');
                prevTriggers1 = get(stackBin1,'trig');
                stackBin1 = set(stackBin1,'waveforms',[prevWaveforms1;tempWaveforms(k)]);
                stackBin1 = set(stackBin1,'trig',[prevTriggers1;tempTriggers(k)]);
            else
                prevWaveforms2 = get(stackBin2,'waveforms');
                prevTriggers2 = get(stackBin2,'trig');
                stackBin2 = set(stackBin2,'waveforms',[prevWaveforms2;tempWaveforms(k)]);
                stackBin2 = set(stackBin2,'trig',[prevTriggers2;tempTriggers(k)]);
            end
        end
    end
    
    % link = linkage(eventDistance,'average');
    % clust = cluster(link,'maxclust',nclust);

    stackBin1 = crop(stackBin1, STFcorrBoudaries * 4);
    stackBin2 = crop(stackBin2, STFcorrBoudaries * 4);
    stackBin1 = detrend(stackBin1);
    stackBin2 = detrend(stackBin2);
    stackBin1 = taper(stackBin1);
    stackBin2 = taper(stackBin2);
    
    stackBin1 = xcorr(stackBin1, 0.5 * STFcorrBoudaries);
    stackBin2 = xcorr(stackBin2, 0.5 * STFcorrBoudaries);
    stackBin1 = adjusttrig(stackBin1,'MMM', td/2);
    stackBin2 = adjusttrig(stackBin2,'MMM', td/2);
    
    stackBin1 = xcorr(stackBin1, STFcorrBoudaries);
    stackBin2 = xcorr(stackBin2, STFcorrBoudaries);
    stackBin1 = norm(stackBin1);
    stackBin2 = norm(stackBin2);
    
    if weightBeforeStacking
        % Put different weights on waveforms in stack
        % 1. weight by average CC with all the other waveforms
        ccmatrix = get(stackBin1,'corr');
        wavs = get(stackBin1,'waveforms');
        meanCC = zeros(numel(wavs),1);
        for k=1:1:numel(wavs)
            meanCC(k) = mean(ccmatrix([1:k-1,k:end],k));
        end
        wavs = times(wavs, meanCC);
        stackBin1 = set(stackBin1, 'waveforms', wavs);
        %
        ccmatrix = get(stackBin2,'corr');
        wavs = get(stackBin2,'waveforms');
        meanCC = zeros(numel(wavs),1);
        for k=1:1:numel(wavs)
            meanCC(k) = mean(ccmatrix([1:k-1,k:end],k));
        end
        wavs = times(wavs, meanCC);
        stackBin2 = set(stackBin2, 'waveforms', wavs);

        % 2. reduce to 1 % if below SNR limit
        wavs = get(stackBin1,'waveforms');
        for k=1:1:numel(wavs)
            wav = wavs(k);
            wavchanneltag = get(wav,'channeltag');
            if any(wavchanneltag == badSNRchannelTags)
                wavs(k) = times(wav, 0.01);
            end
        end
        stackBin1 = set(stackBin1,'waveforms', wavs);
        %
        wavs = get(stackBin2,'waveforms');
        for k=1:1:numel(wavs)
            wav = wavs(k);
            wavchanneltag = get(wav,'channeltag');
            if any(wavchanneltag == badSNRchannelTags)
                wavs(k) = times(wav, 0.01);
            end
        end
        stackBin2 = set(stackBin2,'waveforms', wavs);

    %     for n=1:1:numel(wavs)
    %         if any(isnan(get(wavs(n),'data')))
    %             wavs(n)
    %         end
    %     end
    end
    
    % and stack each bin!    
    stackBin1 = stack(stackBin1);
    stackBin2 = stack(stackBin2);
    
    % and then stack the the two bins 
    %get the last, stacked trace of each stack
    STFstack1 = subset(stackBin1,numel(get(stackBin1,'waveforms')));
    STFstack2 = subset(stackBin2,numel(get(stackBin2,'waveforms')));
    % flip the second one around
    STFstack2 = set(STFstack2,'waveforms',times(get(STFstack2,'waveforms'),-1));
    STFstack = correlation();
    STFstack = set(STFstack,'waveforms',[get(STFstack1,'waveforms'); get(STFstack2,'waveforms')]);
    STFstack = set(STFstack,'trig',[get(STFstack1,'trig'); get(STFstack2,'trig')]);
    STFstack = xcorr(STFstack, STFcorrBoudaries);
    STFstack = adjusttrig(STFstack,'MMM', td/2);
    %STFstack = norm(STFstack);
    STFstack = stack(STFstack);
    
    STFstack = crop(STFstack, STFcorrBoudaries .* [0.25, 1.25]);
    STFstack = taper(STFstack);
    STFstack = xcorr(STFstack, STFcorrBoudaries);
    ccmtx = get(STFstack,'corr');

    % and plot
    if doplot
        fig = figure;
        plot(stackBin1,'wig')
        fig1 = gcf;
        title({['Event ', datestr(workcatalog.otime(j),'yyyy-mm-dd HH:MM:SS'),...
            ' Magnitude: ', num2str(workcatalog.mag(j))];...
            'Deconvolution stack bin 1'})
        %xlim([-1 2])
        xlim( 4 * STFcorrBoudaries)    
        ax1 = copyobj(gca,fig);
        ax1.Position = [0.2 0.65 0.75 0.3];
        close(fig1)

        plot(stackBin2,'wig')
        fig2 = gcf;
        title('Deconvolution stack bin 2')
        %xlim([-1 2])
        xlim( 4 * STFcorrBoudaries)
        ax2 = copyobj(gca,fig);
        ax2.Position = [0.2 0.28 0.75 0.3];
        close(fig2)

        plot(STFstack,'wig')
        title({'Sum of the stacks, with bin 2 flipped upside down';...
            ['Correlation ', num2str(ccmtx(1,2),2)]})
        fig3 = gcf;
        %xlim([-1 2])
        xlim( 4 * STFcorrBoudaries)
        ax3 = copyobj(gca,fig);
        close(fig3)
        ax3.Position = [0.2 0.05 0.75 0.15];
        set(fig, 'Position', [500 44 464 1061]);

    %     rep = waitforbuttonpress;
    %     close(fig)
        uiwait(fig); disp('figure closed')
        close(figSTF)
    end
    
    
    STFestimate(p,1) = get(subset(STFstack,...
        numel(get(STFstack,'waveforms'))),'waveform');
    p = p+1;
end