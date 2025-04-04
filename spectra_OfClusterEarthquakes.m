clear all
% close all

loadFromPreviousState = true;
loadCorrelationsFromPreviousState = false;
catalogFile = 'greekevents_629_2018-10-24.mat';
corObjFile = 'greekevents_CorrelationObject_483_2018-08-01.mat';

geocBin = ['~/Documents2/SIMULR_5Datasets_runs_nConv/',...
    'GeosphereCorrection/bin/LL2carteGEOC_v3'];
SEISAN_TOP='/Volumes/Stilton/SEI';


seisanREAdir{1} = fullfile(SEISAN_TOP, 'REA', 'MEDUT');
seisanREAdir{2} = fullfile(SEISAN_TOP, 'REA', 'MEDET');
seisanREAdir{3} = fullfile(SEISAN_TOP, 'REA', 'GRNW6');
seisanWAVdir = fullfile('/Volumes/nasdata/SEISAN/WAV');

SIMULTOP = '/Volumes/Stilton/Documents/Greece_MWcluster/Relocation3D/';
simulrRELOCdir{1} = fullfile(SIMULTOP, 'MEDUT_Reloc');
simulrRELOCdir{2} = fullfile(SIMULTOP, 'MEDET_Reloc');
simulrRELOCdir{3} = fullfile(SIMULTOP, 'GRNW_Reloc');

detTimeWindow=seconds(10);

loadDD=true;
DDfolder = ['/Volumes/Stilton/Documents/Greece_MWcluster/Relocation_hypoDD/',...
    '682Events/MinLink_8_6/DD3D_wCC_NoAbs_NoTaper_Corr0.7/']; 
nclusters=14; zCor=11; DDiteration=25;

nDataBases=length(seisanREAdir);


cstation = {'S009','S010','S011','S012','S013','S014','S015','S016',...
    'PE02','PE05','PE07'};


%% now load in all seisan s-files and waveforms and combine different 
% Seisan-databases to one catalog "greekevents"
if loadFromPreviousState
    load(catalogFile, 'greekevents');
else
    greekevents = readMEDUSAseisanDatabases(SEISAN_TOP, seisanREAdir,...
        loadWaveforms, seisanWAVdir, startTime, endTime, simulrRELOCdir,...
        loadDD, DDfolder, loadJackknifeErrorsFromFile, detTimeWindow,...
        requestStations);
end

%% find earthquakes in subclusters

% !!southern Tripoli cluster!!
% minY = -35;
% maxY = 10;
% as wide as possible, but still rather vertical incoming rays
% maxY = 50;
% minY = -50;
maxY = 25;
minY = -35;
minZ = 20;
maxZ = 200;

%get earthquakes within a region around the GRT cross section image
%interface fault cluster
I_IFCluster = find(greekevents.table.y > -15 & greekevents.table.y < -4.9 &...
    greekevents.table.x > 83.5 & greekevents.table.x < 93.5 &...
    greekevents.table.depth > 54 & greekevents.table.depth < 60 );
%within-MW-cluster
% I_MWCluster = find(greekevents.table.y > -35 & greekevents.table.y < 11 &...
%     greekevents.table.x > 69 & greekevents.table.x < 83 &...
%     greekevents.table.depth < 56);
I_MWCluster = find(greekevents.table.y > -20 & greekevents.table.y < 0 &...
    greekevents.table.x > 69 & greekevents.table.x < 83 &...
    greekevents.table.depth < 56);
%intraslabclsuter
I_ISCluster = find(greekevents.table.y > -10 & greekevents.table.y < 0 &...
    greekevents.table.x > 65 & greekevents.table.x < 113 );
% I_ISCluster=find(clusterEvents.table.x > 65 & clusterEvents.table.x < 115 );
I_ISCluster = I_ISCluster(~ismember(I_ISCluster,I_IFCluster) &...
    ~ismember(I_ISCluster,I_MWCluster));

eventI = [I_MWCluster; I_IFCluster; I_ISCluster]';

%other events below the array
I_ArrayCluster = find(...
    greekevents.table.x > 65 & greekevents.table.x < 125 &...
    greekevents.table.y > minY & greekevents.table.y < maxY &...
    greekevents.table.z > minZ & greekevents.table.z < maxZ);
I_ArrayCluster = I_ArrayCluster(~ismember(I_ArrayCluster,eventI));


%northern Kalavrita cluster
% lat=38.15;
% lon=22.03;

%Interface Fault Cluster
% eventI = I_MWCluster';
% eventI = I_IFCluster';
% eventI = I_ISCluster';

% eventI = [I_MWCluster; I_IFCluster; I_ISCluster]';
eventI = [I_MWCluster; I_IFCluster; I_ISCluster; I_ArrayCluster]';
% eventI = I_IFCluster;

%all:
% eventI=[1:1:greekevents.numberOfEvents];
eventI = sort(eventI);

% selected high magnitude/quality-events:
% eventI = [276, 313, 10, ...
%     4, 12, 144,...
%     232, 468, 387,...
%     337, 120 , 225];

%maybe MW: 514, 481, IS: 566, 538, deeper: 164, 582

nSelectedEvents = length(eventI);

%Interface Fault Cluster
% clusterI{1} = I_MWCluster';
% clusterI{2} = I_IFCluster';
% clusterI{3} = I_ISCluster';

% select different sets based on magnitude
binBound = [0 0.5 1 1.5 2 2.5];
% binBound = [0.5 1.5 2.5];
clusterI = cell(numel(binBound)+1,1);
for j=1:1:numel(binBound)+1
    if j==1
        foundThem = find(greekevents.mag < binBound(j));
    elseif j==numel(binBound)+1
        foundThem = find(greekevents.mag >= binBound(j-1));
    else
        foundThem = find(greekevents.mag >= binBound(j-1) &...
            greekevents.mag < binBound(j));
    end
    clusterI{j} = foundThem(ismember(foundThem, eventI));
end


%% MAKE some FIGUREs 
% Show 3D model:
% [model3DFig] = plotEarthquakeCatalog3D(greekevents,eventI);

% Make another figure on top of GRT images:
% plotGRTfigureWithEarthquakes(greekevents,-30,20,true)


%% now get me all waveforms for the events that I want to compare into one
% structure of correlatable events
%that's the structure c: c.(component).(stationname) is a correlation object
% if loadCorrelationsFromPreviousState
%     load(corObjFile, 'c');
%     cstation = fieldnames(c.Z)';
% else
%     % CHeck whether greekevents still contains waveforms - if there are
%     % none, reload the whole database from .mat-file
%     if all(cellfun(@isempty,greekevents.waveforms))
%         clear greekevents;
%         load(catalogFile, 'greekevents');
%     end
%     c = buildCorrelationCatalog(greekevents, eventI, cstation, corcomp,...
%         100);
% end

%% spectra
specFig = figure;
specAx = axes(specFig);
hold(specAx,'on')
cmap = flipud(cbrewer('div','PiYG',11));
darkgreen = cmap(1,:);
darkmagenta = cmap(end,:);
orange = [255/255 130/255 0];
cmap = [darkmagenta; orange; darkgreen];

% the colormap with 3 different colors is for the specific case of
% splitting the dataset by mantle wedge/interface/intraslab earthquakes
if numel(clusterI) > 3
    cmap = jet(numel(clusterI));
end

clear crossStationSignalAmp  crossStationNoiseAmp crossStationSNRAmp;
clear averageSignalSpectrumHandle averageNoiseSpectrumHandle
totalnRemoved = 0;
%now compare for every cluster
for ncl=1:1:numel(clusterI)
    ntraces = 0;
    eventI = clusterI{ncl};
    nSelectedEvents = length(eventI);

    
    corcomp={'Z','E','N'};
    % now get me all waveforms for the events that I want to compare
    c = buildCorrelationCatalog(greekevents, eventI, cstation, corcomp,...
        100);
    
    %define signal and noise
    ncc = 1; %Z-comp
    for k=1:1:length(cstation)
        c_signal.(corcomp{ncc}).(cstation{k}).corr = crop(c.(...
            corcomp{ncc}).(cstation{k}).corr, 0, 10);
        c_signal.(corcomp{ncc}).(cstation{k}).cat = c.(corcomp{ncc}).(...
            cstation{k}).cat;
        c_noise.(corcomp{ncc}).(cstation{k}).corr = crop(c.(...
            corcomp{ncc}).(cstation{k}).corr, -11, -1);
        c_noise.(corcomp{ncc}).(cstation{k}).cat = c.(corcomp{ncc}).(...
            cstation{k}).cat;
    end
    
    %clear crossStationSignalAmp  crossStationNoiseAmp crossStationSNRAmp;
%     k=1;
    for s=[1:1:numel(cstation)]
        signal = get(c_signal.(corcomp{ncc}).(cstation{s}).corr,'waveforms');
        signal = detrend(demean(signal));
        noise = get(c_noise.(corcomp{ncc}).(cstation{s}).corr,'waveforms');
        noise = detrend(demean(noise));
        
        % check for spikes
        keepWavIdx = true(size(noise));
        for k=1:1:numel(noise)
            if max(diff(noise(k))) > 10000
                keepWavIdx(k) = false;
            end
        end
        signal = signal(keepWavIdx);
        noise = noise(keepWavIdx);
        nRemoved = length(keepWavIdx) - sum(keepWavIdx);
        disp(['spike check: removed ', num2str(nRemoved), ' of ',...
            num2str(length(keepWavIdx)),' traces'])
        totalnRemoved = totalnRemoved + nRemoved;
        ntraces = ntraces + sum(keepWavIdx);
        
        if ~isempty(signal) && ~isempty(noise) 
            signalSpec = amplitude_spectrum(signal);
            noiseSpec = amplitude_spectrum(noise);
        %    spec=plot_spectrum(get(c.Z.LTK.corr,'waveforms'));
    %        close(gcf);
            clear normSNRAmp normSignalAmp normNoiseAmp signalAmp...
                noiseAmp SNRAmp freqCell;
            freqCell={signalSpec(:).f};
            if any(cellfun(@length,freqCell) ~= length(signalSpec(1).f))
               Freq = linspace(0,25,10000);
            else
              Freq = signalSpec(1).f;
            end

            signalAmp = zeros(length(signalSpec), length(Freq));
            noiseAmp = signalAmp; SNRAmp = signalAmp;
            for sp=1:1:length(signalSpec)
               %normalized the amplitudes and store them in a matrix for
               %cross-event averaging. each one is 
               normSNRAmp = signalSpec(sp).amp ./ noiseSpec(sp).amp;
               normSignalAmp = signalSpec(sp).amp;%/max(signalSpec(sp).amp);
               normNoiseAmp = noiseSpec(sp).amp;%/max(noiseSpec(sp).amp);

               if any(cellfun(@length,freqCell) ~= length(signalSpec(1).f))
                   signalAmp(sp,:) = interp1(signalSpec(sp).f,...
                       normSignalAmp, Freq);
                   noiseAmp(sp,:) = interp1(noiseSpec(sp).f,...
                       normNoiseAmp, Freq);
                   SNRAmp(sp,:) = interp1(signalSpec(sp).f, normSNRAmp,...
                       Freq);
               else
                  signalAmp(sp,:) = normSignalAmp;
                  noiseAmp(sp,:) = normNoiseAmp;
                  SNRAmp(sp,:) = normSNRAmp;
               end
               %averageSignalSpectrumHandle(ncl) = plot(specAx,Freq,...
               %    signalAmp(sp,:),'color',cmap(ncl,:),'linewidth',0.1);
               %averageNoiseSpectrumHandle(ncl) = plot(specAx,Freq,...
               %    noiseAmp(sp,:),':','color',cmap(ncl,:),'linewidth',0.1);
               
               if nanmean(nanmean(noiseAmp)) > 10
                   disp(['cluster ',num2str(ncl)])
                   disp(['station ',cstation{s}])
                   disp(['sp= ', num2str(sp)])
                   error('noise too high?!')
               end
            end
            
            averageSignalSpectrumPerStation.(cstation{s}).Amp = ...
                nanmean(signalAmp);
            averageSignalSpectrumPerStation.(cstation{s}).Freq = Freq;
            averageNoiseSpectrumPerStation.(cstation{s}).Amp = ...
                nanmean(noiseAmp);
            averageNoiseSpectrumPerStation.(cstation{s}).Freq = Freq;
            averageSNRSpectrumPerStation.(cstation{s}).Amp = ...
                nanmean(SNRAmp);
            averageSNRSpectrumPerStation.(cstation{s}).Freq = Freq;
%             averageSignalSpectrumHandle(ncl)=plot(specAx,Freq,nanmean(signalAmp),'--','color',cmap(ncl,:));
%             averageNoiseSpectrumHandle(ncl)=plot(specAx,Freq,nanmean(noiseAmp),':','color',cmap(ncl,:));
%             averageSNRSpectrumHandle(ncl)=plot(specAx,Freq,nanmean(SNRAmp),'color',cmap(ncl,:));
        end
    %         if any(cellfun(@length,averageSNRSpectrumPerStation.(cstation{:}).Freq)~=...
    %             length(averageSNRSpectrumPerStation.(cstation{1}).Freq))
    %        crossStationSignalAmp(k,:)=interp1(signalSpec(sp).f,normSignalAmp,Freq);
    %        crossStationNoiseAmp(k,:)=interp1(noiseSpec(sp).f,normNoiseAmp,Freq);
    %        crossStationSNRAmp(k,:)=interp1(signalSpec(sp).f,normSNRAmp,Freq);
    %     else
           crossStationSignalAmp{ncl}(s,:) = averageSignalSpectrumPerStation.(...
               cstation{s}).Amp;
           crossStationNoiseAmp{ncl}(s,:) = averageNoiseSpectrumPerStation.(...
               cstation{s}).Amp;
           crossStationSNRAmp{ncl}(s,:) = averageSNRSpectrumPerStation.(...
               cstation{s}).Amp;
%            k=k+1;
    %     end
    end
    
%     averageSNRSpectrumHandle(ncl)=plot(specAx,Freq,nanmean(...
%         crossStationSNRAmp{ncl}),'color',cmap(ncl,:));
    ntracesPerCluster(ncl) = ntraces;
end

noisePerCluster = zeros(numel(clusterI), size(crossStationNoiseAmp{1},2));
for ncl=1:1:numel(clusterI)
    averageSignalSpectrumHandle(ncl) = plot(specAx,Freq,nanmean(...
        crossStationSignalAmp{ncl}),'color',cmap(ncl,:),'linewidth',0.5);
%     averageNoiseSpectrumHandle(ncl) = plot(specAx,Freq,nanmean(...
%        crossStationNoiseAmp{ncl}),'--','color',cmap(ncl,:),'linewidth',1);
    noisePerCluster(ncl,:) = nanmean(crossStationNoiseAmp{ncl});        
end

% noisePerCluster = zeros(numel(clusterI)-1, size(crossStationNoiseAmp{1},2));
% for ncl=1:1:numel(clusterI)-1
%     noisePerCluster(ncl,:) = nanmean(crossStationNoiseAmp{ncl});
% end

averageNoiseSpectrumHandle = plot(specAx,Freq,nanmean(...
        crossStationNoiseAmp{ncl}),'-','color','k','linewidth',0.5);


xlabel('Frequency')
ylabel('Relative Amplitude')
box on
specAx.XScale='log';
specAx.YScale='log';
xlim([0.1 25])
% legend([averageSignalSpectrumHandle,averageNoiseSpectrumHandle(1)],...
%     'Mantle Wedge','Interface','Intraslab','Noise Spectra');

% legend([averageSignalSpectrumHandle,averageNoiseSpectrumHandle(1)],...
%     '          M < 0.0', '0.0 <= M < 0.5', '0.5 <= M < 1.0',...
%     '1.0 <= M < 1.5', '1.5 <= M < 2.0', '2.0 <= M < 2.5', '2.5 <= M',...
%     'Noise Spectra');


legendText =  {'2.5 <= M < 3.0'; '2.0 <= M < 2.5'; '1.5 <= M < 2.0';...
     '1.0 <= M < 1.5'; '0.5 <= M < 1.0'; '0.0 <= M < 0.5';...
     '-0.5 <= M < 0.0'; 'Noise Spectrum'};
for ncl=1:1:numel(clusterI)
    legendText{ncl} = [legendText{ncl}, ', n=',...
        num2str(ntracesPerCluster(ncl)) ];
end
legendText{end} = [legendText{end}, ', n=', num2str(sum(ntracesPerCluster))];

legend([fliplr(averageSignalSpectrumHandle),averageNoiseSpectrumHandle(1)],legendText)

% legend([fliplr(averageSignalSpectrumHandle),averageNoiseSpectrumHandle(1)],...
%     '2.5 <= M', '2.0 <= M < 2.5', '1.5 <= M < 2.0', '1.0 <= M < 1.5',...
%     '0.5 <= M < 1.0', '0.0 <= M < 0.5',...
%     '           M < 0.0', 'Noise Spectrum');


% legend([averageSignalSpectrumHandle,averageNoiseSpectrumHandle(1)],...
%     '          M < 0.5', '0.5 <= M < 1.5', '1.5 <= M < 2.0', '2.5 <= M',...
%     'Noise Spectra');

% legend([fliplr(averageSignalSpectrumHandle),averageNoiseSpectrumHandle(1)],...
%      '2.5 <= M', '1.5 <= M < 2.0','0.5 <= M < 1.5','          M < 0.5',...
%     'Noise Spectra');

set(gcf,'color','w')

% print(gcf,'spectra_clusterEarthquakes','-dpdf','-painters','-bestfit')
print(gcf,'spectra_Earthquakes_perMagnitudeBin_mean','-dpdf','-painters','-bestfit')

