function eventCout = threeComponentProcessingForEvObj(eventC, n, J, K,...
    dt, width, adjustAlignment)

% this function creates three-component objects for each of the records in
% the network-correlation object, applies rotation and a polarization 
% filter to enhance body wave arrivals.
% c2 is a structure with the substructure: Channel, cstation,
% correlation/catalog object. 
% cstation are the stations to be processed. n, J, K, dt, width are the
% parameters for polarizationfilter.

% now parallelizing over events, but function is still not the quickest for
% larger datasets

secPerDay = 60*60*24;
dayPerSec = 1/secPerDay;

eventC.Zp = eventC.Z;
eventC.R = eventC.Z;
eventC.Rp = eventC.Z;
eventC.T = eventC.Z;
eventC.Tp = eventC.Z;
eventC.TC = eventC.Z;


%for each event
for j=1:1:eventC.cat.numberOfEvents
    eventC.TC{j} = repmat(threecomp, size(get(eventC.Z{j},'waveforms')));

    Z = get(eventC.Z{j}, 'waveforms');
    N = get(eventC.N{j}, 'waveforms');
    E = get(eventC.E{j}, 'waveforms');
    Ztrigger = get(eventC.Z{j}, 'trig');
    nw = length(Z);
    
    % if no waveforms are left, then set all components to empty waveform
    % arrays
    if nw==0
        eventC.R{j} = set(eventC.R{j},'waveforms',Z);
        eventC.T{j} = set(eventC.T{j},'waveforms',Z);
        eventC.Zp{j} = set(eventC.Zp{j},'waveforms',Z);
        eventC.Rp{j} = set(eventC.Rp{j},'waveforms',Z);
        eventC.Tp{j} = set(eventC.Tp{j},'waveforms',Z);
        continue
    end
    
    % check that all waveforms have same length across components
    for k=1:1:nw
        lZ = length(get(Z(k),'data'));
        lN = length(get(N(k),'data'));
        lE = length(get(E(k),'data'));
        if lZ ~= lN || lZ ~= lE
            warning(['Components for event', datestr(eventC.cat.otime(j),...
                'yyyymmddHHMM'), ' station ', get(Z(k),'station'),...
                ' are not the same length, cutting all three to the ',...
                'shortest waveform.'])
            [shortest, shortestI] = min([lZ, lN, lE]);
            if shortestI == 1
                starttime = get(Z(k),'start');
                endtime = get(Z(k),'end');
            elseif shortestI == 2
                starttime = get(N(k),'start');
                endtime = get(N(k),'end');
            elseif shortestI == 3
                starttime = get(E(k),'start');
                endtime = get(E(k),'end');
            end
            Z(k) = extract(Z(k),'TIME', starttime, endtime);
            Z(k) = set(Z(k),'start',starttime);
            N(k) = extract(N(k),'TIME', starttime, endtime);
            N(k) = set(N(k),'start',starttime);
            E(k) = extract(E(k),'TIME', starttime, endtime);
            E(k) = set(E(k),'start',starttime);
        end     
    end
    

    % BackAzimuth from great circle arc between station and
    % hypoDD-location
    stations = get(Z,'station');
    evcat = eventC.cat.table(j,:);
    backazimuth = backazmiuthStationToEvent(stations, evcat);
    
    % BackAzimuth from Seisan Calculations
    % Zazimuth = c2.Z.(cstation{k}).cat.arrivals{j}.p_caz;
    % baz = Zazimuth-180;
    % backazimuth2(j,1) = baz;
    
    % incident = c2.Z.S010.cat.arrivals{1}.p_aofinc;
        
    R = Z;
    T = Z;
    Zp = Z;
    Rp = Z;
    Tp = Z;
    
    % create dummy array of three-comp objects
    zrttc = threecomp([Z(1),N(1),E(1)],0,0,[0 0 0 90 90 90]);
    zrttc = repmat(zrttc,numel(Z),1);
    %for each event
    % parfor k=1:1:nw
    for k=1:1:nw
        % Check if any channel is NaN or all zero, then don't rotate
        if any(isnan(get(Z(k),'data'))) || any(isnan(get(N(k),'data'))) ||...
                any(isnan(get(E(k),'data'))) ||...
                sum(get(Z(k),'data'))==0 || sum(get(N(k),'data'))==0 ||...
                sum(get(E(k),'data'))==0
            warning(['One component contains NaN or all zero, cannot ',...
                'rotate (setting R and T to NaN'])
            % fill all non-rotated/filtered components with zeros
            Zp(k) = Z(k);
            Zp(k) = set(Zp(k),'data',zeros(size((get(Zp(k),'data')))));
            
            R(k) = Z(k);
            chn = get(R(k),'ChannelTag');
            chn.channel = [chn.channel(1:end-1), 'R'];
            R(k) = set(R(k),'ChannelTag',chn);
            R(k) = set(R(k),'data',zeros(size((get(R(k),'data')))));           
            
            T(k) = Z(k);
            chn = get(T(k),'ChannelTag');
            chn.channel = [chn.channel(1:end-1), 'T'];
            T(k) = set(T(k),'ChannelTag',chn);
            T(k) = set(T(k),'data',zeros(size((get(T(k),'data')))));
            
            Rp(k) = set(R(k),'data',zeros(size((get(R(k),'data')))));
            Tp(k) = set(T(k),'data',zeros(size((get(T(k),'data')))));
        else
            W = [Z(k), N(k), E(k)];
            % Check for sub-sample time differences of channels
            st = get(W,'start');
            if st(1) ~= st(2) || st(1) ~= st(3)
                % Check that the difference is less than 10 % of one sample
                if st(1) - st(2) <= 0.1*1/get(W(1)*dayPerSec,'freq') &&...
                        st(1) - st(3) <= 0.1*1/get(W(1)*dayPerSec,'freq')
                    starttim = get(W(1),'start');
                    W(2) = set(W(2), 'start', starttim);
                    W(3) = set(W(3), 'start', starttim);
                end
            end 
            
            % now fill the required information on backazimuth, trigger
            % etc.
            trigger = [Ztrigger(k), Ztrigger(k), Ztrigger(k)]; 
                      
            orientation = [0 0 0 90 90 90];
            baz = backazimuth(k);
            tc = threecomp(W, trigger, baz, orientation);
            % set rotated records
            zrt = rotate(tc);
%             zrttc(j) = zrt;
            
%             zrt = particlemotion(zrt, 0.02, 0.8);
            % chose n = 0.5 or n = 1; dunno how much a difference it makes
            % zrt = polarizationfilter(zrt, 0.5, 1, 2, 0.02, 0.8);
            % J= 1, K=2
            % dt here as the sampling frequency. and width as twice the 
            % longest period
            % (or maybe better once to twice the DOMINANT period? -->)
            % looks like 1.5 times the longest period is a good compromise!
%             zrt = polarizationfilter(zrt, 1, 1, 2, 0.01, 0.6); % THIS
%             APPEARS LIKE SOME GOOD VALUES
%             zrt = polarizationfilter(zrt, 0.5, 1, 2, 0.01, 0.8);
            
            ww = get(zrt,'waveform');
            R(k) = ww(2);
            T(k) = ww(3);
            
            % set polarization-filtered records
            zrt = polarizationfilter(zrt, n, J, K, dt, width);
            ww = get(zrt,'waveform');
            Zp(k) = ww(1);
            Rp(k) = ww(2);
            Tp(k) = ww(3);
        end
    end
    eventC.TC{j} = zrttc;
    
    % compare back azimuths from Seisan vs. great-circle arc calculated ones
    % diff = backazimuth-backazimuth2;
    % histogram(diff,[-90:1:90]);
    
    
    eventC.R{j} = set(eventC.R{j},'waveforms',R);
    eventC.T{j} = set(eventC.T{j},'waveforms',T);
    
    eventC.Zp{j} = set(eventC.Zp{j},'waveforms',Zp);
    eventC.Rp{j} = set(eventC.Rp{j},'waveforms',Rp);
    eventC.Tp{j} = set(eventC.Tp{j},'waveforms',Tp);
    
    if adjustAlignment
        % Adjust triggers according to the polarizationfiltered-vertical
        % Component.
        % 1. run xcorr on first ~1.5 s of Zp-component
        echo off;
        eventC.Zp{j} = xcorr(eventC.Zp{j}, [-width/3 width]);
        if ~isreal(get(eventC.Zp{j},'corr'))
            warning(['Correlation matrix contains irrational values, ',...
                'taking only real'])
            eventC.Zp{j} = set(eventC.Zp{j},'corr',...
                real(get(eventC.Zp{j},'corr')));
        end
        
        % abort if there is only one waveform
        if length(get(eventC.Zp{j},'trig')) == 1
            continue
        end
        
        eventC.Zp{j} = linkage(eventC.Zp{j});
        % 2. Correct all components with the same adjustments from Zp-
        %    component
        Zplags = get(eventC.Zp{j},'lag');
        comp = fieldnames(eventC);
        % Don't try this with the catalog in the first field of eventC
        comp = comp(2:end);
        for k=1:1:numel(comp)
            % don't try this on the three-component object
            if ~strcmp(comp{k},'TC')
                eventC.(comp{k}){j} = set(eventC.(comp{k}){j},'lag',...
                    Zplags);
                eventC.(comp{k}){j} = adjusttrig(eventC.(comp{k}){j},...
                    'MMM', width/3);
                %0.06
        % 3. run xcorr again to update xcorr/lag fields etc.
                eventC.(comp{k}){j} = xcorr(eventC.(comp{k}){j},...
                    [-0.2 10]); 
                % check for irrational values...
                if ~isreal(get(eventC.(comp{k}){j},'corr'))
                    warning(['Correlation matrix contains irrational values, ',...
                        'taking only real'])
                    eventC.(comp{k}){j} = set(eventC.(comp{k}){j},'corr',...
                        real(get(eventC.(comp{k}){j},'corr')));
                end
                eventC.(comp{k}){j} = linkage(eventC.(comp{k}){j});
            end
        end
        echo on;
    end
    
    disp(['Three-component processing: ', num2str(...
        j/eventC.cat.numberOfEvents*100), ' % completed']);
end

eventCout = eventC;

% Sample code to plot waveforms and rectilinearity/weighting functions
% www = [get(zrt,'waveform'), get(zrt,'rectilinearity')];
% fig = figure();
% for j=1:1:length(www)
%     ax(j) = subplot(length(www),1,j);
%     plot(www(j),'axeshandle',ax(j))
% end
