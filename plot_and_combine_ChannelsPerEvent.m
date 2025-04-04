function  plot_and_combine_ChannelsPerEvent(eventC, eventIDs, p,...
    printFigure, plotEnvelope, baseScale, fileNameAddition0, arrivals)

    if plotEnvelope
        fileNameAddition0 = [fileNameAddition0, 'Envelope_'];
    end
    
    if size(p,1) > size(p,2)
        p = p';
    end
    
    for p=p       
        scale = baseScale;
        % Plot the three-component figures 3 times: 
        % 1. with raw ZRT channels,
        % 2. with polarization-filtered ZRT channels,
        % 3. with polarization-filtered ZRT channels plus theoretical arrivals
        % 4. with polarization-filtered ZRT channels, plus arrivals and labels
        switch p
            case 1
                plotArrivals = false;
                labelArrivals = false;
                plotComp = {'Z','R','T'};
                fileNameAddition0 = [fileNameAddition0, ''];
                plotEnv = false;
                plotType = 'wig'; linewidth = 0.5;
            case 2
                plotArrivals = true; labelArrivals = true;
                plotComp = {'Zp','Rp','Tp'};
                fileNameAddition0 = [fileNameAddition0, ''];
                plotEnv = false;
                plotType = 'wigbyy'; linewidth = 0.5; scale = 0;
            case 3
                plotArrivals = false; labelArrivals = false;
                plotComp = {'Zp','Rp','Tp'};
                fileNameAddition0 = [fileNameAddition0, ''];
                plotEnv = plotEnvelope; linewidth = 0.1;
                if plotEnvelope; plotType='bwig';else; plotType='wig'; end
            case 4
                plotArrivals = true; labelArrivals = false;
                plotComp = {'Zp','Rp','Tp'};
                fileNameAddition0 = [fileNameAddition0, ''];
                plotEnv = plotEnvelope; linewidth = 0.1;
                if plotEnvelope; plotType='bwig';else; plotType='wig'; end
            case 5
                plotArrivals = true; labelArrivals = true;
                plotComp = {'Zp','Rp','Tp'};
                fileNameAddition0 = [fileNameAddition0, ''];
                plotEnv = plotEnvelope; linewidth = 0.1;
                if plotEnvelope; plotType='bwig';else; plotType='wig'; end
        end
    
        if plotEnv
            scale = baseScale * 1.65;
        end

        eventNums = find(eventC.cat.table.EventID == eventIDs);

        for j=1:1:length(eventNums)
            % give a name based on the location of the event
            DFST = eventC.cat.table.DistFromSlabTop(eventNums(j));
            if DFST > 1
                fileNameAddition1 = ['MW_',fileNameAddition0];
                upperXLimit = 20;
            elseif abs(DFST) < 1
                fileNameAddition1 = ['IF_',fileNameAddition0];
                upperXLimit = 22;
            elseif DFST < -1 && DFST > -8
                fileNameAddition1 = ['SC_',fileNameAddition0];
                upperXLimit = 24;
            elseif DFST < -30
                fileNameAddition1 = ['WBZ_lowerPlane_',fileNameAddition0];
                upperXLimit = 30;
            else
                fileNameAddition1 = ['SM_',fileNameAddition0];
                upperXLimit = 25;
            end

            
            if eventC.cat.table.DistFromSlabTop(eventNums(j)) > 1
                maxPhases = 4;
            else
                maxPhases = 4;
            end
            
            for k=1:1:numel(plotComp)
                plotWavesOfEventSortedBy(eventC, eventNums(j), plotComp{k},...
                    'X', plotType,'scale', scale, 'linewidth',...
                    linewidth, 'markArrivalTimes', plotArrivals,...
                    'labelarrivals', labelArrivals, 'arrivals', arrivals,...
                    'plotEnvelope', plotEnv, 'maxconversions', 1,...
                    'maxreflections', 1, 'maxPhaseLength', maxPhases);

                xlim([-2 upperXLimit]);
                oldPos = get(gcf,'Position');
                set(gcf,'Position',[1+(j-1)*230, oldPos(2), 350, oldPos(4)])
                ax(k) = gca();
            end

            if printFigure
                %eventString = 'LowerPlatenDSZ';
                %eventString = ['eventID_', num2str(...
                %    eventC4.cat.table.EventID(j))];
                eventString = eventC.cat.table.yyyy_mm_dd(eventNums(j),:);
                formatThreeComponentWaveformFigures(ax, eventString,...
                    plotComp, plotType, plotArrivals, labelArrivals, ...
                    fileNameAddition1);
                close all;
            end
        end
    end