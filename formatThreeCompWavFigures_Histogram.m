function formatThreeCompWavFigures_Histogram(ax, stationORevent, plotComp,...
    plotType, plotArrivals, labelArrivals, fileNameAddition0)

    imagePaperA4 = false;
    for k=1:1:numel(plotComp)
        box(ax(k),'on');
        hAll = findall(ax(k));
        for idx = 1 : length(hAll)
          try
            set(hAll(idx),'LineWidth',0.02);
          catch
            % never mind...
          end
          try
           set(hAll(idx),'fontsize',7);
          catch
           % never mind...
          end
          
          try
              if ~contains(stationORevent,'event') &&...
                      ~contains(stationORevent,'200') && (contains(hAll(...
                      idx).String,'P') || contains(hAll(idx).String,'S'))
                  imagePaperA4 = true;
                      set(hAll(idx),'BackgroundColor','w');
                      set(hAll(idx),'EdgeColor','k');
                      set(hAll(idx),'Margin',1);
                      set(hAll(idx),'Margin',0);
              end
          catch
            % never mind
          end
          try
              if hAll(idx).Color == [0 0 1]
                set(hAll(idx),'Color','k');
              end
          catch
            % never mind...
          end
        end
    end
    for k=1:1:numel(plotComp)
        ax(k).Title.FontSize = 9;
        ax(k).XLabel.FontSize = 8;
        % make sure the title is in the middle
        ax(k).Title.Position(1) = sum(ax(k).XLim)/2;
        
        if k==1
            ax(k).YLabel.FontSize = 8;
            ax(k).YLabel.BackgroundColor = 'w';
            ax(k).YLabel.EdgeColor = 'k';
            ax(k).YLabel.Margin = 1;
        else
            ax(k).YTickLabel = {};
            ax(k).YLabel.String = '';
            
        end
    end
    drawnow;

    np = numel(plotComp);
    newax = ax;
    if np>1
        fig = figure;
        for k=1:1:numel(plotComp)
            newax(k) = copyobj(ax(k),fig);
            % with ylabels kept on each axis:
            % newax(k).Position(1) = 1/np * (k-1)+0.035;
            % newax(k).Position(3) = 1/np - 0.04;
            
            % with y label only kept on leftmost axis:
            newax(k).Position(1) = (1-0.085)/np * (k-1)+0.08;
            newax(k).Position(3) = (1-0.085)/np - 0.01;
            
%                     newax(k).Position(3) = 1/np - 0.1;

            newax(k).Title.FontSize = 9;
            newax(k).XLabel.FontSize = 8;
            newax(k).YLabel.FontSize = 8;
            newax(k).YLabel.FontSize = 8;
        end
    else
        fig = gcf;
    end

    set(fig, 'Color', 'w');
    set(fig,'Units','centimeters');
    set(fig,'PaperUnits','centimeters');
    pos = get(fig,'Position');
    
    %make a a4-page plot only for the station gathers
    if imagePaperA4
        set(fig,'Position',[pos(1), pos(2), 21.0, 29.7]);
        % set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3),pos(4)])
        set(fig,'PaperPositionMode','Auto','PaperSize', [21.0, 29.7] )
    else
        for j=1:1:numel(newax)
            %oldpos = newax(j).Position;
            %newpos = [0.08 oldpos(2) oldpos(3) oldpos(4)-0.08];
            %newax(j).Position = newpos;
            newax(j).Position(2) = 0.05;
            newax(j).Position(4) = newax(j).Position(4)-0.09;
            ybounds = get(newax(j),'YLim');
            newax(j).Title.Position(2) = ybounds(1) - diff(ybounds)*0.11;
        end
        set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3),pos(4)])
    end

    % chose the file name by station/event, components, processing steps
    % etc.
    compName = char(plotComp);
    compName = reshape(compName', 1, numel(compName));
    drawnow;

    if contains(compName,'p')
        fileNameAddition1 = 'PolFiltered_';
    else
        fileNameAddition1 = 'Filt_';
    end
    
    if plotArrivals && labelArrivals
        fileNameAddition2 = 'wArrivalz_Labels_';
    elseif plotArrivals
        fileNameAddition2 = 'wArrivals_';
    else
        fileNameAddition2 = '';
    end

    print(fig,['Figures/', fileNameAddition0, fileNameAddition1,...
        fileNameAddition2, '_', plotType, '_', stationORevent, '_', compName, ...
        '.pdf'], '-dpdf', '-painters', '-r600')
    
    

% for j=1:1:numel(ax)
%     if mod(j,2)==0
%         fig = figure;
%         ax1 = copyobj(ax(j-1),fig);
%         ax2 = copyobj(ax(j),fig);
%         ax2.Position(1) = 0.55;
%         
%         ax1.Position(3) = 0.42;
%         ax2.Position(3) = 0.42;
%     end
% end