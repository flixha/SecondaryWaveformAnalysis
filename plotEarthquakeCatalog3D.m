function [figurehandle]=plotEarthquakeCatalog3D(seisanCatalog,markedEvents)
    %this function plots the supplied earthquakes from a catalog object in
    %3D
    modelFig=figure;
    figurehandle=modelFig;
    
    mF_ax1=axes;
    hold on;
    scatter3(greekevents.lon,greekevents.lat,greekevents.depth,2.^(greekevents.mag+4),'b','linewidth',1.5);
    scatter3(greekevents.lon(markedEvents),greekevents.lat(markedEvents),greekevents.depth(markedEvents),...
        2.^(greekevents.mag(markedEvents)+4),'r','linewidth',3);
    % scatter3(eqksLLD(:,5),eqksLLD(:,4),-eqksLLD(:,6),'k');
    % scatter3(greekevents.lon,greekevents.lat,greekevents.depth,'b','linewidth',1);
    % scatter3(greekevents.lon(eventI),greekevents.lat(eventI),greekevents.depth(eventI),'r','linewidth',2);

    mF_ax1.ZDir='reverse';
    mF_ax1.Box='on';mF_ax1.BoxStyle='full';
    mF_ax1.XGrid='on';mF_ax1.YGrid='on';mF_ax1.ZGrid='on';
    mF_ax1.DataAspectRatio=[1.1 1 111];
    xlim([20 25]); ylim([36 40]); zlim([0 180]);
    xlim([21.8 23]); ylim([37.4 38.4]); zlim([30 80]);
    xlim([22.38 22.9]); ylim([37.4 37.8]); zlim([30 80]);

    view(60,0)
end