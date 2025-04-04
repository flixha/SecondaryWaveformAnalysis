

% difference in traveltimes for two fm3d grids of different resolution:
arr1 = loadFM3Darrivals('/Volumes/nasdata2/Documents2/Greece_MWcluster/FastMarching/Runs/RunSetup_0.05degGrid');
arr2 = loadFM3Darrivals('/Volumes/nasdata2/Documents2/Greece_MWcluster/FastMarching/Runs/RunSetup_0.02degGrid');

sta = fieldnames(arr1.arrivals);
diff = zeros(0,0);
for j = [1:1:numel(sta)]
    try
        thisdiff = arr1.arrivals.(char(sta(j))).time...
            - arr2.arrivals.(char(sta(j))).time;
        diff = [diff; thisdiff];
    catch
        continue
    end
end
    
figure
histogram(diff,[-20:0.05:20])
    
    