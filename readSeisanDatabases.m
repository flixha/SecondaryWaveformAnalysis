function seisCatalog = readSeisanDatabases(SEISAN_TOP, seisanREAdir,...
    loadWaveforms, seisanWAVdir, startTime, endTime, varargin)

    if ~isempty(varargin) > 0
        cstation = varargin(1);
    end

    nDataBases=length(seisanREAdir);
    
    SECPERDAY = 60*60*24;

    %% list of all possible stations
    % l=0;
    % for k=[1:1:42, 104, 124, 126, 129]
    %     l=l+1;
    %     cstation{l}=['S',sprintf('%03d',k)];
    % end
    % for k=[1:5]
    %     l=l+1;
    %     cstation{l}=['AT',sprintf('%02d',k)];
    % end
    % for k=[1:11]
    %     l=l+1;
    %     cstation{l}=['PE',sprintf('%02d',k)];
    % end
    % cstation=[cstation, 'IDHR', 'ITM', 'KEAI','KEK', 'LIT', 'LKD', 'LTK', 'SOH', 'KYTH', 'SYRO', 'SERI'];
%         cstation={'S010','S011','S012','S013','S014','S015','S016','PE02','PE05','PE07','KLV','LTK','VLX'};

    %Or only the most suitable stations for conversion/reflection analysis
    % l=0;
    % for k=[10:1:14,16]
    %     l=l+1;
    %     cstation{l}=['S',sprintf('%03d',k)];
    % end
    % for k=[2,5,7]
    %     l=l+1;
    %     cstation{l}=['PE',sprintf('%02d',k)];
    % end



    % allStaCode=cell(0,1);
    % for j=1:1:seisCatalog.numberOfEvents
    %     for k=1:1:length(seisCatalog.arrivals{j,1})
    %         allStaCode=[allStaCode; seisCatalog.arrivals{j,1}(k).stacode];
    %     end
    % end
    % cstation=unique(allStaCode)';

    %% now load in all seisan s-files and waveforms and combine different 
    % Seisan-databases to one catalog "seisCatalog"

    for j=1:1:nDataBases
        eventDB{j} = Catalog();
        eventDB{j} = Catalog.retrieve('seisan', ...
            'dbpath', seisanREAdir{j}, ...
            'startTime', startTime, ....
            'endTime', endTime)
        %load in the waveform files
        % Check whether eventDB has any fields at all:
        % if ~isfield(eventDB{j},'numberOfEvents')
        %    error(['No catalog read in, no events found in ', seisanREAdir{j}]);
        %end

        if loadWaveforms
            for k=1:1:eventDB{j}.numberOfEvents
                year = datestr(eventDB{j}.otime(k),'yyyy');
                month = datestr(eventDB{j}.otime(k),'mm');

                if length(eventDB{j}.wavfiles{k}) == 0
                    warning('No waveform file given in event %s', datestr(eventDB{j}.otime(k), 'yyyy-mm-dd HH:MM:ss'))
                    continue
                end
                wavfilepath = fullfile(seisanWAVdir, year, month, char(eventDB{...
                    j}.wavfiles{k}));
                ds2 = datasource('miniseed', wavfilepath);
                % ds2 = datasource('seisan', wavfilepath);
                scnl = scnlobject('*', '*');
                starttime = datenum(startTime);
                endtime = datenum(endTime);
                %  try
                %      eventDB{j}.waveforms{k,1} = waveform(ds1, scnl,...
                %      starttime, endtime)'; 
                %  catch ME

                % Correct the channel tags to three characters (for some reason
                % they are all 4...)
                wav = waveform(ds2, scnl, starttime, endtime);

                % Throw away all waveforms that do not belong to requested
                % correlation-stations
                sta = get(wav,'station');
                
                % Check if only a subset of stations was requested
                if exist('cstation','var')
                    isRequestedStation = contains(sta,cstation);
                    wav = wav(isRequestedStation);
                end

                % shorten the traces a bit
                cut_start = eventDB{j}.otime(k);
                cut_end = eventDB{j}.otime(k) + 150/SECPERDAY;
                % Only if there is more than 120 /XX seconds of waveform
                for m=1:1:numel(wav)
                   if get(wav(m),'start') < cut_start &&...
                           get(wav(m),'end') > cut_end
                       wav(m) = extract(wav(m), 'TIME', cut_start, cut_end);
                   end
                end

                chn = get(wav,'channel');
                for l=1:1:length(chn)
                    if length(chn{l}) == 4
                        setChn = [chn{l}(1:2), chn{l}(4)];
                        wav(l) = set(wav(l),'channel',setChn);
                    end
                end
                %chn = char(chn);
                %wav = set(wav,'channel',chn);
                eventDB{j}.waveforms{k,1} = wav;
                fclose all;
            end
        end
        
        
        eventDB{j}.table.DateString = datestr(eventDB{j}.otime);
        eventDB{j}.table.ExistsInCatalogs = false(eventDB{j}.numberOfEvents,1);
    end


    % now add the events to one single database, but check whether they already
    % appear in that complete one before accepting them
    for j=1:1:nDataBases
        if j==1
            seisCatalog = eventDB{1};
        else
            %check whether in catalog
            for e=1:1:eventDB{j}.numberOfEvents
                %if in time period of one of the catalogs or itself
                timeBeforeEvent=seisCatalog.otime - detTimeWindow/2;
                timeAfterEvent=seisCatalog.otime + detTimeWindow/2;

                if any(eventDB{j}.otime(e) > timeBeforeEvent & eventDB{...
                        j}.otime(e) < timeAfterEvent)
                    eventDB{j}.table.ExistsInCatalogs(e)=true;
                end
            end

            selectEventsFromCatalogToCopy = (eventDB{...
                j}.table.ExistsInCatalogs==false);
            appendDB = Catalog();
            appendDB.table = eventDB{j}.table(...
                selectEventsFromCatalogToCopy,:);
            if loadWaveforms
                appendDB.waveforms = eventDB{j}.waveforms(...
                    selectEventsFromCatalogToCopy);
            end
            seisCatalog.table = [seisCatalog.table; appendDB.table];
            seisCatalog.waveforms = [seisCatalog.waveforms;...
                appendDB.waveforms];
        end
    end

 

    seisCatalog.table.Date = datetime(seisCatalog.otime, 'ConvertFrom', 'datenum');
    seisCatalog.table.DateString = datestr(seisCatalog.otime);
    seisCatalog.table.DateNumber = datenum(seisCatalog.otime);
    %sort the catalog by event origin time, and order waveforms in the same
    % way
    [seisCatalog.table, rowIndex] = sortrows(seisCatalog.table,...
        'otime','ascend');
    if loadWaveforms
        seisCatalog.waveforms = seisCatalog.waveforms(rowIndex);
    end
    %add an event index for simple identification of events, ascending with
    %time starting from event 1
    seisCatalog.table.EventID = [1:1:seisCatalog.numberOfEvents]';

    %% Calculate the actual distance from the slab top, not the depth!
    %loads longrid, latgrid and slab_top_surf
    %tic
    % load('slab_top_surf_v2.0.mat');
    slabtop_csv = load("slab_top_full3_Antilles_latlondep.csv");
    longrid = slabtop_csv(:, 1);
    latgrid = slabtop_csv(:, 2);
    slab_top_surf = slabtop_csv(:, 3);
    meanlon = mean(longrid);
    meanlat = mean(latgrid);
    disty = ll2utm(meanlat, meanlon) - ll2utm(meanlat + 1, meanlon);
    distx = ll2utm(meanlat, meanlon) - ll2utm(meanlat, meanlon + 1);
    distx = abs(distx(2)) / 10;
    disty = abs(disty(1)) / 10;
    %toc
    %tic
    for p=1:1:seisCatalog.numberOfEvents
        % 89 and 111 is the correction for distance in lat and lon-direction
        %calculate the distance from each lat/lon-gridnode grid
        dist_from_eq = sqrt(((longrid - seisCatalog.lon(p)) * distx) .^ 2+...
                         ((latgrid - seisCatalog.lat(p)) * disty) .^ 2+...
                         (slab_top_surf - seisCatalog.depth(p)) .^ 2).* ...
                    sign(slab_top_surf - seisCatalog.depth(p));
        [minEqDistFromSlabTop minEqDistFromSlabTopI] = min(abs(dist_from_eq(:)));
        eqDistFromSlabTop(p) = dist_from_eq(minEqDistFromSlabTopI);
    end
    %toc
    seisCatalog.table.DistFromSlabTop = eqDistFromSlabTop';
    
end
