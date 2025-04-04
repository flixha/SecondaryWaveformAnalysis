function firstArrival = getFirstArrivalTimeFromFM3D(c, arrivals,...
    sta, comp, traceOrder, evID, arrivalLine)

    firstArrivalPhase = char(arrivals(evID).arrivals.(sta).phase(1));
    firstArrivalPhases = firstArrivalPhase(1:3:end)';
    if strcmp(firstArrivalPhases,repmat('P',...
            length(firstArrivalPhases),1))
        firstArrival = arrivals(evID).arrivals.(sta).time(1);
    else
        % find an alternative first arrival from neighboring
        % events
        warning(['Earliest arrival is NOT all P-phases,',... 
            'did FM3D find no path?']);

        jm1 = arrivalLine;
        firstArrM1 = nan;
        while jm1 > 1
            jm1 = jm1 - 1;
            evIDm1 = c.(comp).(sta).cat.table.EventID(traceOrder(jm1));

            if ~isempty(arrivals(evIDm1).arrivals)
                if ~isempty(arrivals(evIDm1).arrivals.(sta).time)
                    firstArrivalPhase = arrivals(...
                        evIDm1).arrivals.(sta).phase(1);
                    firstArrivalPhases = firstArrivalPhase(1:3:end)';
                    if strcmp(firstArrivalPhases,repmat('P',...
                            length(firstArrivalPhases),1))
                        firstArrM1 = arrivals(evIDm1).arrivals.(sta).time(1);
                        continue
                    end
                end
            end
        end

        jp1 = arrivalLine;
        firstArrP1 = nan;
        while jp1 < c.(comp).(sta).cat.numberOfEvents - 1
            jp1 = jp1 + 1;
            evIDp1 = c.(comp).(sta).cat.table.EventID(traceOrder(jp1));

            if ~isempty(arrivals(evIDp1).arrivals)
                if ~isempty(arrivals(evIDp1).arrivals.(sta).time)
                    firstArrivalPhase = arrivals(...
                        evIDp1).arrivals.(sta).phase(1);
                    firstArrivalPhases = firstArrivalPhase(1:3:end)';
                    if strcmp(firstArrivalPhases,repmat('P',...
                            length(firstArrivalPhases),1))
                        firstArrP1 = arrivals(evIDp1).arrivals.(sta).time(1);
                        continue
                    end
                end
            end
        end                                      

        firstArrival = nanmean([firstArrM1, firstArrP1]);
    end
    
    
    
    %%%
%     myPhases = cell(0,0);
%     for k=1:1:length(cstation)
%         for j=1:1:length(arrivals)
%             if ~isempty(arrivals(j).arrivals)
%                 if ~isempty(arrivals(j).arrivals.(sta).phase)
%                     sta = char(cstation(k));
%                     myPhases(j,k) = arrivals(j).arrivals.(sta).phase(1);
%                 end
%             end
%         end
%     end
%     descLengths = cellfun(@length,myPhases);
%     [sortedLengths, order] = sort(descLengths);
%     myPhasesSorted = myPhases(order(:,1),:);
