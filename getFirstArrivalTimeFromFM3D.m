function firstArrival = getFirstArrivalTimeFromFM3D(c, arrivals,...
    sta, comp, traceOrder, evID, arrivalLine, varargin)

    p = inputParser;

%     addRequired(p,'c',@isstruct);
%     addRequired(p,'arrivals',@isstruct);
%     addRequired(p,'sta',@ischar);
%     addRequired(p,'comp',@ischar);
%     addRequired(p,'traceOrder',@isnumeric);
%     addRequired(p,'evID',@isnumeric);
%     addRequired(p,'arrivalLine',@isnumeric);
    
    defaultPhase = 'P';
    validPhases = {'P','S'};
    checkPhases = @(x) any(validatestring(x,validPhases));
    addOptional(p,'Phase', defaultPhase, checkPhases);
    
%     parse(p, c, arrivals, sta, comp, traceOrder, evID, arrivalLine, varargin{:});
    parse(p, varargin{:});
    res = p.Results;
    Phase = res.Phase;
    
    [firstArrivalPhase, idx] = getFirstArrivalPhase(Phase, arrivals, evID, sta);
        
    firstArrivalPhases = firstArrivalPhase(1:3:end)';
    if strcmp(firstArrivalPhases,repmat(Phase,...
            length(firstArrivalPhases),1))
        firstArrival = arrivals(evID).arrivals.(sta).time(idx);
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
                    [firstArrivalPhase, idx] = getFirstArrivalPhase(Phase,...
                        arrivals, evIDm1, sta);
                    firstArrivalPhases = firstArrivalPhase(1:3:end)';
                    if strcmp(firstArrivalPhases,repmat(Phase,...
                            length(firstArrivalPhases),1))
                        firstArrM1 = arrivals(evIDm1).arrivals.(sta).time(idx);
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
                    [firstArrivalPhase, idx] = getFirstArrivalPhase(Phase,...
                        arrivals, evIDp1, sta);
                    firstArrivalPhases = firstArrivalPhase(1:3:end)';
                    if strcmp(firstArrivalPhases,repmat(Phase,...
                            length(firstArrivalPhases),1))
                        firstArrP1 = arrivals(evIDp1).arrivals.(sta).time(idx);
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

function [firstArrivalPhase, idx] = getFirstArrivalPhase(Phase, arrivals,...
    evID, sta)

    if strcmp(Phase,'P')
        firstArrivalPhase = char(arrivals(evID).arrivals.(sta).phase(1));
    elseif strcmp(Phase,'S')
        noP = ~contains(arrivals(evID).arrivals.(sta).phase,'P');
        phasesWithoutP = arrivals(evID).arrivals.(sta).phase(noP);
        firstArrivalPhase = char(phasesWithoutP(1));
    end
    idx = find(ismember(arrivals(evID).arrivals.(sta).phase,...
        firstArrivalPhase));
    
    