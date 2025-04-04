function arrivals = duplicateArrivalsForDoubleInstrumentSites(arrivals)



    for j=1:1:numel(arrivals)
        if isfield(arrivals(j).arrivals,'S028')
            if  ~isfield(arrivals(j).arrivals,'S036')
	            arrivals(j).arrivals.S036 = arrivals(j).arrivals.S028;
            end
            if  ~isfield(arrivals(j).arrivals,'S036')
                    arrivals(j).arrivals.S036 = arrivals(j).arrivals.S028;
            end
            if  ~isfield(arrivals(j).arrivals,'S041')
                    arrivals(j).arrivals.S041 = arrivals(j).arrivals.S028;
            end
            if  ~isfield(arrivals(j).arrivals,'S042')
                    arrivals(j).arrivals.S042 = arrivals(j).arrivals.S028;
            end
            
        end
        if isfield(arrivals(j).arrivals,'S038') &&...
                ~isfield(arrivals(j).arrivals,'ANDR')
            arrivals(j).arrivals.ANDR = arrivals(j).arrivals.S038;
        end
        if isfield(arrivals(j).arrivals,'ITHO')
            if ~isfield(arrivals(j).arrivals,'PE06')
                arrivals(j).arrivals.PE06 = arrivals(j).arrivals.ITHO;
            end
            if  ~isfield(arrivals(j).arrivals,'ITM')
                arrivals(j).arrivals.ITM = arrivals(j).arrivals.ITHO;
            end
        end
        if isfield(arrivals(j).arrivals,'S024') &&...
                ~isfield(arrivals(j).arrivals,'S124')
            arrivals(j).arrivals.S124 = arrivals(j).arrivals.S024;
        end
        
        if isfield(arrivals(j).arrivals,'N037') &&...
                ~isfield(arrivals(j).arrivals,'SOH')
            arrivals(j).arrivals.SOH = arrivals(j).arrivals.N037;
        end
        
        if isfield(arrivals(j).arrivals,'N039') &&...
                ~isfield(arrivals(j).arrivals,'SRS')
            arrivals(j).arrivals.SRS = arrivals(j).arrivals.N039;
        end
        
        if isfield(arrivals(j).arrivals,'AT01') &&...
                ~isfield(arrivals(j).arrivals,'KAT')
            arrivals(j).arrivals.KAT = arrivals(j).arrivals.AT01;
        end

    end
