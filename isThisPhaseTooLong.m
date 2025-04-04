function itis = isThisPhaseTooLong(thisPhase, maxPhaseLength,...
    distFromSlabTop)

    itis = true;
    
    % set as too long when it's too long
    if length(thisPhase) <= (maxPhaseLength * 3) - 2
        itis = false;
    end
    
    % be a bit less strict if it is very close by the subduction interface
    % - then allow one more phase in length
    if maxPhaseLength==4
        if distFromSlabTop > 0 && distFromSlabTop < 1
            if length(thisPhase) <= ((maxPhaseLength + 1) * 3) - 2
                itis = false;
            end
        end
        
        % also plot the S-reflection from the subd. Moho (contains no P)
        if length(thisPhase) <= ((maxPhaseLength + 1) * 3) - 2 &&...
                ~contains(thisPhase,'P')
            itis = false;
        end
    end
        
        