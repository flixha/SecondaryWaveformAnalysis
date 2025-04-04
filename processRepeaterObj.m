function c = processRepeaterObj(c)

comp = fieldnames(c);
sta = fieldnames(c.(comp{1}));

for nc=1:1:length(comp)
    for s=1:1:length(sta)
        if isa(c.(comp{nc}).(sta{s}),'threecomp')
            % nothing
        elseif length(get(c.(comp{nc}).(sta{s}).corr, 'waveform')) > 1
            c.(comp{nc}).(sta{s}).corr = adjusttrig(c.(comp{nc}).(sta{s}...
                ).corr,'MMM',0.2);
            c.(comp{nc}).(sta{s}).corr = stack(c.(comp{nc}).(sta{s}).corr);
            
            %for t=1:1:width(c.(comp{nc}).(sta{s}).cat.table)
            fields = fieldnames(c.(comp{nc}).(sta{s}).cat.table);
            %c.(comp{nc}).(sta{s}).cat.table = [c.(comp{nc}).(sta{s}...
            %    ).cat.table; 
            newTableLine = c.(comp{nc}).(sta{s}).cat.table(1,:);
            for f=1:1:length(fields)
                try
                    % only if this is a non-integer value, not an index!
                    if any(round(c.(comp{nc}).(sta{s}).cat.table.(fields{...
                            f})) - c.(comp{nc}).(sta{s}).cat.table.(...
                            fields{f}) == 0)
                        continue
                    else
                        newTableLine.(fields{f}) = mean(c.(comp{nc}).(sta{s}...
                            ).cat.table.(fields{f}));
                    end
                catch
                    %
                end
            end
            c.(comp{nc}).(sta{s}).cat.table = [c.(comp{nc}).(sta{s}...
                ).cat.table; newTableLine];
                    

                    
        end
    end
end

