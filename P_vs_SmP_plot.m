% p vs s-p time

ptriggers = zeros(0,0);
striggers = zeros(0,0);
otimes = zeros(0,0);
m=1;
for j=1:1:greekevents.numberOfEvents
    arr = greekevents.arrivals{j};
    for k=1:1:numel(arr)
        ptrig = arr(k).p_time;
        strig = arr(k).s_time;
        origtime = greekevents.otime(j);
        if ~isempty(ptrig) & ~isempty(strig)
            ptriggers(m,1) = ptrig;
            striggers(m,1) = strig;
            origtimes(m,1) = origtime;
            m=m+1;
        end
    end
end

figure;
hold on;
secondsInDay = 60*60*24;
ptimes = (ptriggers - origtimes)*secondsInDay;
smptimes = (striggers - ptriggers)*secondsInDay;
plot(ptimes, smptimes,'.k')
