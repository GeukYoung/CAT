function [idx] = FxFindIdx(timeinfo,tp)
if isdatetime(timeinfo(1))
    if isdatetime(tp)
        [~, idx] = min(abs(timeinfo-tp));
    elseif tp <= 1
        [h,m,s] = hms(timeinfo);
        timeinfo_n = (h*60*60 + m*60 + s)/(24*60*60); % 24h normalize
        [~, idx] = min(abs(timeinfo_n-tp));
    end
else
    [~, idx] = min(abs(timeinfo-tp));
end