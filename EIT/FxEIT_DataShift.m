function [idx_shift] = FxEIT_DataShift(nshift,ch,idx_mask)
if nargin < 3
    mask = FxEIT_mask(ch);
elseif nargin < 2
    ch = 16;
    mask = FxEIT_mask(ch);
else
    mask = FxEIT_mask(ch,idx_mask);
end

idx_shift = reshape(1:ch^2,ch,ch);
for i = 1:nshift
    idx_shift = idx_shift([ch 1:ch-1],[ch 1:ch-1]);
end

idx_shift = reshape(idx_shift,ch^2,1);
for i = 1:length(mask)
    idx_shift(idx_shift == mask(i)) = [];
end

cnt = 1;
for i = 1:ch^2
    if ismember(i,idx_shift) == true
        idx_shift(idx_shift == i) = cnt;
        cnt = cnt + 1;       
    end
end

