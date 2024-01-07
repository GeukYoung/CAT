function [RE] = FxEIT_RE(Data,mask)
% input
%   Data -> [vv x timeseries]
% output
%   RE.raw -> whole ch RE data
%     .min -> minimum RE of time series
%     .avg -> average RE of time series
%     .max -> maximum RE of time series

if size(Data,1) <= 8^2
    ch = 8;
elseif size(Data,1) <= 16^2
    ch = 16;
elseif size(Data,1) <= 32^2
    ch = 32;
end

if nargin < 2 
    mask = FxEIT_mask(ch); % dafault : adjacent
end

if (size(Data,1) == 8^2) || (size(Data,1) == 16^2) || (size(Data,1) == 32^2)
    temp_full = Data;
else
    idx_full = 1:ch^2;
    idx_full(mask) = [];
    temp_full = zeros(ch^2,size(Data,2));
    temp_full(idx_full,:) = Data;
end
    
for cnt = 1:size(Data,2)
    temp_data = temp_full(:,cnt);
    temp_data = reshape(temp_data,ch,ch);
    temp_RE = zeros(ch,ch);
    for j = 1:16
        for k = 1:16
            temp_RE(j,k) = abs(((temp_data(j,k) - temp_data(k,j))/((temp_data(j,k) + temp_data(k,j))/2)) * 100);
        end
    end
    temp_RE2 = temp_RE;
    temp_RE = reshape(temp_RE,ch^2,1);
    temp_RE(mask) = [];
    RE.RE_raw(:,cnt) = temp_RE;
   
    % visualize
    if 0
    temp_data2 = temp_full(:,cnt);
    temp_data2(mask) = 0;
    temp_data2 = reshape(temp_data2,ch,ch);
%     plot(temp_data2);
    for i = 1:16
        temp_data2(:,i) = circshift(temp_data2(:,i),-(i-1));
    end
    
    close;
    fig = figure; set(fig,'units','normalized' ,'outerposition' ,[0.15 0.1 0.5 0.35]);
    subplot(121);
    errorbar(mean(temp_data2'),std(temp_data2'),'-o');
    xlim([2 16]); title('Voltage Overlay Plot');
    xlabel('num'); ylabel('Demodulation data');
    
    subplot(122);
    temp_RE2 = reshape(temp_RE2,ch^2,1);
    temp_RE2(mask) = NaN;
    result_avgRE = nanmean(temp_RE2);
    temp_RE2 = reshape(temp_RE2,ch,ch);
    
    pcolor(temp_RE2); title('RE map'); xlabel('Inj num'); ylabel('Meas num');
    hold on; plot([0 ch],[0 ch],'-.r','linewidth',2); set(gca, 'YDir', 'reverse'); axis equal; axis([1 ch 1 ch]);
    colorbar; colormap(parula);
    disp(num2str(result_avgRE));
    end
end

RE.RE_min = min(RE.RE_raw);
RE.RE_max = max(RE.RE_raw);
RE.RE_avg = mean(RE.RE_raw);
