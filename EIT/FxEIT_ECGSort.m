function [ECG_data2, ECG_index] = FxEIT_ECGSort(ECG_data,EIT_stat,intp_num)
% input
%   ECG_data : ECG data from biopac
%   EIT_stat : EIT trigger signal from biopac
%   intp_num : Interpolation number
% output
%   ECG_data2 : matched data with EIT sampling rate

if size(ECG_data,1)>size(ECG_data,2)
    ECG_data = ECG_data';
end

flag = 1;
cnt = 0;
% EIT_stat = EIT_stat-min(EIT_stat);
EIT_stat = EIT_stat./max(EIT_stat);
for i = 1:length(EIT_stat)
    if (EIT_stat(i) > 0.5) && (flag == 1)
        cnt = cnt + 1;
        flag = 0;
        roc(cnt) = i;
    elseif (EIT_stat(i) < 0.5) && (flag == 0)
        cnt = cnt + 1;
        flag = 1;
        roc(cnt) = i;
    end
        EIT_index(i) = cnt;
end
roc(end-10:end) = [];

if nargin < 3
    intp_num = 1;
end

if intp_num > roc(2)-roc(1)
    disp('intpnum is lager than fs');
end

temp = roc(2:end) - roc(1:end-1);
plot(temp); ylim([0 max(temp)]);
% figure; plot(EIT_stat); hold on;
% plot(roc+5,EIT_stat(roc+5),'r.');
% ECG_index = linspace(roc(1),roc(end),length(roc)*intp_num);
ECG_index = roc;
ECG_data2 = ECG_data(:,round(ECG_index));