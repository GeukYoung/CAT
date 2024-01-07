function [locs_Rwave,mRRI,ECG_freq,ref] = FxEIT_findRpeak(ECG_data,Fs,th_peak)
% input
%   ECG_data : ECG raw data
%   Fs : ECG sampling rate
%   th_peak(optional) : thresh hold level of ECG peak detection
% output
%   locs_Rwave : R peak location
%   RR_interval : RR interval
%   ECG_freq : fundamental ECG frequncy
%   ref : generate ECG reference signal
if size(ECG_data,1) == 1
    ECG_data = ECG_data';
end

% detrend data
% [p,~,mu] = polyfit((1:numel(ECG_data))',ECG_data,6);
% f_y = polyval(p,(1:numel(ECG_data))',[],mu);
% Wn = 0.5;
% [b,a] = butter(5, 2*Wn/Fs,'high');
% ECG_data2 = filtfilt(b,a,ECG_data);
% 
% plot(ECG_data); hold on;
% plot(ECG_data2,'r');
% temp = movmean(ECG_data,Fs); plot(temp);
% plot(ECG_data2)
[up, ~] = envelope(ECG_data,Fs*2,'rms');
Detrend_ECG = ECG_data./up;
figure;
subplot(311);
plot(ECG_data); ylim([-500 500]); title('Raw ECG')
subplot(312);
plot(up); ylim([-150 150]); title('RMS Envelope');
subplot(313); 
plot(Detrend_ECG); title('After Normalized')
% Detrend_ECG = ECG_data - f_y;        
% Detrend_ECG = Detrend_ECG.^2;
if nargin < 3
    th_peak = 2*std(Detrend_ECG);
end

% find peak
[~,locs_Rwave] = findpeaks(Detrend_ECG,'MinPeakHeight',th_peak,...
    'MinPeakDistance',round(Fs/3));
mRRI = mean(locs_Rwave(2:end)-locs_Rwave(1:end-1));

figure; plot(Detrend_ECG);
hold on; plot(locs_Rwave,Detrend_ECG(locs_Rwave),'rv'); hold off;

% generate rectangle pulse wave
ref = zeros(length(Detrend_ECG),1);
ref(locs_Rwave) = 1;

temp = ref;
for i = 1:round(mRRI/2)
    temp = temp + circshift(ref,i);
end
ref = temp';
ref(1:locs_Rwave(1)) = 0;
ref(ref>1) = 1;

figure; 
plot(Detrend_ECG);
hold on; plot(locs_Rwave,Detrend_ECG(locs_Rwave),'kv');
% plot(ref*(max(Detrend_ECG)-min(Detrend_ECG))+min(Detrend_ECG),'r'); hold off;
ECG_freq = length(locs_Rwave)/length(ECG_data)*Fs;