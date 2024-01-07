function [qrs_amp_raw,qrs_i_raw,delay,HR,ecg_SNR]=pan_tompkin2(ecg,fs,gr)
% real_time_draw = 0;
%% function [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(ecg,fs)
% Complete implementation of Pan-Tompkins algorithm

%% Inputs
% ecg : raw ecg vector signal 1d signal
% fs : sampling frequency e.gv 200Hz, 400Hz and etc
% gr : flag to plot or not plot (set it 1 to have a plot or set it zero not
% to see any plots
%% Outputs
% qrs_amp_raw : amplitude of R waves amplitudes
% qrs_i_raw : index of R waves
% delay : number of samples which the signal is delayed due to the
% filtering
%% Method :

%% PreProcessing
% 1) Signal is preprocessed , if the sampling frequency is higher then it is downsampled
% and if it is lower upsampled to make the sampling frequency 200 Hz
% with the same filtering setups introduced in Pan
% tompkins paper (a combination of low pass and high pass filter 5-15 Hz)
% to get rid of the baseline wander and muscle noise.

% 2) The filtered signal
% is derivated using a derivating filter to high light the QRS complex.

% 3) Signal is squared.4)Signal is averaged with a moving window to get rid
% of noise (0.150 seconds length).

% 5) depending on the sampling frequency of your signal the filtering
% options are changed to best match the characteristics of your ecg signal

% 6) Unlike the other implementations in this implementation the desicion
% rule of the Pan tompkins is implemented completely.

%% Decision Rule
% At this point in the algorithm, the preceding stages have produced a roughly pulse-shaped
% waveform at the output of the MWI . The determination as to whether this pulse
% corresponds to a QRS complex (as opposed to a high-sloped T-wave or a noise artefact) is
% performed with an adaptive thresholding operation and other decision
% rules outlined below;

% a) FIDUCIAL MARK - The waveform is first processed to produce a set of weighted unit
% samples at the location of the MWI maxima. This is done in order to localize the QRS
% complex to a single instant of time. The w[k] weighting is the maxima value.

% b) THRESHOLDING - When analyzing the amplitude of the MWI output, the algorithm uses
% two threshold values (THR_SIG and THR_NOISE, appropriately initialized during a brief
% 2 second training phase) that continuously adapt to changing ECG signal quality. The
% first pass through y[n] uses these thresholds to classify the each non-zero sample
% (CURRENTPEAK) as either signal or noise:
% If CURRENTPEAK > THR_SIG, that location is identified as a QRS complex
% candidate and the signal level (SIG_LEV) is updated:
% SIG _ LEV = 0.125 CURRENTPEAK + 0.875 SIG _ LEV

% If THR_NOISE < CURRENTPEAK < THR_SIG, then that location is identified as a
% noise peak and the noise level (NOISE_LEV) is updated:
% NOISE _ LEV = 0.125CURRENTPEAK + 0.875 NOISE _ LEV
% Based on new estimates of the signal and noise levels (SIG_LEV and NOISE_LEV,
% respectively) at that point in the ECG, the thresholds are adjusted as follows:
% THR _ SIG = NOISE _ LEV + 0.25  (SIG _ LEV ? NOISE _ LEV )
% THR _ NOISE = 0.5 (THR _ SIG)
% These adjustments lower the threshold gradually in signal segments that are deemed to
% be of poorer quality.


% c) SEARCHBACK FOR MISSED QRS COMPLEXES - In the thresholding step above, if
% CURRENTPEAK < THR_SIG, the peak is deemed not to have resulted from a QRS
% complex. If however, an unreasonably long period has expired without an abovethreshold
% peak, the algorithm will assume a QRS has been missed and perform a
% searchback. This limits the number of false negatives. The minimum time used to trigger
% a searchback is 1.66 times the current R peak to R peak time period (called the RR
% interval). This value has a physiological origin - the time value between adjacent
% heartbeats cannot change more quickly than this. The missed QRS complex is assumed
% to occur at the location of the highest peak in the interval that lies between THR_SIG and
% THR_NOISE. In this algorithm, two average RR intervals are stored,the first RR interval is
% calculated as an average of the last eight QRS locations in order to adapt to changing heart
% rate and the second RR interval mean is the mean
% of the most regular RR intervals . The threshold is lowered if the heart rate is not regular
% to improve detection.

% d) ELIMINATION OF MULTIPLE DETECTIONS WITHIN REFRACTORY PERIOD - It is
% impossible for a legitimate QRS complex to occur if it lies within 200ms after a previously
% detected one. This constraint is a physiological one  due to the refractory period during
% which ventricular depolarization cannot occur despite a stimulus[1]. As QRS complex
% candidates are generated, the algorithm eliminates such physically impossible events,
% thereby reducing false positives.

% e) T WAVE DISCRIMINATION - Finally, if a QRS candidate occurs after the 200ms
% refractory period but within 360ms of the previous QRS, the algorithm determines
% whether this is a genuine QRS complex of the next heartbeat or an abnormally prominent
% T wave. This decision is based on the mean slope of the waveform at that position. A slope of
% less than one half that of the previous QRS complex is consistent with the slower
% changing behaviour of a T wave  otherwise, it becomes a QRS detection.
% Extra concept : beside the points mentioned in the paper, this code also
% checks if the occured peak which is less than 360 msec latency has also a
% latency less than 0,5*mean_RR if yes this is counted as noise

% f) In the final stage , the output of R waves detected in smoothed signal is analyzed and double
% checked with the help of the output of the bandpass signal to improve
% detection and find the original index of the real R waves on the raw ecg
% signal

%% References :

%[1]PAN.J, TOMPKINS. W.J,"A Real-Time QRS Detection Algorithm" IEEE
%TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. BME-32, NO. 3, MARCH 1985.

%% Author : Hooman Sedghamiz
% Linkoping university
% email : hoose792@student.liu.se
% hooman.sedghamiz@medel.com

% Any direct or indirect use of this code should be referenced
% Copyright march 2014
%%
if ~isvector(ecg)
    error('ecg must be a row or column vector');
end


if nargin < 3
    gr = 1;   % on default the function always plots
end
ecg = ecg(:); % vectorize

%% Initialize

ax = zeros(1,6);

% figure;

%% Noise cancelation(Filtering) % Filters (Filter in between 5-15 Hz)
if fs == 200
    %% start figure
    %% remove the mean of Signal
    ecg = ecg - mean(ecg);
    %% Low Pass Filter  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
    %%It has come to my attention the original filter doesnt achieve 12 Hz
    %    b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
    %    a = [1 -2 1];
    %    ecg_l = filter(b,a,ecg);
    %    delay = 6;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Wn = 12*2/fs;
    N = 3; % order of 3 less processing
    [a,b] = butter(N,Wn,'low'); %bandpass filtering
    ecg_l = filtfilt(a,b,ecg);
    ecg_l = ecg_l/ max(abs(ecg_l));
    if gr
        ax(1)=subplot(321);plot(ecg);axis tight;title('Raw signal');
        ax(2)=subplot(322);plot(ecg_l);axis tight;title('Low pass filtered');
    end
    %% High Pass filter H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))
    %%It has come to my attention the original filter doesn achieve 5 Hz
    %    b = zeros(1,33);
    %    b(1) = -1; b(17) = 32; b(33) = 1;
    %    a = [1 1];
    %    ecg_h = filter(b,a,ecg_l);    % Without Delay
    %    delay = delay + 16;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Wn = 5*2/fs;
    N = 3; % order of 3 less processing
    [a,b] = butter(N,Wn,'high'); %bandpass filtering
    ecg_h = filtfilt(a,b,ecg_l);
    ecg_h = ecg_h/ max(abs(ecg_h));
    if gr
        ax(3)=subplot(323);plot(ecg_h);axis tight;title('High Pass Filtered');
    end
else
    %% bandpass filter for Noise cancelation of other sampling frequencies(Filtering)
    ecg(isnan(ecg)) = 0;
    f1=5; %cuttoff low frequency to get rid of baseline wander
    % f2=15; %cuttoff frequency to discard high frequency noise
    f2=25; %cuttoff frequency to discard high frequency noise
    Wn=[f1 f2]*2/fs; % cutt off based on fs
    N = 3; % order of 3 less processing
    [a,b] = butter(N,Wn); %bandpass filtering
    ecg_h = filtfilt(a,b,ecg);
    ecg_h = ecg_h/ max( abs(ecg_h));
    if gr
        ax(1) = subplot(3,2,[1 2]);plot(ecg);axis tight;title('Raw Signal');
        
        ax(3)=subplot(323);plot(ecg_h);axis tight;title('Band Pass Filtered');
    end
end
%% derivative filter H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))
if fs ~= 200
    int_c = (5-1)/(fs*1/40);
    b = interp1(1:5,[1 2 0 -2 -1].*(1/8)*fs,1:int_c:5);
else
    b = [1 2 0 -2 -1].*(1/8)*fs;
end

ecg_d = filtfilt(b,1,ecg_h);
ecg_d = ecg_d/max(ecg_d);

if gr
    ax(4)=subplot(324);plot(ecg_d);
    axis tight;
    title('Filtered with the derivative filter');
end
%% Squaring nonlinearly enhance the dominant peaks
ecg_s = ecg_d.^2;
if gr
    ax(5)=subplot(325);plot(ecg_s);axis tight;title('Squared');
end

%% Moving average Y(nt) = (1/N)[x(nT-(N - 1)T)+ x(nT - (N - 2)T)+...+x(nT)]
ecg_m = conv(ecg_s ,ones(1 ,round(0.150*fs))/round(0.150*fs));
delay = round(0.150*fs)/2;

if gr
    ax(6)=subplot(326);plot(ecg_m);
    axis tight;
    title('Averaged with 30 samples length,Black noise,Green Adaptive Threshold,RED Sig Level,Red circles QRS adaptive threshold');
    axis tight;
end

%% Main
op.mov_minmax = .50;
op.mov_mean = .2;
op.mov_th = .2*2;

% remove baseline
temp = movmean(movmin(ecg_m,round([fs fs]*op.mov_minmax)),[fs fs]*op.mov_mean); % figure; plot(ecg_m); hold on; plot(ecg_m-temp);
ecg_m = ecg_m - temp;
ecg_m = movmean(ecg_m,[3 3]);
ecg_m = ecg_m.^2;
% fine ecg R threshold
ecg_th = movmean(movmax(ecg_m,round([fs fs]*op.mov_minmax)),[fs fs]*op.mov_mean)*op.mov_th;
figure; plot(ecg_m); hold on; plot(ecg_th,'r');

% find peaks
% [pks,locs] = findpeaks(ecg_m);
% locs_orig = locs;
% locs = locs(pks>ecg_th(locs)); % remove small peak

tp_le = find(diff(ecg_m>ecg_th)>0);
tp_fe = find(-diff(ecg_m>ecg_th)>0);
if tp_le(1) > tp_fe(1)
    tp_fe(1) = [];
end
if length(tp_le) ~= length(tp_fe)
    tp_le(end) = [];
end

% figure; plot(ecg_m); hold on; plot(ecg_th); plot(tp_le,ecg_th(tp_le),'g^'); plot(tp_fe,ecg_th(tp_fe),'vr'); 
for cnt = 1:length(tp_le)
    [pks(cnt), locs(cnt)] = max(ecg_m(tp_le(cnt):tp_fe(cnt)));
    locs(cnt) = locs(cnt) + tp_le(cnt) - 1;
end
figure; plot(ecg_m); hold on; plot(locs, ecg_m(locs), 'rv','markerfacecolor','r');

% R peak selection
cnt = 4;
stat.hr = 80/60*fs;
stat.en_loc = false(length(locs),1);
stat.en_loc(1) = true;
% op.th_hr = .6;
op.th_hr = .2;
op.min_hr = 15/60*fs;
op.avg_hr = 4;
op.th_sdhr = 5 / 60 * fs; % 5 bpm variation
% figure; plot(ecg_m); hold on; plot(locs,ecg_m(locs),'kv');
while 1
    cnt = cnt + 1;
%     plot(locs(cnt),ecg_m(locs(cnt)),'rv','markerfacecolor','r');
    if sum(diff([locs(find(stat.en_loc,1,'last')) locs(cnt:cnt+1)]) < stat.hr(end)*op.th_hr)
        [~, tp] = min(abs([locs(cnt:cnt+1)-locs(find(stat.en_loc,1,'last'))]-stat.hr(end)));
        stat.en_loc(cnt + tp - 1) = true;
        stat.hr(cnt + tp - 1) = mean(diff(locs(find(stat.en_loc,op.avg_hr,'last'))));
%         plot(locs(cnt + tp - 1),ecg_m(locs(cnt + tp - 1)),'gv','markerfacecolor','g');
        cnt = cnt + 1;
    elseif diff([locs(find(stat.en_loc,1,'last')) locs(cnt)]) < op.min_hr
        % ignore
    else
        stat.en_loc(cnt) = true;
        stat.hr(cnt) = mean(diff(locs(find(stat.en_loc,op.avg_hr,'last'))));
%         plot(locs(cnt),ecg_m(locs(cnt)),'gv','markerfacecolor','g');
    end
    
    if std(diff(locs(cnt-4:cnt))) < op.th_sdhr
        stat.hr(cnt) = mean(diff(locs(cnt-4:cnt)));
    end
    
    if cnt >= length(locs)-1
        break;
    end
end

figure; 
h(1) = subplot(211);
plot(ecg_m); hold on; % plot(locs_orig,ecg_m(locs_orig),'kv','markerfacecolor','k');
plot(locs,ecg_m(locs),'rx','markerfacecolor','r');
plot(locs(stat.en_loc),ecg_m(locs(stat.en_loc)),'gv','markerfacecolor','g');
plot(ecg_th);
h(2) = subplot(212);
plot(locs(stat.en_loc(2:end)),diff(locs(stat.en_loc)),'.');
linkaxes(h,'x');

qrs_i_raw = locs(stat.en_loc) - round(delay);
qrs_i_raw(qrs_i_raw<fs) = [];
figure; plot(ecg); hold on; plot(qrs_i_raw,ecg(qrs_i_raw),'rx','markerfacecolor','r');
delay_window = round(delay*1.2);
for cnt = 1:length(qrs_i_raw)
%     [~, tp] = max(ecg(qrs_i_raw(cnt)-delay_window:qrs_i_raw(cnt)+delay_window));
%     qrs_i_raw(cnt) = qrs_i_raw(cnt) + (tp - delay_window - 1);
    
    [~, tp] = max(ecg(qrs_i_raw(cnt)-delay_window:qrs_i_raw(cnt)+delay_window));
    [~, tp2] = max(ecg(qrs_i_raw(cnt)+tp-2*delay_window-1:qrs_i_raw(cnt)+tp)-1);
    qrs_i_raw(cnt) = qrs_i_raw(cnt)+tp-2*delay_window-1 + tp2-1 ;
end
% figure; plot(ecg(qrs_i_raw(cnt)+tp-2*delay_window:qrs_i_raw(cnt)+tp))
qrs_amp_raw = ecg(qrs_i_raw);
plot(qrs_i_raw,ecg(qrs_i_raw),'gv','markerfacecolor','g');
% figure; hist(locs_new - qrs_i_raw)

HR = zeros(size(ecg));
HR_raw = zeros(size(ecg));
temp_HR = 1./(diff(qrs_i_raw)/fs)*60;

mov_size = 30;
th_lv = 0.2;
figure; plot(temp_HR); hold on;
plot(movmedian(temp_HR,[mov_size mov_size])*(1+th_lv));
plot(movmedian(temp_HR,[mov_size mov_size])*(1-th_lv));
upper_limit = movmedian(temp_HR,[mov_size mov_size])*(1+th_lv);
lower_limit = movmedian(temp_HR,[mov_size mov_size])*(1-th_lv);
for i = 1:length(temp_HR)
    if (temp_HR(i) > upper_limit(i)) || (temp_HR(i) < lower_limit(i))
        temp_HR(i) = 0.5 * (upper_limit(i) + lower_limit(i));
    end
end
temp_HR_mean = movmean(temp_HR, [5 5]);
plot(temp_HR_mean,'r');

for i = 1:length(temp_HR_mean)-1
    HR(qrs_i_raw(i):qrs_i_raw(i+1)) = temp_HR_mean(i);
    HR_raw(qrs_i_raw(i):qrs_i_raw(i+1)) = temp_HR(i);
end
HR(1:qrs_i_raw(1)) = temp_HR_mean(1);
HR_raw(1:qrs_i_raw(1)) = temp_HR(1);
HR = HR';
HR_raw = HR_raw';

ecg_SNR = abs(20*log10((abs((HR-HR_raw)+1)./HR)));
ecg_SNR = ecg_SNR([qrs_i_raw(1) qrs_i_raw(1:end-1)])';
qrs_i_raw = qrs_i_raw';

figure; 
h(1) = subplot(311);
plot(ecg); hold on; plot(qrs_i_raw,ecg(qrs_i_raw),'gv','markerfacecolor','g'); title('ECG R peak')
h(2) = subplot(312);
plot(HR); title('HR');
h(3) = subplot(313);
bar(qrs_i_raw,ecg_SNR); title('SNR');
linkaxes(h,'x');

end









