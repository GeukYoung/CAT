function [ data_filt ] = FxEIT_Filter( data, fs, N, Wn, type1, type2 )
% FxEIT_Filter( EIT_data, sampling, order, cut_off_freq, type1, type2 )
% type1 1:butter 2:cheby1
% type2 1:lowpass 2:highpass 3:bandpass 4.notch
% DataSet.EIT_filt = FxEIT_Filter(temp, fs1, 1, 6, 2, 1); % data, fs, N, Wn, type1, type2
% N = 4;
% Wn = [0.6 4];
% Wn = 4;
% type1 = 1;
% type2 = 2;
%% finite check
idx_notfin = ~isfinite(data);
data(idx_notfin) = NaN;
data = fillmissing(data','nearest','EndValues','nearest')';

switch type1
    case 1
        switch type2
            case 1
                [b,a] = butter(2, 2*Wn/fs,'low');
%                 freqs(b,a)
            case 2
                [b,a] = butter(N, 2*Wn/fs,'high');
            case 3
                [b,a] = butter(N, 2*Wn/fs,'bandpass');
            case 4
                [b,a] = butter(N, 2*Wn/fs,'stop');
            case 5 % filter design for remov swan (2.5 Hz)
                Fpass = 2;         % Passband Frequency
                Fstop = 2.5;         % Stopband Frequency
                Apass = 10;          % Passband Ripple (dB)
                Astop = 80;          % Stopband Attenuation (dB)
                match = 'stopband';  % Band to match exactly
                
                % Construct an FDESIGN object and call its BUTTER method.
                h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, fs);
                Hd = design(h, 'butter', 'MatchExactly', match);
        end
        if type2 == 5
            data_filt = filtfilt(Hd.sosMatrix, Hd.ScaleValues, data')';
        else
            data_filt = filtfilt(b,a,data')';
        end
    case 2
        [b a] = cheby1(N, 10, 2*Wn/fs);
        data_filt = filtfilt(b,a,data')';
    case 3
        ref1 = mean(data);
        ref1 = detrend(ref1);       % remove trend
        [ref1,]=remmean(ref1);      % remove mean
        ref1 = filtfilt(ones(30,1),1,ref1);  % mov avg filt
end

data_filt(idx_notfin) = NaN;

% if fs==50
%     data_length = 3000;
% else
%     data_length = 600;
% end
% if size(data,1)==1
%     figure(1); subplot(4,1,1); plot(data(1:data_length));
%     subplot(4,1,2); plot(data_filt(1:data_length));
%     subplot(4,1,3);
%     N = size(data,2);
%     k=0:N-1;                %create a vector from 0 to N-1
%     n_pad = 1;              % n' times zero padding
%     N = N*n_pad;            %get the number of points
%     freq = k/(N/fs);    %create the frequency range
%     for i = 1:size(data,1)
%         temp=fft(data(i,:),N)/N*2; % normalize the data
%         temp = temp(1:N/2); %take only the first half of the spectrum
%         stem(freq(1:N/2),abs(temp));
%         xlabel('Freq (Hz)'); ylabel('Amplitude');
%         axis([0 fs/2 0 0.15]); grid on;
%     end
%     clear i N k n_pad freq temp;
% 
%     subplot(4,1,4);
%     N = size(data_filt,2);
%     k=0:N-1;
%     n_pad = 1;
%     N = N*n_pad;
%     freq = k/(N/fs);
%     for i = 1:size(data_filt,1)
%         temp=fft(data_filt(i,:),N)/N*2;
%         temp = temp(1:N/2);
%         stem(freq(1:N/2),abs(temp));
%         xlabel('Freq (Hz)'); ylabel('Amplitude');
%         axis([0 fs/2 0 0.15]); grid on;
%     end
%     clear i N k n_pad freq temp;
% else
%     figure(1); subplot(4,1,1); plot(sum(data(:,1:data_length),1));
%     subplot(4,1,2); plot(sum(data_filt(:,1:data_length),1));
%     subplot(4,1,3);
%     x = sum(data,1);
%     N = size(x,2);
%     k=0:N-1;
%     n_pad = 1;
%     N = N*n_pad;
%     freq = k/(N/fs);
%     temp=fft(x,N)/N*2;
%     temp = temp(1:N/2);
%     stem(freq(1:N/2),abs(temp)); xlabel('Freq (Hz)'); ylabel('Amplitude');
%     axis([0 fs/2 0 500]); grid on;
%     clear i N k n_pad freq temp;
% 
%     subplot(4,1,4);
%     x = sum(data_filt,1);
%     N = size(x,2);
%     k=0:N-1;
%     n_pad = 1;
%     N = N*n_pad;
%     freq = k/(N/fs);
%     temp=fft(x,N)/N*2;
%     temp = temp(1:N/2);
%     stem(freq(1:N/2),abs(temp)); xlabel('Freq (Hz)'); ylabel('Amplitude');
%     axis([0 fs/2 0 500]); grid on;
%     clear i N k n_pad freq temp;
% end


end

