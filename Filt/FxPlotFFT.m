function [freq,output] = FxPlotFFT(input,fs,option)
N=length(input); %get the number of points
k=0:N-1;     %create a vector from 0 to N-1
T=N/fs;      %get the frequency interval
freq=k/T;    %create the frequency range
X=fft(input)/N*2; % normalize the data

%only want the first half of the FFT, since it is redundant
cutOff = ceil(N);
% if nargin > 2
    cutOff = round(cutOff*0.5);
% end
%take only the first half of the spectrum
X = X(1:cutOff);
freq = freq(1:cutOff);

% figure;
% subplot(1,2,1);
% plot(T,input)
% subplot(1,2,2);

if nargin > 2
%     stem(freq,abs(X),'color',[0.85 0.33 0.1],'Markersize',4);
%     plot(freq,abs(X),'Markersize',4,'linewidth',2);
    switch option
        case 1
            plot(freq,abs(X),'color','k','Markersize',4,'linewidth',2);
        case 2
            semilogy(freq(2:end),abs(X(2:end)),'color','k','Markersize',4,'linewidth',2);
%             figure; plot(freq,log10(abs(X)),'color','k','Markersize',4,'linewidth',2);
    end
%     plot(freq,abs(X),'k','linewidth',2);
else
%     stem(freq,abs(X),'Markersize',4);
%     grid on
end

% xlabel('Freq (Hz)')
% ylabel('Amplitude')
% title('Using the positiveFFT function')


output = X;
